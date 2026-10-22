#!/usr/bin/env python3
"""Whether a pair-HMM likelihood prefers the truth where edit distance prefers the consensus.

usage: python3 -I experiments/msa/likelihood.py <karan> <set dir> <config> <jobs> [max letters]
       (PLEGMA_CONTEXT=1 for the homopolymer-aware model, see log_likelihood)

For every set whose consensus `<set>.fa.<config>.cons.fa` differs from its truth (and whose reads hold at most
`max letters` letters, 2,000,000 by default), prints the reads' summed log-likelihood under each and the summed edit
distance to each. The model sums over every alignment of a read to a template: before each template letter and
after the last, k inserted letters (each uniform) with probability pi^k (1 - pi); each template letter deleted with
probability pd, else read as itself with probability 1 - ps, as each other letter with ps / 3. Its rates are
estimated from the reads' alignments to the consensus. Each read's forward matrix is filled in a band of 32 letters
around its alignment to the template (karan, unit costs), in probability space, every column rescaled.
"""
import glob
import math
import os
import re
import subprocess
import sys
from multiprocessing import Pool

import numpy as np

BAND = 32
CONTEXT = os.environ.get("PLEGMA_CONTEXT") == "1"


def read_fasta(path):
    seqs = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    return ["".join(s) for s in seqs]


def alignments(karan, template_path, reads_path):
    out = subprocess.run([karan, template_path, reads_path], capture_output=True, text=True).stdout
    return [line.split("\t")[9] for line in out.strip().splitlines()[1:]]


def context_counts(cigar, template, read, tally):
    """Adds to `tally` the deletions and letters in runs and alone, and inserted letters next to a neighbour's."""
    i = j = 0
    for n, op in re.findall(r"(\d+)([=XID])", cigar):
        for _ in range(int(n)):
            if op == "I":
                near = {template[x] for x in (i - 1, i) if 0 <= x < len(template)}
                tally["ins_near"] += read[j] in near
                tally["ins"] += 1
                j += 1
                continue
            run = in_run(template, i)
            tally["run" if run else "alone"] += 1
            if op == "D":
                tally["del_run" if run else "del_alone"] += 1
            else:
                j += 1
            i += 1


def path_of(cigar, length):
    """Per template boundary i, the read offsets the alignment enters and leaves it at; and the edit counts."""
    enter = [0] * (length + 1)
    leave = [0] * (length + 1)
    counts = {"=": 0, "X": 0, "I": 0, "D": 0}
    i = j = 0
    for n, op in re.findall(r"(\d+)([=XID])", cigar):
        n = int(n)
        counts[op] += n
        for _ in range(n):
            if op == "I":
                j += 1
                continue
            leave[i] = j
            if op != "D":
                j += 1
            i += 1
            enter[i] = j
    leave[length] = j
    return enter, leave, counts


def in_run(template, i):
    return (i > 0 and template[i - 1] == template[i]) or (i + 1 < len(template) and template[i + 1] == template[i])


def log_likelihood(template, read, enter, leave, pi, pd, ps, context=None):
    """`context`, when given, is (pd in a run, pd alone, q): letters in a run of two or more are deleted at the
    first rate, others at the second, and an inserted letter is one of the slot's neighbours' letters with
    probability q (shared among them), any other letter else."""
    n = len(read)
    codes = np.frombuffer(read.encode(), dtype=np.uint8)
    lo = [max(0, enter[i] - BAND) for i in range(len(template) + 1)]
    hi = [min(n, leave[i] + BAND) for i in range(len(template) + 1)]
    log_scale = 0.0
    f = np.zeros(hi[0] - lo[0] + 1)
    f[0 - lo[0]] = 1.0
    stay = pi / 4.0
    for i in range(len(template) + 1):
        # The slot before letter i (after the last at i = len): inserted letters, then 1 - pi to leave it.
        h = f.copy()
        if context is None:
            for k in range(1, len(h)):
                h[k] += h[k - 1] * stay
        else:
            near = {template[x] for x in (i - 1, i) if 0 <= x < len(template)}
            q = context[2]
            for k in range(1, len(h)):
                letter = read[lo[i] + k - 1]
                weight = q / len(near) if letter in near else (1 - q) / (4 - len(near))
                h[k] += h[k - 1] * pi * weight
        g = h * (1.0 - pi)
        if i == len(template):
            break
        t = ord(template[i])
        deleted = pd if context is None else (context[0] if in_run(template, i) else context[1])
        nlo, nhi = lo[i + 1], hi[i + 1]
        nxt = np.zeros(nhi - nlo + 1)
        # Deleted: same read offset.
        a, b = max(lo[i], nlo), min(hi[i], nhi)
        if a <= b:
            nxt[a - nlo:b - nlo + 1] += g[a - lo[i]:b - lo[i] + 1] * deleted
        # Read as read[j - 1].
        a, b = max(lo[i] + 1, nlo, 1), min(hi[i] + 1, nhi)
        if a <= b:
            emitted = np.where(codes[a - 1:b] == t, 1.0 - ps, ps / 3.0)
            nxt[a - nlo:b - nlo + 1] += g[a - 1 - lo[i]:b - lo[i]] * (1.0 - deleted) * emitted
        top = nxt.max()
        if top <= 0:
            return -math.inf
        nxt /= top
        log_scale += math.log(top)
        f = nxt
    end = n - lo[len(template)]
    if end < 0 or end >= len(g) or g[end] <= 0:
        return -math.inf
    return log_scale + math.log(g[end])


def check(task):
    karan, fasta, config = task
    truth_path = fasta[:-3] + ".truth.fa"
    cons_path = fasta + "." + config + ".cons.fa"
    truth, cons = read_fasta(truth_path)[0], read_fasta(cons_path)[0]
    if truth == cons:
        return None
    reads = read_fasta(fasta)
    # Rates from the reads' alignments to the consensus.
    totals = {"=": 0, "X": 0, "I": 0, "D": 0}
    tally = {"ins_near": 0, "ins": 0, "run": 0, "alone": 0, "del_run": 0, "del_alone": 0}
    paths_cons = []
    for cigar, read in zip(alignments(karan, cons_path, fasta), reads):
        context_counts(cigar, cons, read, tally)
        enter, leave, counts = path_of(cigar, len(cons))
        paths_cons.append((enter, leave))
        for op in totals:
            totals[op] += counts[op]
    letters = len(cons) * len(reads)
    pd = totals["D"] / letters
    ps = totals["X"] / max(1, totals["="] + totals["X"])
    per_slot = totals["I"] / ((len(cons) + 1) * len(reads))
    pi = per_slot / (1 + per_slot)
    context = (tally["del_run"] / max(1, tally["run"]), tally["del_alone"] / max(1, tally["alone"]),
               tally["ins_near"] / max(1, tally["ins"])) if CONTEXT else None
    paths_truth = [path_of(c, len(truth))[:2] for c in alignments(karan, truth_path, fasta)]
    ll_cons = sum(log_likelihood(cons, r, e, l, pi, pd, ps, context) for r, (e, l) in zip(reads, paths_cons))
    ll_truth = sum(log_likelihood(truth, r, e, l, pi, pd, ps, context) for r, (e, l) in zip(reads, paths_truth))
    return os.path.basename(fasta), pi, pd, ps, ll_cons, ll_truth, context


def main():
    karan, folder, config, jobs = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    max_letters = int(sys.argv[5]) if len(sys.argv) > 5 else 2_000_000
    sets = []
    for f in sorted(glob.glob(os.path.join(folder, "*.fa"))):
        if f.endswith(".truth.fa") or ".cons." in f or not os.path.exists(f + "." + config + ".cons.fa"):
            continue
        if os.path.getsize(f) <= max_letters * 1.05:
            sets.append(f)
    with Pool(jobs) as pool:
        rows = [r for r in pool.map(check, [(karan, f, config) for f in sets]) if r]
    print("| set | pi | pd | ps | context (pd run, pd alone, q) | log L consensus | log L truth | prefers |")
    print("| :-- | --: | --: | --: | :-- | --: | --: | :-- |")
    better = 0
    for name, pi, pd, ps, lc, lt, context in rows:
        better += lt > lc
        shown = "%.4f, %.4f, %.3f" % context if context else "-"
        print("| %s | %.4f | %.4f | %.4f | %s | %.2f | %.2f | %s |" % (name, pi, pd, ps, shown, lc, lt,
                                                                        "truth" if lt > lc else "consensus"))
    print("\n%d of %d sets: the likelihood prefers the truth" % (better, len(rows)))


if __name__ == "__main__":
    main()
