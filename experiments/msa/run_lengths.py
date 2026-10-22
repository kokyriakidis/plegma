#!/usr/bin/env python3
"""Whether calling each homopolymer run's length from the reads' run lengths mends a consensus's run errors.

usage: python3 -I experiments/msa/run_lengths.py <karan> <set dir> <config> <jobs> [rounds] [max bytes]

For every set with a truth and a `<set>.fa.<config>.cons.fa`: the reads are aligned to the consensus (karan, unit
costs, indels left), and each run of one letter in it, h letters long, is read in each read as the letters the
alignment puts against it, k of them when they are all that letter (a read with another letter there says nothing).
P(k - h | h), from every run of every read (h over 8 as 8, the difference within -4..4, each count one more), gives
each run the length h' near h (h - 2 ..= h + 2, at least 1) that makes the reads' run lengths likeliest. The
consensus so mended is aligned to again for `rounds` rounds (3 by default). Prints each set's errors before and
after (errors.py's kinds), and the totals.
"""
import collections
import glob
import math
import os
import re
import subprocess
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = 8
SPREAD = 4


def read_fasta(path):
    seqs = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    return ["".join(s) for s in seqs]


def runs_of(seq):
    out = []
    i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        out.append((i, j - i))
        i = j
    return out


def observed(karan, cons_path, cons, reads_path, reads):
    """Per run of the consensus, the run lengths the reads give it (None where a read holds another letter there):
    the read letters from where its alignment enters the run to where it enters the next letter, and those of the
    run's letter inserted right after it."""
    lines = subprocess.run([karan, cons_path, reads_path], capture_output=True, text=True).stdout.strip().splitlines()
    runs = runs_of(cons)
    seen = [[] for _ in runs]
    for line, read in zip(lines[1:], reads):
        cigar = line.split("\t")[9]
        enter = [0] * (len(cons) + 1)
        leave = [0] * (len(cons) + 1)
        i = j = 0
        for n, op in re.findall(r"(\d+)([=XID])", cigar):
            for _ in range(int(n)):
                if op == "I":
                    j += 1
                    continue
                leave[i] = j
                if op != "D":
                    j += 1
                i += 1
                enter[i] = j
        leave[len(cons)] = j
        for r, (start, h) in enumerate(runs):
            letter = cons[start]
            lo, hi = enter[start], enter[start + h]
            while hi < leave[start + h] and read[hi] == letter:
                hi += 1
            piece = read[lo:hi]
            if not piece:
                seen[r].append(0)
            elif piece.count(letter) == len(piece):
                seen[r].append(len(piece))
            else:
                seen[r].append(None)
    return runs, seen


def errors(karan, truth, cons):
    out = subprocess.run(["python3", "-I", os.path.join(HERE, "errors.py"), karan, truth, cons],
                         capture_output=True, text=True).stdout.split()
    return [int(x) for x in out] if len(out) == 5 else None


def mend(task):
    karan, fasta, config, rounds = task
    truth_path = fasta[:-3] + ".truth.fa"
    cons_path = fasta + "." + config + ".cons.fa"
    before = errors(karan, truth_path, cons_path)
    reads = read_fasta(fasta)
    cons = read_fasta(cons_path)[0]
    work = fasta + ".runs.cons.fa"
    open(work, "w").write(">c\n%s\n" % cons)
    for _ in range(rounds):
        runs, seen = observed(karan, work, cons, fasta, reads)
        table = collections.defaultdict(lambda: [1] * (2 * SPREAD + 1))
        for (start, h), ks in zip(runs, seen):
            for k in ks:
                if k is not None and abs(k - h) <= SPREAD:
                    table[min(h, TOP)][k - h + SPREAD] += 1

        def log_p(k, h):
            row = table[min(h, TOP)]
            d = k - h
            if abs(d) > SPREAD:
                return math.log(0.5 / sum(row))
            return math.log(row[d + SPREAD] / sum(row))

        out = []
        changed = 0
        for (start, h), ks in zip(runs, seen):
            ks = [k for k in ks if k is not None]
            best_h, best = h, sum(log_p(k, h) for k in ks)
            for other in range(max(1, h - 2), h + 3):
                score = sum(log_p(k, other) for k in ks)
                if score > best + 1e-9:
                    best_h, best = other, score
            changed += best_h != h
            out.append(cons[start] * best_h)
        cons = "".join(out)
        open(work, "w").write(">c\n%s\n" % cons)
        if changed == 0:
            break
    after = errors(karan, truth_path, work)
    return os.path.basename(fasta), before, after


def main():
    karan, folder, config, jobs = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    rounds = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    max_bytes = int(sys.argv[6]) if len(sys.argv) > 6 else 400_000
    sets = [f for f in sorted(glob.glob(os.path.join(folder, "*.fa")))
            if not f.endswith(".truth.fa") and ".cons." not in f and os.path.getsize(f) <= max_bytes
            and os.path.exists(f + "." + config + ".cons.fa")]
    with Pool(jobs) as pool:
        rows = pool.map(mend, [(karan, f, config, rounds) for f in sets])
    totals = [[0] * 5, [0] * 5]
    print("| set | errors before | after | homopolymer before | after |")
    print("| :-- | --: | --: | --: | --: |")
    for name, before, after in rows:
        if before is None or after is None:
            print("| %s | failed |" % name)
            continue
        for k in range(5):
            totals[0][k] += before[k]
            totals[1][k] += after[k]
        if before[0] or after[0]:
            print("| %s | %d | %d | %d | %d |" % (name, before[0], after[0], before[2], after[2]))
    print("\ntotal %d -> %d (substitution %d -> %d, homopolymer %d -> %d, tandem %d -> %d, other %d -> %d)" % (
        totals[0][0], totals[1][0], totals[0][1], totals[1][1], totals[0][2], totals[1][2], totals[0][3],
        totals[1][3], totals[0][4], totals[1][4]))


if __name__ == "__main__":
    main()
