#!/usr/bin/env python3
"""The consensus the reads support when each is placed against the truth: the errors no MSA can avoid.

usage: python3 -I experiments/msa/oracle.py <karan> <set dir> <jobs>

For every `<set>.fa` with a `<set>.truth.fa` in `set dir`, each read is aligned to the truth with karan (unit
costs) and the truth's columns vote: a truth letter is kept as the reads' majority letter, or dropped when most
reads delete it; between two truth letters the j-th inserted letter is kept when most reads insert at least j
there, as the majority letter among them. Writes `<set>.fa.oracle.cons.fa`, so accuracy.py and error_kinds.py
can score it as the configuration `oracle`.
"""
import collections
import glob
import os
import re
import subprocess
import sys
from multiprocessing import Pool


def read_fasta(path):
    names, seqs = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    return names, ["".join(s) for s in seqs]


def vote(fasta, karan):
    truth_path = fasta[:-3] + ".truth.fa"
    truth = read_fasta(truth_path)[1][0]
    _, reads = read_fasta(fasta)
    lines = subprocess.run([karan, truth_path, fasta], capture_output=True, text=True).stdout.strip().splitlines()
    # Per truth letter the letters reads give it ("" a deletion); per slot before it the reads' insertions.
    letters = [collections.Counter() for _ in truth]
    inserted = [[] for _ in range(len(truth) + 1)]
    for line, read in zip(lines[1:], reads):
        cigar = line.split("\t")[9]
        i = j = 0
        insertion = []
        for n, op in re.findall(r"(\d+)([=XIDM])", cigar):
            n = int(n)
            if op == "I":
                insertion.extend(read[j:j + n])
                j += n
                continue
            if insertion:
                inserted[i].append(insertion)
                insertion = []
            for _ in range(n):
                if op == "D":
                    letters[i][""] += 1
                    i += 1
                else:
                    letters[i][read[j]] += 1
                    i += 1
                    j += 1
        if insertion:
            inserted[i].append(insertion)
    half = len(reads) / 2
    out = []
    for i in range(len(truth) + 1):
        runs = inserted[i]
        k = 0
        while sum(1 for r in runs if len(r) > k) > half:
            out.append(collections.Counter(r[k] for r in runs if len(r) > k).most_common(1)[0][0])
            k += 1
        if i < len(truth):
            out.append(letters[i].most_common(1)[0][0])
    with open(fasta + ".oracle.cons.fa", "w") as f:
        f.write(">oracle\n%s\n" % "".join(out))


def main():
    karan, folder, jobs = sys.argv[1], sys.argv[2], int(sys.argv[3])
    sets = [f for f in sorted(glob.glob(os.path.join(folder, "*.fa")))
            if not f.endswith(".truth.fa") and ".cons." not in f and os.path.exists(f[:-3] + ".truth.fa")]
    with Pool(jobs) as pool:
        pool.starmap(vote, [(f, karan) for f in sets])


if __name__ == "__main__":
    main()
