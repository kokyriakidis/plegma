#!/usr/bin/env python3
"""Sorts a consensus's errors against its truth into kinds.

usage: python3 -I experiments/msa/errors.py <karan> <truth.fa> <consensus.fa> [k] [distance]

The consensus is aligned to the truth with karan (unit costs) and each edit classified at its truth position:
- `substitution`;
- `homopolymer`: an indel in or beside a run of three or more of one letter;
- `tandem`: any other indel inside a tandem repeat (a k-mer through the position occurs again in the truth
  within `distance` letters; k 12, distance 500 by default);
- `other`: any other indel.
Prints `total substitution homopolymer tandem other` on one line.
"""
import re
import subprocess
import sys


def read_one(path):
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip().upper())
    return "".join(seq)


def run_length(s, i):
    """The length of the run of one letter covering position i (0 when i is outside s)."""
    if i < 0 or i >= len(s):
        return 0
    lo = i
    while lo > 0 and s[lo - 1] == s[i]:
        lo -= 1
    hi = i
    while hi + 1 < len(s) and s[hi + 1] == s[i]:
        hi += 1
    return hi - lo + 1


def in_tandem(s, i, k, distance, positions):
    for start in range(max(0, i - k + 1), min(i, len(s) - k) + 1):
        kmer = s[start:start + k]
        for p in positions.get(kmer, ()):
            if p != start and abs(p - start) <= distance:
                return True
    return False


def main():
    karan, truth_path, cons_path = sys.argv[1], sys.argv[2], sys.argv[3]
    k = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    distance = int(sys.argv[5]) if len(sys.argv) > 5 else 500
    truth = read_one(truth_path)
    out = subprocess.run([karan, truth_path, cons_path], capture_output=True, text=True).stdout.strip().splitlines()
    fields = out[-1].split("\t")
    cigar = fields[9]
    positions = {}
    for p in range(len(truth) - k + 1):
        positions.setdefault(truth[p:p + k], []).append(p)
    counts = {"substitution": 0, "homopolymer": 0, "tandem": 0, "other": 0}
    i = 0
    for n, op in re.findall(r"(\d+)([=XIDM])", cigar):
        n = int(n)
        if op in "=M":
            i += n
            continue
        for t in range(n):
            if op == "X":
                counts["substitution"] += 1
                i += 1
                continue
            at = i if op == "D" else i  # an insertion sits before truth letter i
            hp = max(run_length(truth, at), run_length(truth, at - 1)) >= 3
            if hp:
                counts["homopolymer"] += 1
            elif in_tandem(truth, min(at, len(truth) - 1), k, distance, positions):
                counts["tandem"] += 1
            else:
                counts["other"] += 1
            if op == "D":
                i += 1
    total = sum(counts.values())
    print(total, counts["substitution"], counts["homopolymer"], counts["tandem"], counts["other"])


if __name__ == "__main__":
    main()
