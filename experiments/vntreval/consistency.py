#!/usr/bin/env python3
"""How consistent the haplotypes' optimal pairwise alignments are with each other: whether one MSA can hold them all.

usage: python3 -I experiments/vntreval/consistency.py <vntr-eval> <karan> <regions> [triples] [karan costs]

For each region (comma list), every pair of distinct hap32 sequences is aligned optimally (karan, unit costs,
indels placed left). For sampled triples (a, b, c) (500 by default), every letter of a paired with a letter of b
that is in turn paired with a letter of c: the share where a's alignment to c pairs the same two letters
(transitive). An MSA keeps every pairing only where all are transitive; where they are not, some pair's optimal
alignment must give way, which is the floor under the sum of pairs over the optima.
"""
import itertools
import os
import random
import re
import subprocess
import sys
import tempfile


def read_fasta(path):
    seqs = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    return ["".join(s) for s in seqs]


def pairing(cigar, n):
    """For each query letter, the reference letter it is paired with (=, X) or -1."""
    out = [-1] * n
    i = j = 0
    for k, op in re.findall(r"(\d+)([=XID])", cigar):
        for _ in range(int(k)):
            if op == "D":
                i += 1
            elif op == "I":
                j += 1
            else:
                out[j] = i
                i += 1
                j += 1
    return out


def main():
    vntr_eval, karan, regions = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
    triples = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    costs = sys.argv[5] if len(sys.argv) > 5 else "edit"
    rng = random.Random(1)
    for rid in regions:
        seqs = list(dict.fromkeys(read_fasta(os.path.join(vntr_eval, "regions", rid, "hap32.fa"))))
        n = len(seqs)
        pairs = list(itertools.permutations(range(n), 2))
        with tempfile.TemporaryDirectory() as tmp:
            tsv = os.path.join(tmp, "pairs.tsv")
            with open(tsv, "w") as f:
                for k, (a, b) in enumerate(pairs):
                    f.write("%d\t%s\t%s\n" % (k, seqs[b], seqs[a]))
            out = subprocess.run([karan, "--pairs", tsv, "--costs", costs], capture_output=True, text=True).stdout
        maps = {}
        for line in out.strip().splitlines()[1:]:
            f = line.split("\t")
            a, b = pairs[int(f[0])]
            maps[(a, b)] = pairing(f[9], len(seqs[a]))
        chosen = [rng.sample(range(n), 3) for _ in range(triples)]
        total = kept = 0
        for a, b, c in chosen:
            ab, bc, ac = maps[(a, b)], maps[(b, c)], maps[(a, c)]
            for i, j in enumerate(ab):
                if j < 0:
                    continue
                k = bc[j]
                if k < 0:
                    continue
                total += 1
                kept += ac[i] == k
        print("%s (%s): %d distinct haplotypes, %d triples: %.4f of chained pairings transitive" % (
            rid, costs, n, triples, kept / total if total else 1.0))


if __name__ == "__main__":
    main()
