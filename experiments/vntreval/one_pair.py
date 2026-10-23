#!/usr/bin/env python3
"""One haplotype pair, the alignment an MSA implies for it against its optimal one, along the sequences.

usage: python3 -I experiments/vntreval/one_pair.py <karan> <msa.fa> <name a> <name b> [step] [period]

Prints, every `step` letters of a (200 by default), the letter of b each alignment pairs it with (or the last one
before it) and the diagonal (b's offset minus a's): where the two differ by a multiple of the repeat's period, the
MSA pairs different copies than the optimum does. Then the stretches of a where the implied alignment's cost runs
ahead of the optimal one, with their excess, and both costs.
"""
import re
import subprocess
import sys
import tempfile


def read_rows(path):
    names, rows = [], []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            rows.append([])
        elif line:
            rows[-1].append(line.upper())
    return dict(zip(names, ["".join(r) for r in rows]))


def main():
    karan, msa, na, nb = sys.argv[1:5]
    step = int(sys.argv[5]) if len(sys.argv) > 5 else 200
    period = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    rows = read_rows(msa)
    ra, rb = rows[na], rows[nb]
    a, b = ra.replace("-", ""), rb.replace("-", "")
    # The implied pairing and per-letter-of-a cost.
    implied_pair = [-1] * len(a)
    implied_cost = [0] * (len(a) + 1)
    i = j = 0
    cost = 0
    for x, y in zip(ra, rb):
        if x == "-" and y == "-":
            continue
        if x != y:
            cost += 1
        if x != "-":
            implied_pair[i] = j if y != "-" else -1
            i += 1
            implied_cost[i] = cost
        if y != "-":
            j += 1
    implied_cost[len(a)] = cost
    with tempfile.TemporaryDirectory() as tmp:
        fa = tmp + "/a.fa"
        fb = tmp + "/b.fa"
        open(fa, "w").write(">a\n%s\n" % a)
        open(fb, "w").write(">b\n%s\n" % b)
        line = subprocess.run([karan, fa, fb], capture_output=True, text=True).stdout.strip().splitlines()[-1]
    cigar = line.split("\t")[9]
    optimal_pair = [-1] * len(a)
    optimal_cost = [0] * (len(a) + 1)
    i = j = 0
    cost = 0
    for k, op in re.findall(r"(\d+)([=XID])", cigar):
        for _ in range(int(k)):
            if op == "D":
                cost += 1
                i += 1
                optimal_cost[i] = cost
            elif op == "I":
                cost += 1
                j += 1
            else:
                if op == "X":
                    cost += 1
                optimal_pair[i] = j
                i += 1
                j += 1
                optimal_cost[i] = cost
    optimal_cost[len(a)] = cost
    print("a %s %d letters, b %s %d letters; optimal cost %d, implied %d" % (
        na, len(a), nb, len(b), optimal_cost[-1], implied_cost[-1]))
    print("%8s %10s %10s %10s %10s %8s %8s" % ("a pos", "opt b", "opt diag", "impl b", "impl diag", "shift",
                                               "excess"))

    def last(pair, p):
        while p >= 0 and pair[p] < 0:
            p -= 1
        return (pair[p], p) if p >= 0 else (-1, p)

    for p in range(0, len(a), step):
        ob, op_ = last(optimal_pair, p)
        ib, ip = last(implied_pair, p)
        od = ob - op_ if ob >= 0 else None
        idg = ib - ip if ib >= 0 else None
        shift = (idg - od) if od is not None and idg is not None else None
        units = " (%+.2f units)" % (shift / period) if period and shift else ""
        print("%8d %10s %10s %10s %10s %8s%s %8d" % (p, ob, od, ib, idg, shift, units,
                                                     implied_cost[p] - optimal_cost[p]))


if __name__ == "__main__":
    main()
