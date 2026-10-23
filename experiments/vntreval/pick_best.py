#!/usr/bin/env python3
"""Per region, the candidate MSA whose refinement objective is lowest, among several methods' (several starts).

usage: python3 -I experiments/vntreval/pick_best.py <vntr-eval> <out method> <gap open> <method> ...

The objective is the one plegma's refinement (-I) lowers: over every pair of distinct sequences, the columns where
their rows differ (a column where both are gaps costing nothing), plus `gap open` for every gap the pair opens. For
each region, the method whose work/plegma/<method>/<id>.msa.fa scores lowest gives candidates/<out method>/<id>.gfa
(a copy of its graph). Prints how often each method was picked.
"""
import collections
import glob
import os
import shutil
import sys


def read_rows(path):
    rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            rows.append([])
        elif line:
            rows[-1].append(line.upper())
    return ["".join(r) for r in rows]


def objective(rows, gap_open):
    distinct = list({r.replace("-", ""): r for r in rows}.values())
    total = 0
    for a in range(len(distinct)):
        ra = distinct[a]
        for b in range(a + 1, len(distinct)):
            rb = distinct[b]
            state = 0
            for x, y in zip(ra, rb):
                gx, gy = x == "-", y == "-"
                if gx and gy:
                    continue
                if x != y:
                    total += 1
                now = 1 if gx else (2 if gy else 0)
                if now and now != state:
                    total += gap_open
                state = now
    return total


def main():
    vntr_eval, out, gap_open = sys.argv[1], sys.argv[2], int(sys.argv[3])
    methods = sys.argv[4:]
    os.makedirs(os.path.join(vntr_eval, "candidates", out), exist_ok=True)
    picked = collections.Counter()
    regions = sorted(os.path.basename(p)[:-7] for p in glob.glob(
        os.path.join(vntr_eval, "work", "plegma", methods[0], "*.msa.fa")))
    for rid in regions:
        best = None
        for m in methods:
            msa = os.path.join(vntr_eval, "work", "plegma", m, rid + ".msa.fa")
            gfa = os.path.join(vntr_eval, "candidates", m, rid + ".gfa")
            if not (os.path.exists(msa) and os.path.exists(gfa)):
                continue
            score = objective(read_rows(msa), gap_open)
            if best is None or score < best[0]:
                best = (score, m, gfa)
        if best:
            shutil.copyfile(best[2], os.path.join(vntr_eval, "candidates", out, rid + ".gfa"))
            picked[best[1]] += 1
    for m in methods:
        print("%s picked for %d regions" % (m, picked[m]))


if __name__ == "__main__":
    main()
