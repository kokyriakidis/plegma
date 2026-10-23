#!/usr/bin/env python3
"""Sum-of-pairs costs of MSAs over their distinct rows, under the three costs that matter here.

usage: python3 -I experiments/vntreval/sp_score.py <msa.fa> [<msa.fa> ...]

- unit: the columns where two rows differ (both gaps cost nothing): vntr-eval's all_cost, as the MSA implies it;
- refine: unit plus 1 for every gap the pair opens: what plegma's refinement lowers (-I's gap 1);
- affine: 4 per mismatch, 6 + 2k per gap of k letters: vntr-eval's affine cost, as the MSA implies it.
Every pair of distinct rows is counted once; prints one line per MSA, and the columns.
"""
import sys


def read_rows(path):
    names, rows = [], []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            rows.append([])
        elif line:
            rows[-1].append(line.upper().replace(".", "-"))
    return names, ["".join(r) for r in rows]


def distinct(rows):
    seen, out = set(), []
    for r in rows:
        s = r.replace("-", "")
        if s not in seen:
            seen.add(s)
            out.append(r)
    keep = [c for c in range(len(out[0])) if any(r[c] != "-" for r in out)]
    return ["".join(r[c] for c in keep) for r in out]


def pair(x, y):
    """(differing columns, mismatches, gap letters, gap opens) of rows x and y."""
    diff = mis = gaps = opens = 0
    state = 0
    for a, b in zip(x, y):
        ga, gb = a == "-", b == "-"
        if ga and gb:
            continue
        if a != b:
            diff += 1
            if not ga and not gb:
                mis += 1
            else:
                gaps += 1
        now = 1 if ga else (2 if gb else 0)
        if now and now != state:
            opens += 1
        state = now
    return diff, mis, gaps, opens


def main():
    print("%-60s %6s %8s %10s %10s %10s %8s" % ("msa", "rows", "columns", "unit", "refine", "affine", "opens"))
    for path in sys.argv[1:]:
        rows = distinct(read_rows(path)[1])
        unit = refine = affine = opens = 0
        for a in range(len(rows)):
            for b in range(a + 1, len(rows)):
                d, m, g, o = pair(rows[a], rows[b])
                unit += d
                refine += d + o
                affine += 4 * m + 2 * g + 6 * o
                opens += o
        print("%-60s %6d %8d %10d %10d %10d %8d" % (path[-60:], len(rows), len(rows[0]), unit, refine, affine, opens))


if __name__ == "__main__":
    main()
