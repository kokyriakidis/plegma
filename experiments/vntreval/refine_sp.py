#!/usr/bin/env python3
"""A prototype: an MSA refined on its sum of pairs, each distinct sequence realigned in turn to the others' profile.

usage: python3 -I experiments/vntreval/refine_sp.py <in.msa.fa> <out.msa.fa> [passes] [band]

The sum of pairs (unit costs): over every pair of distinct sequences, the columns where their rows differ, a column
where both are gaps costing nothing; vntr-eval's all_cost_over_opt divides it by the pairs' optima. With the other
rows fixed, a sequence's best rows for it is a DP over the columns: its letter against column c costs the others
that do not hold that letter there (gaps included), its gap there costs the others holding a letter, and its
letter in a column of its own costs every other. The DP keeps within `band` columns of the sequence's present
place (64 by default), so the present alignment is always one of its choices: a pass never raises the sum.
Identical sequences are refined once and share a row.
"""
import sys


def read(path):
    names, rows = [], []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            rows.append([])
        elif line:
            rows[-1].append(line.upper())
    return names, ["".join(r) for r in rows]


def sum_of_pairs(rows):
    n = len(rows)
    total = 0
    for c in range(len(rows[0])):
        counts = {}
        for r in rows:
            counts[r[c]] = counts.get(r[c], 0) + 1
        same = sum(k * (k - 1) // 2 for k in counts.values())
        total += n * (n - 1) // 2 - same
    return total


def realign(rows, s, band):
    """Row s realigned against the others; returns the new rows (columns all gaps dropped)."""
    others = [r for k, r in enumerate(rows) if k != s]
    m = len(others)
    columns = len(rows[0])
    letters = rows[s].replace("-", "")
    counts = []
    for c in range(columns):
        col = {}
        for r in others:
            col[r[c]] = col.get(r[c], 0) + 1
        counts.append(col)
    # The present column of each of the sequence's letters, for the band.
    place = [c for c, x in enumerate(rows[s]) if x != "-"]
    n = len(letters)
    INF = float("inf")
    # D[c][j]: the best cost of columns ..c-1 against letters ..j-1; band on j around the letters near column c.
    lo = [0] * (columns + 1)
    hi = [0] * (columns + 1)
    import bisect
    for c in range(columns + 1):
        k = bisect.bisect_left(place, c)  # letters placed before column c
        lo[c] = max(0, k - band)
        hi[c] = min(n, k + band)
    D = [dict() for _ in range(columns + 1)]
    back = [dict() for _ in range(columns + 1)]
    for j in range(lo[0], hi[0] + 1):
        if j == 0:
            D[0][0] = 0
        elif (j - 1) in D[0]:
            D[0][j] = D[0][j - 1] + m
            back[0][j] = "I"
    for c in range(1, columns + 1):
        col = counts[c - 1]
        gap_cost = m - col.get("-", 0)
        for j in range(lo[c], hi[c] + 1):
            best, how = INF, None
            v = D[c - 1].get(j)
            if v is not None and v + gap_cost < best:
                best, how = v + gap_cost, "D"
            if j > 0:
                v = D[c - 1].get(j - 1)
                if v is not None:
                    cost = v + m - col.get(letters[j - 1], 0)
                    if cost < best:
                        best, how = cost, "M"
                v = D[c].get(j - 1)
                if v is not None and v + m < best:
                    best, how = v + m, "I"
            if how:
                D[c][j] = best
                back[c][j] = how
    # Trace back: per column, the sequence's letter or gap, and its own columns inserted.
    c, j = columns, n
    new_row = []  # entries: (column index or None for an inserted column, symbol)
    while c > 0 or j > 0:
        how = back[c].get(j)
        if how == "M":
            new_row.append((c - 1, letters[j - 1]))
            c, j = c - 1, j - 1
        elif how == "D":
            new_row.append((c - 1, "-"))
            c -= 1
        else:
            new_row.append((None, letters[j - 1]))
            j -= 1
    new_row.reverse()
    out = [[] for _ in rows]
    for col, x in new_row:
        for k in range(len(rows)):
            if k == s:
                out[k].append(x)
            else:
                out[k].append(rows[k][col] if col is not None else "-")
    out = ["".join(r) for r in out]
    keep = [c for c in range(len(out[0])) if any(r[c] != "-" for r in out)]
    return ["".join(r[c] for c in keep) for r in out], D[columns][n]


def main():
    src, dst = sys.argv[1], sys.argv[2]
    passes = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    band = int(sys.argv[4]) if len(sys.argv) > 4 else 64
    names, rows = read(src)
    seqs = [r.replace("-", "") for r in rows]
    first = {}
    for k, s in enumerate(seqs):
        first.setdefault(s, k)
    distinct = sorted(first.values())
    work = [rows[k] for k in distinct]
    # Drop the columns gapped in every distinct row.
    keep = [c for c in range(len(work[0])) if any(r[c] != "-" for r in work)]
    work = ["".join(r[c] for c in keep) for r in work]
    print("start", sum_of_pairs(work))
    for p in range(passes):
        for s in range(len(work)):
            work, _ = realign(work, s, band)
        print("pass", p + 1, sum_of_pairs(work), flush=True)
    row_of = {seqs[k]: work[i] for i, k in enumerate(distinct)}
    with open(dst, "w") as f:
        for name, s in zip(names, seqs):
            f.write(">%s\n%s\n" % (name, row_of[s]))


if __name__ == "__main__":
    main()
