#!/usr/bin/env python3
"""Tracing plegma's refinement: whether its DP's costs are the objective's, and whether better MSAs exist near it.

usage: python3 -I experiments/vntreval/trace_refine.py exact <msa.fa> [gap] [sequences] [band]
       python3 -I experiments/vntreval/trace_refine.py force <msa.fa> <hap32.fa> <name a> <name b> [gap] [k]

The objective (what refine lowers): over every pair of distinct rows, the columns where the two differ (both gaps
cost nothing), plus `gap` (1) for every gap the pair opens.

`exact`: for the first `sequences` distinct rows (all by default), the refinement's three-state DP (a port of
plegma/refine.mojo's _realign, unbanded) realigns the row against the others; prints the DP's own cost of the
present row and of its best one, and the true objective before and after the new row is put in. Where the DP's
gain differs from the objective's, the DP's gap-open counts are not the objective's.

`force`: the seeds k-mers unique within both rows' array pieces (as one_region.py finds them) fix b's letters to
a's columns, and b is realigned by the same DP through them; prints the true objective before and after, and the
pair's own cost, so whether the seeds' placement is better for the whole MSA or only for this pair.
"""
import bisect
import sys

INF = float("inf")


def read_rows(path):
    names, rows = [], []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            rows.append([])
        elif line:
            rows[-1].append(line.upper())
    return names, ["".join(r) for r in rows]


def distinct_rows(names, rows):
    seen, out_names, out = set(), [], []
    for n, r in zip(names, rows):
        s = r.replace("-", "")
        if s not in seen:
            seen.add(s)
            out_names.append(n)
            out.append(r)
    keep = [c for c in range(len(out[0])) if any(r[c] != "-" for r in out)]
    return out_names, ["".join(r[c] for c in keep) for r in out]


def pair_cost(x, y, gap):
    cost, state = 0, 0
    for a, b in zip(x, y):
        ga, gb = a == "-", b == "-"
        if ga and gb:
            continue
        if a != b:
            cost += 1
        now = 1 if ga else (2 if gb else 0)
        if now and now != state:
            cost += gap
        state = now
    return cost


def objective(rows, gap):
    return sum(pair_cost(rows[a], rows[b], gap) for a in range(len(rows)) for b in range(a + 1, len(rows)))


def realign(rows, s, gap, forced=None, band=None):
    """The refinement's DP for row s against the others: its cost for the present row, its best cost and the new
    rows. `forced`: {letter index: column} the new row must keep. `band`: letters either side of the centre (the
    forced letters' places, interpolated, when forced; else the present place), unbanded when None."""
    others = [r for k, r in enumerate(rows) if k != s]
    m = len(others)
    C = len(rows[0])
    counts = [dict() for _ in range(C)]
    held = [0] * C
    starts = [0] * C
    ends = [0] * C
    for r in others:
        was_letter = True
        for c, x in enumerate(r):
            counts[c][x] = counts[c].get(x, 0) + 1
            if x != "-":
                held[c] += 1
                if not was_letter:
                    ends[c] += 1
            elif was_letter:
                starts[c] += 1
            was_letter = x != "-"
    letters = rows[s].replace("-", "")
    n = len(letters)
    # Present cost under the DP's model.
    present, in_gap = 0, False
    for c, x in enumerate(rows[s]):
        if x != "-":
            present += m - counts[c].get(x, 0) + gap * starts[c]
            in_gap = False
        else:
            present += held[c] + (gap * ends[c] if in_gap else gap * held[c])
            in_gap = True
    g = gap
    # The band's centre per column: letters before it.
    before = [0] * (C + 1)
    k = 0
    for c in range(C):
        before[c] = k
        k += rows[s][c] != "-"
    before[C] = k
    centre = before
    if forced:
        points = sorted((col, j) for j, col in forced.items())
        pc = [0] + [p[0] for p in points] + [C]
        pj = [0] + [p[1] for p in points] + [n]
        centre = [0] * (C + 1)
        q = 0
        for c in range(C + 1):
            while q + 1 < len(pc) - 1 and pc[q + 1] <= c:
                q += 1
            c0, c1, j0, j1 = pc[q], pc[q + 1], pj[q], pj[q + 1]
            centre[c] = j0 if c1 <= c0 else j0 + (c - c0) * (j1 - j0) // (c1 - c0)
    def inside(c, j):
        return band is None or abs(j - centre[c]) <= band
    M = [[INF] * (n + 1) for _ in range(C + 1)]
    G = [[INF] * (n + 1) for _ in range(C + 1)]
    O = [[INF] * (n + 1) for _ in range(C + 1)]
    BM = [[0] * (n + 1) for _ in range(C + 1)]
    BG = [[0] * (n + 1) for _ in range(C + 1)]
    BO = [[0] * (n + 1) for _ in range(C + 1)]
    M[0][0] = 0
    for c in range(C + 1):
        x = c - 1
        held_here = held[x] if c > 0 else m
        lo = 0 if band is None else max(0, centre[c] - band)
        hi = n if band is None else min(n, centre[c] + band)
        for j in range(lo, hi + 1):
            if c > 0:
                if j > 0 and (forced is None or forced.get(j - 1, x) == x):
                    unit = m - counts[x].get(letters[j - 1], 0)
                    best, how = INF, 0
                    for v, h in ((M[c - 1][j - 1] + g * starts[x], 0), (G[c - 1][j - 1] + g * starts[x], 1),
                                 (O[c - 1][j - 1], 2)):
                        if v < best:
                            best, how = v, h
                    M[c][j], BM[c][j] = best + unit, how
                best, how = INF, 0
                for v, h in ((M[c - 1][j] + g * held[x], 0), (G[c - 1][j] + g * ends[x], 1),
                             (O[c - 1][j] + g * held[x], 2)):
                    if v < best:
                        best, how = v, h
                G[c][j], BG[c][j] = best + held[x], how
            if j > 0 and (forced is None or (j - 1) not in forced):
                best, how = INF, 0
                for v, h in ((M[c][j - 1] + g * held_here, 0), (G[c][j - 1] + g * held_here, 1), (O[c][j - 1], 2)):
                    if v < best:
                        best, how = v, h
                O[c][j], BO[c][j] = best + m, how
    finals = [(M[C][n], 0), (G[C][n], 1), (O[C][n], 2)]
    best_cost, state = min(finals)
    steps, mine = [], []
    c, j = C, n
    while c > 0 or j > 0:
        if state == 0:
            came = BM[c][j]
            steps.append(c - 1)
            mine.append(letters[j - 1])
            c, j = c - 1, j - 1
        elif state == 1:
            came = BG[c][j]
            steps.append(c - 1)
            mine.append("-")
            c -= 1
        else:
            came = BO[c][j]
            steps.append(-1)
            mine.append(letters[j - 1])
            j -= 1
        state = came
    steps.reverse()
    mine.reverse()
    new = []
    for k, r in enumerate(rows):
        if k == s:
            new.append("".join(mine))
        else:
            new.append("".join(r[o] if o >= 0 else "-" for o in steps))
    keep = [t for t in range(len(new[0])) if any(r[t] != "-" for r in new)]
    new = ["".join(r[t] for t in keep) for r in new]
    return present, best_cost, new


def seeds_in_window(a, b, k):
    """Chained k-mers unique within each whole sequence's array piece and shared (as one_region.py): (i, j)."""
    def kpos(s):
        p = {}
        for i in range(len(s) - k + 1):
            p.setdefault(s[i:i + k], []).append(i)
        return p
    pa, pb = kpos(a), kpos(b)
    seeds = sorted((ia[0], pb[w][0]) for w, ia in pa.items() if len(ia) == 1 and len(pb.get(w, ())) == 1)
    idx, back = [], [-1] * len(seeds)
    js = []
    for t, (_, j) in enumerate(seeds):
        q = bisect.bisect_left(js, j)
        if q > 0:
            back[t] = idx[q - 1]
        if q == len(idx):
            idx.append(t)
            js.append(j)
        else:
            idx[q] = t
            js[q] = j
    out, t = [], idx[-1] if idx else -1
    while t >= 0:
        out.append(seeds[t])
        t = back[t]
    return out[::-1]


def main():
    mode = sys.argv[1]
    if mode == "exact":
        names, rows = distinct_rows(*read_rows(sys.argv[2]))
        gap = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        count = int(sys.argv[4]) if len(sys.argv) > 4 else len(rows)
        base = objective(rows, gap)
        print("objective %d over %d distinct rows, %d columns" % (base, len(rows), len(rows[0])))
        print("%-28s %10s %10s %10s %12s %12s" % ("row", "DP present", "DP best", "DP gain", "true before", "true gain"))
        band = int(sys.argv[5]) if len(sys.argv) > 5 else None
        for s in range(count):
            present, best, new = realign(rows, s, gap, band=band)
            after = objective(new, gap)
            print("%-28s %10d %10d %10d %12d %12d%s" % (names[s][:28], present, best, present - best, base,
                                                       base - after, "" if present - best == base - after else
                                                       "   <- DP gain differs"))
    elif mode == "force":
        msa, hap, na, nb = sys.argv[2:6]
        gap = int(sys.argv[6]) if len(sys.argv) > 6 else 1
        k = int(sys.argv[7]) if len(sys.argv) > 7 else 21
        names, rows = distinct_rows(*read_rows(msa))
        ia, ib = names.index(na), names.index(nb)
        a, b = rows[ia].replace("-", ""), rows[ib].replace("-", "")
        col_of_a = [c for c, x in enumerate(rows[ia]) if x != "-"]
        seeds = seeds_in_window(a, b, k)
        forced = {}
        for i, j in seeds:
            for t in range(k):
                forced[j + t] = col_of_a[i + t]
        base = objective(rows, gap)
        present, best, new = realign(rows, ib, gap, forced=forced, band=300)
        after = objective(new, gap)
        print("seeds %d (%d letters of b fixed); objective %d -> %d (%+d); DP present %d, forced best %d" % (
            len(seeds), len(forced), base, after, after - base, present, best))
        print("pair %s / %s: cost %d -> %d" % (na, nb, pair_cost(rows[ia], rows[ib], gap), pair_cost(new[ia], new[ib], gap)))
        p2, b2, free = realign(rows, ib, gap, band=300)
        print("free (unforced, band 300 around its place) realignment of b: objective %d -> %d; pair cost %d" % (
            base, objective(free, gap), pair_cost(free[ia], free[ib], gap)))


if __name__ == "__main__":
    main()
