#!/usr/bin/env python3
"""An MSA's cost as vntr-eval scores it: each pair aligned optimally between the letters the MSA stacks.

usage: python3 -I experiments/vntreval/implied.py <vntr-eval> <msa.fa> [<msa.fa> ...] [--pairs n]

vntr-eval builds a graph from an MSA (msa_graph.py: in each column the rows with one letter share a node) and
scores a pair of paths by the nodes they share: between consecutive shared nodes the two pieces are aligned
optimally (unit costs, or its affine 4/6/2), so only the identical letters an MSA stacks constrain a pair, each one
an anchor every optimal alignment of the pair must pass through. Per MSA, over its pairs of distinct rows: the
MSA's own induced unit cost (the columns where the two differ), that anchored cost, and the pairs' optima, with the
excess the anchors cost (anchored minus optimal) split into the pairs that carry it. Like vntr-eval, every pair of
input rows counts, so a pair of distinct sequences weighs the product of how many rows spell each. With `--pairs n`, the n pairs
whose anchors cost the most, and for the worst the anchored segment pairs that cost more than their optimum would
have alone (the whole pair is optimal only if the anchors lie on one of its optimal alignments).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp_score import distinct, read_rows  # noqa: E402


def lib(ve):
    sys.path.insert(0, os.path.join(ve, "tools"))
    import evaluate  # noqa: E402
    return evaluate.clib()


def anchored(ra, rb, ed):
    """(induced unit cost, anchored unit cost, segments [(a piece, b piece)]) of rows ra, rb."""
    induced = 0
    a_piece, b_piece = [], []
    segs = []
    cost = 0
    for x, y in zip(ra, rb):
        if x == "-" and y == "-":
            continue
        if x != y:
            induced += 1
        if x == y:
            # A shared node: it closes the open segment.
            if a_piece or b_piece:
                sa, sb = "".join(a_piece).encode(), "".join(b_piece).encode()
                cost += ed(sa, sb)
                segs.append((sa, sb))
                a_piece, b_piece = [], []
            continue
        if x != "-":
            a_piece.append(x)
        if y != "-":
            b_piece.append(y)
    if a_piece or b_piece:
        sa, sb = "".join(a_piece).encode(), "".join(b_piece).encode()
        cost += ed(sa, sb)
        segs.append((sa, sb))
    return induced, cost, segs


def main():
    args = sys.argv[1:]
    show = 0
    if "--pairs" in args:
        k = args.index("--pairs")
        show = int(args[k + 1])
        del args[k:k + 2]
    ve, paths = args[0], args[1:]
    C = lib(ve)
    for path in paths:
        every = read_rows(path)[1]
        rows = distinct(every)
        seqs = [r.replace("-", "").encode() for r in rows]
        count = {}
        for r in every:
            key = r.replace("-", "").upper().encode()
            count[key] = count.get(key, 0) + 1
        weight = [count[s] for s in seqs]
        tot_induced = tot_anchored = tot_opt = 0
        per = []
        for a in range(len(rows)):
            for b in range(a + 1, len(rows)):
                induced, cost, segs = anchored(rows[a], rows[b], C.ed)
                opt = C.ed(seqs[a], seqs[b])
                w = weight[a] * weight[b]
                tot_induced += w * induced
                tot_anchored += w * cost
                tot_opt += w * opt
                per.append((w * (cost - opt), a, b, induced, cost, opt))
        print("%s: %d rows; induced %d, anchored %d, optimal %d: anchored/optimal %.4f, induced/optimal %.4f" % (
            path, len(rows), tot_induced, tot_anchored, tot_opt, tot_anchored / max(1, tot_opt),
            tot_induced / max(1, tot_opt)))
        if show:
            per.sort(reverse=True)
            excess = sum(p[0] for p in per)
            print("  excess %d; the %d pairs costing most carry %d" % (excess, show, sum(p[0] for p in per[:show])))
            for ex, a, b, induced, cost, opt in per[:show]:
                print("    rows %2d-%2d (weights %d x %d): weighted excess %5d (anchored %d, optimal %d, induced %d)" % (
                    a, b, weight[a], weight[b], ex, cost, opt, induced))


if __name__ == "__main__":
    main()
