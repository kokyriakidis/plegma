"""The co-optimal letter pairs of sequence pairs, for qc_loop.py --coopt: which letter pairs some optimal alignment
aligns, so an MSA placing an indel where another optimal alignment would is not counted wrong.

usage: coopt <pairs.tsv> [--two-piece x,o1,e1,o2,e2]
"""

from std.os import getenv
from std.sys import argv

from karan import Aligner, Costs, Mode, Ties

from plegma.graph import DELETION, INSERTION, cigar_runs

comptime FAR = Int32(1 << 29)


def _band_of(cigar: String, n: Int, m: Int, mut lo: List[Int], mut hi: List[Int]) raises:
    """Widen [lo[i], hi[i]] (rows 0..n of the DP over a's letters) to every column the alignment `cigar` of a (the
    reference) and b visits in row i."""
    var runs = cigar_runs(cigar)
    var i = 0
    var j = 0
    lo[0] = min(lo[0], 0)
    hi[0] = max(hi[0], 0)
    for x in range(len(runs.ops)):
        var op = runs.ops[x]
        for _ in range(runs.counts[x]):
            if op == DELETION:
                i += 1
            elif op == INSERTION:
                j += 1
            else:
                i += 1
                j += 1
            lo[i] = min(lo[i], j)
            hi[i] = max(hi[i], j)


def _coopt_set(args: List[String]) raises:
    """`coopt <pairs.tsv> [--two-piece x,o1,e1,o2,e2]`: per line `a b` (two sequences), under two-piece costs
    (mismatch 48, gaps min(32 + 24k, 192 + 10k) by default), every letter pair some optimal alignment aligns. The forward and backward DP run only between the
    pair's leftmost and rightmost optimal alignments (karan, Ties.LEFT and Ties.RIGHT) and MARGIN columns either
    side, MARGIN being 64 letters and the extremes' longest gap (an optimal alignment that splits it may shift by up
    to its length between them); the band is widened (the margin four times over) and the pair done again until its optimum is karan's and a
    band four times wider finds the same co-optimal letter pairs (with two-piece gaps some optimal alignments stray
    past the two extremes, and one that leaves a band is seen by neither pass inside it). COOPT_MARGIN sets
    the first margin (a margin past both lengths is the full DP). Prints per line `best n partners`, `partners` one token per letter of a: `.` for none, else
    its partners in b as intervals `j` or `j1-j2` joined by `:`."""
    var x = Int32(48)
    var o1 = Int32(32)
    var e1 = Int32(24)
    var o2 = Int32(192)
    var e2 = Int32(10)
    for k in range(len(args)):
        if args[k] == "--two-piece" and k + 1 < len(args):
            var parts = args[k + 1].split(",")
            x = Int32(Int(String(parts[0])))
            o1 = Int32(Int(String(parts[1])))
            e1 = Int32(Int(String(parts[2])))
            o2 = Int32(Int(String(parts[3])))
            e2 = Int32(Int(String(parts[4])))
    var costs = Costs.two_piece(Int(x), Int(o1), Int(e1), Int(o2), Int(e2))
    var text: String
    with open(args[1], "r") as f:
        text = f.read()
    var aligner = Aligner()
    for line in text.split("\n"):
        var fields = String(line).split("\t")
        if len(fields) < 2:
            continue
        var a = List[UInt8](String(fields[0]).as_bytes())
        var b = List[UInt8](String(fields[1]).as_bytes())
        var n = len(a)
        var m = len(b)
        var left = aligner.align(Span(a), Span(b), costs, Mode.GLOBAL, ties=Ties.LEFT, max_memory=4 << 30)
        var right = aligner.align(Span(a), Span(b), costs, Mode.GLOBAL, ties=Ties.RIGHT, max_memory=4 << 30)
        var optimum = Int32(left.cost)
        var path_lo = List[Int](length=n + 1, fill=m + 1)
        var path_hi = List[Int](length=n + 1, fill=-1)
        _band_of(left.cigar, n, m, path_lo, path_hi)
        _band_of(right.cigar, n, m, path_lo, path_hi)
        # An optimal alignment may split a gap of the extremes and shift by up to its length between them.
        var longest = 0
        for side in range(2):
            var runs = cigar_runs(left.cigar if side == 0 else right.cigar)
            for r in range(len(runs.ops)):
                if runs.ops[r] == DELETION or runs.ops[r] == INSERTION:
                    longest = max(longest, runs.counts[r])
        var margin = 64 + longest
        var last_cells = -1
        var forced = getenv("COOPT_MARGIN")
        if forced != "":
            margin = Int(forced)
        while True:
            # The band: the two paths' columns per row and `margin` either side, kept monotone so each row's band
            # reaches the next's.
            var lo = List[Int](length=n + 1, fill=0)
            var hi = List[Int](length=n + 1, fill=m)
            for i in range(n + 1):
                lo[i] = max(0, path_lo[i] - margin)
                hi[i] = min(m, path_hi[i] + margin)
            for i in range(1, n + 1):
                if lo[i] < lo[i - 1]:
                    lo[i] = lo[i - 1]
            for i in reversed(range(n)):
                if hi[i] > hi[i + 1]:
                    hi[i] = hi[i + 1]
            var start = List[Int](length=n + 2, fill=0)
            for i in range(n + 1):
                start[i + 1] = start[i] + (hi[i] - lo[i] + 1)
            # Forward, every row's H kept (the band's cells); F1 F2 (gaps in b, down a column) rolling.
            var fwd = List[Int32](length=start[n + 1], fill=FAR)
            var f1 = List[Int32](length=m + 1, fill=FAR)
            var f2 = List[Int32](length=m + 1, fill=FAR)
            for jj in range(lo[0], hi[0] + 1):
                fwd[start[0] + jj - lo[0]] = Int32(0) if jj == 0 else min(o1 + e1 * Int32(jj), o2 + e2 * Int32(jj))
            for ii in range(1, n + 1):
                var ca = a[ii - 1]
                var g1 = FAR
                var g2 = FAR
                var row = start[ii] - lo[ii]
                var prev = start[ii - 1] - lo[ii - 1]
                for jj in range(lo[ii], hi[ii] + 1):
                    var up = fwd[prev + jj] if jj >= lo[ii - 1] and jj <= hi[ii - 1] else FAR
                    var diag = fwd[prev + jj - 1] if jj - 1 >= lo[ii - 1] and jj - 1 <= hi[ii - 1] else FAR
                    f1[jj] = min(up + o1 + e1, f1[jj] + e1)
                    f2[jj] = min(up + o2 + e2, f2[jj] + e2)
                    var h = min(f1[jj], f2[jj])
                    if jj > 0:
                        var west = fwd[row + jj - 1] if jj - 1 >= lo[ii] else FAR
                        g1 = min(west + o1 + e1, g1 + e1)
                        g2 = min(west + o2 + e2, g2 + e2)
                        if diag < FAR:
                            h = min(h, diag + (Int32(0) if ca == b[jj - 1] else x))
                        h = min(h, min(g1, g2))
                    fwd[row + jj] = h
                # A column leaving the band upward resets its gap states (they are FAR above the band).
                if ii < n:
                    for jj in range(lo[ii], lo[ii + 1]):
                        f1[jj] = FAR
                        f2[jj] = FAR
            var best = fwd[start[n] + m - lo[n]] if lo[n] <= m and m <= hi[n] else FAR
            if best != optimum and margin < n + m:
                margin *= 4
                continue
            # Backward, rolling; each match cell (ii, jj) checked against the forward H at (ii - 1, jj - 1).
            var hb_next = List[Int32](length=m + 2, fill=FAR)
            var hb_row = List[Int32](length=m + 2, fill=FAR)
            var fb1 = List[Int32](length=m + 1, fill=FAR)
            var fb2 = List[Int32](length=m + 1, fill=FAR)
            var partners = List[List[Int]]()
            for _ in range(n):
                partners.append(List[Int]())
            for ii in reversed(range(n + 1)):
                var eb1 = FAR
                var eb2 = FAR
                for jj in reversed(range(lo[ii], hi[ii] + 1)):
                    var h: Int32
                    if ii == n and jj == m:
                        h = 0
                    else:
                        h = FAR
                        if ii < n and jj < m and jj + 1 >= lo[ii + 1] and jj + 1 <= hi[ii + 1]:
                            h = hb_next[jj + 1] + (Int32(0) if a[ii] == b[jj] else x)
                        if jj < m:
                            h = min(h, min(eb1 + o1 + e1, eb2 + o2 + e2))
                        if ii < n:
                            h = min(h, min(fb1[jj] + o1 + e1, fb2[jj] + o2 + e2))
                    hb_row[jj] = h
                    eb1 = min(h, eb1 + e1) if jj < m else h
                    eb2 = min(h, eb2 + e2) if jj < m else h
                    fb1[jj] = min(h, fb1[jj] + e1) if ii < n else h
                    fb2[jj] = min(h, fb2[jj] + e2) if ii < n else h
                    if ii >= 1 and jj >= 1 and jj - 1 >= lo[ii - 1] and jj - 1 <= hi[ii - 1]:
                        var f = fwd[start[ii - 1] + jj - 1 - lo[ii - 1]]
                        if f < FAR and h < FAR and f + (Int32(0) if a[ii - 1] == b[jj - 1] else x) + h == best:
                            partners[ii - 1].append(jj - 1)
                # Columns of this row outside the next row's band read as FAR there.
                if ii > 0:
                    for jj in range(hi[ii - 1] + 1, hi[ii] + 1):
                        fb1[jj] = FAR
                        fb2[jj] = FAR
                for jj in range(lo[ii], hi[ii] + 1):
                    hb_next[jj] = hb_row[jj]
                if ii > 0 and lo[ii] > lo[ii - 1]:
                    for jj in range(lo[ii - 1], lo[ii]):
                        hb_next[jj] = FAR
            # Stable across a band four times wider: some optimal alignments stray past the two extremes (a long gap
            # split in two, its pieces before and after theirs), and one that leaves the band is seen by neither
            # pass inside it; so the set is kept once a band four times wider finds the same cells.
            var cells = 0
            for i in range(n):
                cells += len(partners[i])
            if getenv("COOPT_STABLE") == "1" and margin < n + m and cells != last_cells:
                last_cells = cells
                margin *= 4
                continue
            var out = String(best, " ", n, " ")
            for i in range(n):
                if i > 0:
                    out += ","
                ref ps = partners[i]
                if len(ps) == 0:
                    out += "."
                    continue
                # Descending as the backward pass found them: intervals ascending.
                var k = len(ps) - 1
                var first = True
                while k >= 0:
                    var s0 = ps[k]
                    var s1 = s0
                    while k - 1 >= 0 and ps[k - 1] == s1 + 1:
                        k -= 1
                        s1 = ps[k]
                    if not first:
                        out += ":"
                    out += String(s0) if s0 == s1 else String(s0, "-", s1)
                    first = False
                    k -= 1
            print(out)
            break


def main() raises:
    var args = List[String]()
    for a in argv():
        args.append(String(a))
    if len(args) < 2:
        print("usage: coopt <pairs.tsv> [--two-piece x,o1,e1,o2,e2]")
        return
    _coopt_set(args)
