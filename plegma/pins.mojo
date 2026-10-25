"""minimap2's pinned alignment of a sequence to the reference, so `on_reference` writes each sequence's differences
against the reference as minimap2 (2.31, -x asm5 as dipcall runs it: -z200000,10000 -r2k), and so dipcall and the
benchmarks dipcall makes, write them.

The exact global alignment writes a copy-number change inside a tandem array as one event; minimap2 aligns only
between pins (the middles of chained seeds both sequences advanced `MIN_KSW_LEN` letters or more to) and so writes
it at the pins' places, often in pieces. Truth sets made with dipcall hold those pieces, and record-by-record
comparisons (truvari bench) count every other writing of the same haplotype as errors.

Ported line by line from minimap2's sketch.c (mm_sketch), seed.c (mm_seed_mz_flt, mm_seed_collect_all,
mm_seed_select), index.c (mm_idx_cal_max_occ), lchain.c (mg_lchain_rmq, comput_sc_simple, mg_chain_backtrack),
map.c (the long-join re-chain) and align.c (mm_fix_bad_ends, mm_filter_bad_seeds, mm_filter_bad_seeds_alt and
mm_align1's gap filling); experiments/chr20/mm2like.py is the same port in Python, its seeds and chains checked
against minimap2's own (--print-seeds). The RMQ's ties fall to the later anchor (minimap2's to its tree's shape).
The flanks, CHM13's around the region, are context for the seeds and the chain only: the pins inside the sequences
are kept and the stretches between them aligned exactly (karan), so the rows spell the sequences whatever the flanks.
"""

from std.math import log2
from std.os import getenv
from std.sys import stderr

from karan import Aligner, Costs, Mode, Ties

from .graph import DELETION, INSERTION, cigar_runs

comptime K = 19
comptime W = 19
comptime MID_OCC_FRAC = 2e-4
comptime MIN_MID_OCC = 50
comptime MAX_MID_OCC = 500
comptime MAX_MAX_OCC = 4095
comptime OCC_DIST = 500
comptime Q_OCC_FRAC = 0.01
comptime MAX_GAP = 10000
comptime RMQ_INNER_DIST = 1000
comptime BW = 2000
"""dipcall's -r2k."""
comptime BW_LONG = 100000
comptime MAX_CHAIN_SKIP = 25
comptime MIN_CNT = 3
comptime MIN_CHAIN_SCORE = 40
comptime CHAIN_GAP_SCALE = 0.8
comptime RMQ_RESCUE_SIZE = 1000
comptime RMQ_RESCUE_RATIO = 0.1
comptime MIN_KSW_LEN = 200
comptime LONG_JOIN = 1
comptime IGNORE = 2
comptime TANDEM = 4
comptime NONE = UInt64.MAX
comptime POSITION_BITS = 25
"""Sequences past 2^25 letters are not pinned (the sort key's room)."""
comptime PIN_MEMORY = 2 << 30


@fieldwise_init
struct Seed(ImplicitlyCopyable, TrivialRegisterPassable):
    """A chained seed: its reference and query end letters, its span and minimap2's flags (`TANDEM`, `IGNORE`,
    `LONG_JOIN`)."""

    var x: Int
    var q: Int
    var span: Int
    var flags: Int


struct Pins(Movable):
    """The points (reference letter, query letter) an alignment passes through, in order."""

    var r: List[Int]
    var q: List[Int]

    def __init__(out self):
        self.r = List[Int]()
        self.q = List[Int]()


def _hash64(key: UInt64, mask: UInt64) -> UInt64:
    var k = (~key + (key << 21)) & mask
    k = k ^ (k >> 24)
    k = ((k + (k << 3)) + (k << 8)) & mask
    k = k ^ (k >> 14)
    k = ((k + (k << 2)) + (k << 4)) & mask
    k = k ^ (k >> 28)
    k = (k + (k << 31)) & mask
    return k


def _code(c: UInt8) -> Int:
    if c == UInt8(ord("A")) or c == UInt8(ord("a")):
        return 0
    if c == UInt8(ord("C")) or c == UInt8(ord("c")):
        return 1
    if c == UInt8(ord("G")) or c == UInt8(ord("g")):
        return 2
    if c == UInt8(ord("T")) or c == UInt8(ord("t")):
        return 3
    return 4


struct Sketch(Movable):
    """mm_sketch's minimizers: `x` hash << 8 | span, `y` end << 1 | strand."""

    var x: List[UInt64]
    var y: List[UInt64]

    def __init__(out self):
        self.x = List[UInt64]()
        self.y = List[UInt64]()

    def push(mut self, x: UInt64, y: UInt64):
        self.x.append(x)
        self.y.append(y)


def sketch(s: List[UInt8]) -> Sketch:
    """mm_sketch: the (W, K) minimizers, every k-mer tied for a window's least kept, symmetric k-mers skipped."""
    var out = Sketch()
    var shift1 = UInt64(2 * (K - 1))
    var mask = (UInt64(1) << UInt64(2 * K)) - 1
    var fwd = UInt64(0)
    var rev = UInt64(0)
    var buf_x = List[UInt64](length=W, fill=NONE)
    var buf_y = List[UInt64](length=W, fill=NONE)
    var min_x = NONE
    var min_y = NONE
    var l = 0
    var buf_pos = 0
    var min_pos = 0
    for i in range(len(s)):
        var c = _code(s[i])
        var info_x = NONE
        var info_y = NONE
        if c < 4:
            var span = l + 1 if l + 1 < K else K
            fwd = ((fwd << 2) | UInt64(c)) & mask
            rev = (rev >> 2) | (UInt64(3 ^ c) << shift1)
            if fwd == rev:
                continue
            var z = 0 if fwd < rev else 1
            l += 1
            if l >= K and span < 256:
                info_x = (_hash64(fwd if z == 0 else rev, mask) << 8) | UInt64(span)
                info_y = (UInt64(i) << 1) | UInt64(z)
        else:
            l = 0
        buf_x[buf_pos] = info_x
        buf_y[buf_pos] = info_y
        if l == W + K - 1 and min_x != NONE:
            for j in range(buf_pos + 1, W):
                if min_x == buf_x[j] and buf_y[j] != min_y:
                    out.push(buf_x[j], buf_y[j])
            for j in range(buf_pos):
                if min_x == buf_x[j] and buf_y[j] != min_y:
                    out.push(buf_x[j], buf_y[j])
        if info_x <= min_x:
            if l >= W + K and min_x != NONE:
                out.push(min_x, min_y)
            min_x = info_x
            min_y = info_y
            min_pos = buf_pos
        elif buf_pos == min_pos:
            if l >= W + K - 1 and min_x != NONE:
                out.push(min_x, min_y)
            min_x = NONE
            for j in range(buf_pos + 1, W):
                if min_x >= buf_x[j]:
                    min_x = buf_x[j]
                    min_y = buf_y[j]
                    min_pos = j
            for j in range(buf_pos + 1):
                if min_x >= buf_x[j]:
                    min_x = buf_x[j]
                    min_y = buf_y[j]
                    min_pos = j
            if l >= W + K - 1 and min_x != NONE:
                for j in range(buf_pos + 1, W):
                    if min_x == buf_x[j] and min_y != buf_y[j]:
                        out.push(buf_x[j], buf_y[j])
                for j in range(buf_pos + 1):
                    if min_x == buf_x[j] and min_y != buf_y[j]:
                        out.push(buf_x[j], buf_y[j])
        buf_pos += 1
        if buf_pos == W:
            buf_pos = 0
    if min_x != NONE:
        out.push(min_x, min_y)
    return out^


def _key(a: Seed) -> UInt64:
    """The sort key of minimap2's (x, y) order: reference end, then the tandem flag, the span, the query end."""
    var tandem = UInt64(1) if a.flags & TANDEM else UInt64(0)
    return (UInt64(a.x) << 34) | (tandem << 33) | (UInt64(a.span) << UInt64(POSITION_BITS)) | UInt64(a.q)


def _seed_of(key: UInt64) -> Seed:
    var mask = (UInt64(1) << UInt64(POSITION_BITS)) - 1
    return Seed(
        Int(key >> 34),
        Int(key & mask),
        Int((key >> UInt64(POSITION_BITS)) & 0xFF),
        TANDEM if (key >> 33) & 1 else 0,
    )


def _sorted(seeds: List[Seed]) -> List[Seed]:
    var keys = List[UInt64](capacity=len(seeds))
    for a in seeds:
        keys.append(_key(a))
    sort(keys)
    var out = List[Seed](capacity=len(keys))
    for key in keys:
        out.append(_seed_of(key))
    return out^


@fieldwise_init
struct _Hits(Copyable, Movable):
    """A query minimizer's reference hits: its end and strand (`q_pos`), span, the reference ends and strands, the
    tandem flag and mm_seed_select's filter."""

    var q_pos: UInt64
    var span: Int
    var hits: List[UInt64]
    var tandem: Bool
    var flt: Bool


def _seed_select(mut seeds: List[_Hits], qlen: Int, max_occ: Int):
    """mm_seed_select: in each streak of seeds past `max_occ`, up to its length / OCC_DIST of the rarest kept."""
    var n = len(seeds)
    if n <= 1:
        return
    var high = False
    for s in seeds:
        if len(s.hits) > max_occ:
            high = True
    if not high:
        return
    var last0 = -1
    for i in range(n + 1):
        if i == n or len(seeds[i].hits) <= max_occ:
            if i - last0 > 1:
                var ps = 0 if last0 < 0 else Int(seeds[last0].q_pos >> 1)
                var pe = qlen if i == n else Int(seeds[i].q_pos >> 1)
                var st = last0 + 1
                var en = i
                var max_high = Int(Float64(pe - ps) / Float64(OCC_DIST) + 0.499)
                if max_high > 0:
                    var order = List[UInt64]()
                    for j in range(st, en):
                        order.append((UInt64(len(seeds[j].hits)) << 32) | UInt64(j))
                    sort(order)
                    for m in range(min(min(max_high, 128), len(order))):
                        seeds[Int(order[m] & 0xFFFFFFFF)].flt = True
                for j in range(st, en):
                    seeds[j].flt = not seeds[j].flt
                for j in range(st, en):
                    if len(seeds[j].hits) > MAX_MAX_OCC:
                        seeds[j].flt = True
            last0 = i


def _anchors(reference: List[UInt8], query: List[UInt8]) -> List[Seed]:
    """The forward-strand seeds, minimap2's filters applied, in its order."""
    var r = sketch(reference)
    var slot_of = Dict[UInt64, Int]()
    var slots = List[List[UInt64]]()
    for m in range(len(r.x)):
        var h = r.x[m] >> 8
        var at = slot_of.get(h)
        if at:
            slots[at.value()].append(r.y[m])
        else:
            slot_of[h] = len(slots)
            var hits: List[UInt64] = [r.y[m]]
            slots.append(hits^)
    # mm_idx_cal_max_occ, then minimap2's asm clamp.
    var counts = List[Int](capacity=len(slots))
    for s in slots:
        counts.append(len(s))
    sort(counts)
    var max_occ = MIN_MID_OCC
    if len(counts) > 0:
        max_occ = counts[min(Int((1.0 - MID_OCC_FRAC) * Float64(len(counts))), len(counts) - 1)] + 1
    max_occ = min(max(max_occ, MIN_MID_OCC), MAX_MID_OCC)
    var mv = sketch(query)
    # mm_seed_mz_flt: a minimizer the query itself holds too often.
    if len(mv.x) > max_occ:
        var in_query = Dict[UInt64, Int]()
        for x in mv.x:
            in_query[x] = in_query.get(x, 0) + 1
        var kept = Sketch()
        for m in range(len(mv.x)):
            var c = in_query.get(mv.x[m], 0)
            if not (c > max_occ and Float64(c) > Float64(len(mv.x)) * Q_OCC_FRAC):
                kept.push(mv.x[m], mv.y[m])
        mv = kept^
    var seeds = List[_Hits]()
    for i in range(len(mv.x)):
        var h = mv.x[i] >> 8
        var at = slot_of.get(h)
        if not at:
            continue
        var tandem = (i > 0 and mv.x[i - 1] >> 8 == h) or (i + 1 < len(mv.x) and mv.x[i + 1] >> 8 == h)
        seeds.append(_Hits(mv.y[i], Int(mv.x[i] & 0xFF), slots[at.value()].copy(), tandem, False))
    if MAX_MAX_OCC > max_occ:
        _seed_select(seeds, len(query), max_occ)
    else:
        for ref s in seeds:
            s.flt = len(s.hits) > max_occ
    var out = List[Seed]()
    for s in seeds:
        if s.flt:
            continue
        for hit in s.hits:
            if (hit & 1) == (s.q_pos & 1):
                out.append(Seed(Int(hit >> 1), Int(s.q_pos >> 1), s.span, TANDEM if s.tandem else 0))
    return _sorted(out)


struct _MaxTree(Movable):
    """A segment tree over the anchors in (query end, index) order holding each present anchor's RMQ priority:
    the greatest (priority, query end, index) of a range."""

    var size: Int
    var pri: List[Float64]
    var q: List[Int]
    var j: List[Int]

    def __init__(out self, n: Int):
        self.size = 1
        while self.size < n:
            self.size <<= 1
        self.pri = List[Float64](length=2 * self.size, fill=-1e300)
        self.q = List[Int](length=2 * self.size, fill=-1)
        self.j = List[Int](length=2 * self.size, fill=-1)

    def _ge(self, a: Int, b: Int) -> Bool:
        if self.pri[a] != self.pri[b]:
            return self.pri[a] > self.pri[b]
        if self.q[a] != self.q[b]:
            return self.q[a] > self.q[b]
        return self.j[a] >= self.j[b]

    def _take(mut self, m: Int, from_: Int):
        self.pri[m] = self.pri[from_]
        self.q[m] = self.q[from_]
        self.j[m] = self.j[from_]

    def set(mut self, leaf: Int, pri: Float64, q: Int, j: Int):
        var m = leaf + self.size
        self.pri[m] = pri
        self.q[m] = q
        self.j[m] = j
        m >>= 1
        while m > 0:
            self._take(m, 2 * m if self._ge(2 * m, 2 * m + 1) else 2 * m + 1)
            m >>= 1

    def best(self, lo: Int, hi: Int) -> Int:
        """The index of the greatest present anchor of leaves [lo, hi), -1 for none."""
        var best = -1
        var a = lo + self.size
        var b = hi + self.size
        while a < b:
            if a & 1:
                if self.j[a] >= 0 and (best < 0 or self._ge(a, best)):
                    best = a
                a += 1
            if b & 1:
                b -= 1
                if self.j[b] >= 0 and (best < 0 or self._ge(b, best)):
                    best = b
            a >>= 1
            b >>= 1
        return self.j[best] if best >= 0 else -1


def _bisect_left(keys: List[UInt64], target: UInt64) -> Int:
    var lo = 0
    var hi = len(keys)
    while lo < hi:
        var mid = (lo + hi) // 2
        if keys[mid] < target:
            lo = mid + 1
        else:
            hi = mid
    return lo


def _bisect_right(keys: List[UInt64], target: UInt64) -> Int:
    var lo = 0
    var hi = len(keys)
    while lo < hi:
        var mid = (lo + hi) // 2
        if keys[mid] <= target:
            lo = mid + 1
        else:
            hi = mid
    return lo


def _score(ai: Seed, aj: Seed, pen_gap: Float64, mut exact: Bool, mut width: Int) -> Int:
    """comput_sc_simple (no skip penalty)."""
    var dq = ai.q - aj.q
    var dr = ai.x - aj.x
    var dd = dr - dq if dr > dq else dq - dr
    width = dd
    var dg = min(dr, dq)
    var sc = min(aj.span, dg)
    exact = dd == 0 and dg <= aj.span
    if dd != 0 or dq > aj.span:
        var log_pen = log2(Float64(dd + 1)) if dd >= 1 else 0.0
        sc -= Int(pen_gap * Float64(dd) + 0.5 * log_pen)
    return sc


def _chains(a: List[Seed], bw: Int) -> List[List[Seed]]:
    """mg_lchain_rmq and mg_chain_backtrack: the chains, the best first, each its seeds in order."""
    var n = len(a)
    var chains = List[List[Seed]]()
    if n == 0:
        return chains^
    var pen_gap = CHAIN_GAP_SCALE * 0.01 * Float64(K)
    var max_dist = max(MAX_GAP, bw)
    var max_dist_inner = min(RMQ_INNER_DIST, max_dist)
    var f = List[Int](length=n, fill=0)
    var p = List[Int](length=n, fill=-1)
    var t = List[Int](length=n, fill=0)
    var keys = List[UInt64](capacity=n)
    for j in range(n):
        keys.append((UInt64(a[j].q) << 32) | UInt64(j))
    sort(keys)
    var leaf = List[Int](length=n, fill=0)
    for m in range(n):
        leaf[Int(keys[m] & 0xFFFFFFFF)] = m
    var outer = _MaxTree(n)
    var inner = List[Bool](length=n, fill=False)
    var i0 = 0
    var st = 0
    var st_inner = 0
    for i in range(n):
        var max_j = -1
        var max_f = a[i].span
        if i0 < i and a[i0].x != a[i].x:
            for j in range(i0, i):
                outer.set(leaf[j], Float64(f[j]) + 0.5 * pen_gap * Float64(a[j].x + a[j].q), a[j].q, j)
                inner[j] = True
            i0 = i
        while st < i and a[i].x > a[st].x + max_dist:
            if st < i0:
                outer.set(leaf[st], -1e300, -1, -1)
            st += 1
        while st_inner < i and a[i].x > a[st_inner].x + max_dist_inner:
            inner[st_inner] = False
            st_inner += 1
        var yi = a[i].q
        var lo = 0 if yi - max_dist < 0 else _bisect_right(keys, (UInt64(yi - max_dist) << 32) | UInt64(n))
        var j = outer.best(lo, _bisect_left(keys, UInt64(yi) << 32))
        if j >= 0:
            var exact = False
            var width = 0
            var sc = f[j] + _score(a[i], a[j], pen_gap, exact, width)
            if width <= bw and sc > max_f:
                max_f = sc
                max_j = j
            if not exact and yi > 0:
                var n_skip = 0
                var m = _bisect_right(keys, (UInt64(yi - 1) << 32) | UInt64(n)) - 1
                while m >= 0:
                    var k = Int(keys[m] & 0xFFFFFFFF)
                    m -= 1
                    if a[k].q < yi - max_dist_inner:
                        break
                    if not inner[k]:
                        continue
                    var unused = False
                    var w = 0
                    var s = f[k] + _score(a[i], a[k], pen_gap, unused, w)
                    if w <= bw:
                        if s > max_f:
                            max_f = s
                            max_j = k
                            if n_skip > 0:
                                n_skip -= 1
                        elif t[k] == i:
                            n_skip += 1
                            if n_skip > MAX_CHAIN_SKIP:
                                break
                        if p[k] >= 0:
                            t[p[k]] = i
        f[i] = max_f
        p[i] = max_j
    # mg_chain_backtrack, max_drop = bw.
    var z = List[UInt64]()
    for i in range(n):
        if f[i] >= MIN_CHAIN_SCORE:
            z.append((UInt64(f[i]) << 32) | UInt64(i))
    sort(z)
    for i in range(n):
        t[i] = 0
    var scores = List[Int]()
    var k = len(z) - 1
    while k >= 0:
        var zf = Int(z[k] >> 32)
        var zi = Int(z[k] & 0xFFFFFFFF)
        k -= 1
        if t[zi] != 0:
            continue
        # mg_chain_bk_end
        var i = zi
        var end_i = -1
        var max_i = zi
        var max_s = 0
        while True:
            t[i] = 2
            i = p[i]
            end_i = i
            var s = zf if i < 0 else zf - f[i]
            if s > max_s:
                max_s = s
                max_i = i
            elif max_s - s > bw:
                break
            if not (i >= 0 and t[i] == 0):
                break
        i = zi
        while i >= 0 and i != end_i:
            t[i] = 0
            i = p[i]
        end_i = max_i
        var members = List[Int]()
        i = zi
        while i != end_i:
            members.append(i)
            t[i] = 1
            i = p[i]
        var sc = zf if i < 0 else zf - f[i]
        if sc >= MIN_CHAIN_SCORE and len(members) >= MIN_CNT:
            var chain = List[Seed](capacity=len(members))
            for m in range(len(members) - 1, -1, -1):
                chain.append(a[members[m]])
            # Best first: kept in order of score, a later chain of an equal one after it.
            var at = len(scores)
            while at > 0 and scores[at - 1] < sc:
                at -= 1
            scores.insert(at, sc)
            chains.insert(at, chain^)
    return chains^


def _chain(reference: List[UInt8], query: List[UInt8]) -> List[Seed]:
    """The primary chain, as map.c gets it: re-chained at BW_LONG over every chain's seeds when there are several."""
    var a = _anchors(reference, query)
    var chains = _chains(a, BW)
    if len(chains) > 1:
        # map.c tests the chain compact_a puts first: the chains sorted by their first anchor's reference position,
        # not by score ("u[] may not be sorted by score here").
        var first = 0
        for c in range(1, len(chains)):
            if chains[c][0].x < chains[first][0].x:
                first = c
        var st = chains[first][0].q
        var en = chains[first][len(chains[first]) - 1].q
        if len(query) - (en - st) > RMQ_RESCUE_SIZE or Float64(en - st) > Float64(len(query)) * RMQ_RESCUE_RATIO:
            var pooled = List[Seed]()
            for c in chains:
                for s in c:
                    pooled.append(s)
            chains = _chains(_sorted(pooled), BW_LONG)
    if len(chains) == 0:
        return List[Seed]()
    return chains[0].copy()


def _fix_bad_ends(a: List[Seed], bw: Int, min_match: Int) -> Tuple[Int, Int]:
    """mm_fix_bad_ends: the first seed kept and the seeds' count."""
    var cnt = len(a)
    if cnt < 3:
        return (0, cnt)
    var mlen = a[0].span
    for i in range(1, cnt):
        var tl = a[i].x - a[i - 1].x
        var ql = a[i].q - a[i - 1].q
        mlen += a[i].span if tl > a[i].span and ql > a[i].span else min(tl, ql)
    var first = 0
    var m = a[0].span
    var l = m
    for i in range(1, cnt - 1):
        if a[i].flags & LONG_JOIN:
            break
        var lr = a[i].x - a[i - 1].x
        var lq = a[i].q - a[i - 1].q
        var mn = min(lr, lq)
        if max(lr, lq) - mn > l >> 1:
            first = i
        l += mn
        m += min(mn, a[i].span)
        if l >= bw << 1 or (m >= min_match and m >= bw) or m >= mlen >> 1:
            break
    var count = cnt - first
    m = a[cnt - 1].span
    l = m
    var i = cnt - 2
    while i > first:
        if a[i + 1].flags & LONG_JOIN:
            break
        var lr = a[i + 1].x - a[i].x
        var lq = a[i + 1].q - a[i].q
        var mn = min(lr, lq)
        if max(lr, lq) - mn > l >> 1:
            count = i + 1 - first
        l += mn
        m += min(mn, a[i + 1].span)
        if l >= bw << 1 or (m >= min_match and m >= bw) or m >= mlen >> 1:
            break
        i -= 1
    return (first, count)


def _gap(a: List[Seed], i: Int) -> Int:
    return (a[i].q - a[i - 1].q) - (a[i].x - a[i - 1].x)


def _long_gaps(a: List[Seed], min_gap: Int) -> List[Int]:
    var out = List[Int]()
    for i in range(1, len(a)):
        var g = _gap(a, i)
        if g > min_gap or g < -min_gap:
            out.append(i)
    if len(out) <= 1:
        return List[Int]()
    return out^


def _filter_bad_seeds(mut a: List[Seed]):
    """mm_filter_bad_seeds(min_gap 10, diff_thres 40, max_ext_len MAX_GAP / 2, max_ext_cnt 10): the seeds of a
    stretch whose long insertions and deletions mostly cancel are not pinned."""
    var ks = _long_gaps(a, 10)
    var n = len(ks)
    if n == 0:
        return
    var best = 0
    var max_st = -1
    var max_en = -1
    var k = 0
    while True:
        if k == n or k >= max_en:
            if max_en > 0:
                for i in range(ks[max_st], ks[max_en]):
                    a[i].flags |= IGNORE
            best = 0
            max_st = -1
            max_en = -1
            if k == n:
                break
        var i = ks[k]
        var g = _gap(a, i)
        var n_ins = max(g, 0)
        var n_del = max(-g, 0)
        var qs = a[i - 1].q
        var rs = a[i - 1].x
        var max_diff = 0
        var max_diff_l = -1
        for l in range(k + 1, min(n, k + 11)):
            var j = ks[l]
            if a[j].q - qs > MAX_GAP >> 1 or a[j].x - rs > MAX_GAP >> 1:
                break
            g = _gap(a, j)
            if g > 0:
                n_ins += g
            else:
                n_del -= g
            var diff = n_ins + n_del - abs(n_ins - n_del)
            if max_diff < diff:
                max_diff = diff
                max_diff_l = l
        if max_diff > 40 and max_diff > best:
            best = max_diff
            max_st = k
            max_en = max_diff_l
        k += 1


def _filter_bad_seeds_alt(mut a: List[Seed]):
    """mm_filter_bad_seeds_alt(min_gap 30, max_ext MAX_GAP / 2): long gaps close together joined into one."""
    var ks = _long_gaps(a, 30)
    var n = len(ks)
    var k = 0
    while k < n:
        var i = ks[k]
        var gap1 = abs(_gap(a, i))
        var re1 = a[i].x
        var qe1 = a[i].q
        var l = k + 1
        while l < n:
            var j = ks[l]
            if a[j].q - qe1 > MAX_GAP >> 1 or a[j].x - re1 > MAX_GAP >> 1:
                break
            var gap2 = abs(_gap(a, j))
            var pre = a[j - 1].span
            if min(a[j - 1].x + pre - re1, a[j - 1].q + pre - qe1) > gap1 + gap2:
                break
            re1 = a[j].x
            qe1 = a[j].q
            gap1 = gap2
            l += 1
        if l > k + 1:
            var end = ks[l - 1]
            for j in range(ks[k], end):
                a[j].flags |= IGNORE
            a[end].flags |= LONG_JOIN
        k = l


def pins(reference: List[UInt8], query: List[UInt8]) -> Pins:
    """mm_align1's pins for `query` against `reference`: the middles of the primary chain's kept seeds where both
    sequences advanced MIN_KSW_LEN letters or more since the last pin (or a long join), the last seed's always."""
    var out = Pins()
    if len(reference) >= 1 << POSITION_BITS or len(query) >= 1 << POSITION_BITS:
        return out^
    var chain = _chain(reference, query)
    if getenv("PLEGMA_PINS_TRACE") == "1":
        print("pins chain", len(chain), file=stderr)
    if len(chain) == 0:
        return out^
    var kept = _fix_bad_ends(chain, BW, MIN_CHAIN_SCORE * 2)
    var a = List[Seed](capacity=kept[1])
    for i in range(kept[0], kept[0] + kept[1]):
        a.append(chain[i])
    _filter_bad_seeds(a)
    _filter_bad_seeds_alt(a)
    var half = K >> 1
    var rs = a[0].x - half
    var qs = a[0].q - half
    var points_r = List[Int]()
    var points_q = List[Int]()
    points_r.append(rs)
    points_q.append(qs)
    for i in range(1, len(a)):
        var last = i == len(a) - 1
        if a[i].flags & (IGNORE | TANDEM) and not last:
            continue
        var re = a[i].x - half
        var qe = a[i].q - half
        if last or a[i].flags & LONG_JOIN or (qe - qs >= MIN_KSW_LEN and re - rs >= MIN_KSW_LEN):
            points_r.append(re)
            points_q.append(qe)
            rs = re
            qs = qe
    for m in range(len(points_r)):
        if 0 < points_r[m] < len(reference) and 0 < points_q[m] < len(query):
            out.r.append(points_r[m])
            out.q.append(points_q[m])
    return out^


def pins_in_context(reference: List[UInt8], query: List[UInt8], left: List[UInt8], right: List[UInt8]) -> Pins:
    """The pins of `left + reference + right` against `left + query + right` (the flanks context for the seeds and the
    chain only) that fall inside `reference` and `query`, in their own coordinates, each past the last on both."""
    var whole_ref = left.copy()
    whole_ref.extend(reference.copy())
    whole_ref.extend(right.copy())
    var whole_query = left.copy()
    whole_query.extend(query.copy())
    whole_query.extend(right.copy())
    var found = pins(whole_ref, whole_query)
    var out = Pins()
    var n = len(left)
    for m in range(len(found.r)):
        var r = found.r[m] - n
        var q = found.q[m] - n
        if 0 < r < len(reference) and 0 < q < len(query):
            if len(out.r) == 0 or (r >= out.r[len(out.r) - 1] and q >= out.q[len(out.q) - 1]):
                out.r.append(r)
                out.q.append(q)
    return out^


def pinned_pairing(
    mut aligner: Aligner,
    reference: List[UInt8],
    query: List[UInt8],
    left: List[UInt8],
    right: List[UInt8],
    costs: Costs,
) raises -> List[Int]:
    """For each letter of `query`, the letter of `reference` its pinned alignment pairs it with (-1: inserted): the
    stretches between the pins (`pins_in_context`) each aligned exactly at `costs`, gaps left."""
    var found = pins_in_context(reference, query, left, right)
    if getenv("PLEGMA_PINS_TRACE") == "1":
        var text = String("pins")
        for m in range(len(found.r)):
            text += " " + String(found.r[m]) + "," + String(found.q[m])
        print(text, file=stderr)
    var points_r: List[Int] = [0]
    var points_q: List[Int] = [0]
    points_r.extend(found.r.copy())
    points_q.extend(found.q.copy())
    points_r.append(len(reference))
    points_q.append(len(query))
    var out = List[Int](capacity=len(query))
    for m in range(len(points_r) - 1):
        var r0 = points_r[m]
        var r1 = points_r[m + 1]
        var q0 = points_q[m]
        var q1 = points_q[m + 1]
        if q1 == q0:
            continue
        if r1 == r0:
            for _ in range(q0, q1):
                out.append(-1)
            continue
        var window_r = List[UInt8](capacity=r1 - r0)
        for i in range(r0, r1):
            window_r.append(reference[i])
        var window_q = List[UInt8](capacity=q1 - q0)
        for j in range(q0, q1):
            window_q.append(query[j])
        var aligned = aligner.align(
            Span(window_r), Span(window_q), costs, Mode.GLOBAL, ties=Ties.LEFT, max_memory=PIN_MEMORY
        )
        var runs = cigar_runs(aligned.cigar)
        var i = r0
        for x in range(len(runs.ops)):
            var op = runs.ops[x]
            for _ in range(runs.counts[x]):
                if op == DELETION:
                    i += 1
                elif op == INSERTION:
                    out.append(-1)
                else:
                    out.append(i)
                    i += 1
    return out^
