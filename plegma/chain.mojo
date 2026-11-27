"""minimap2 2.31 RMQ seed chaining (`mg_lchain_rmq`, `mg_chain_backtrack`, and long-join re-chaining)."""

from std.memory import bitcast

from .seeds import K, ReferenceIndex, Seed, collect_anchors, radix_sort, sorted_seeds

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
comptime PENALTIES = 4096
"""The gap penalties `rmq_chains` tabulates: every dd the inner window's anchors make, which stays within it."""


@fieldwise_init
struct _Leaf(ImplicitlyCopyable, TrivialRegisterPassable):
    """An anchor as `rmq_chains`'s inner loop reads it, by leaf: its seed, score and predecessor's leaf, minimap2's
    skip mark t[], and the inner window's present leaf before it."""

    var x: Int32
    var q: Int32
    var span: Int32
    var f: Int32
    var j: Int32
    var p: Int32
    var t: Int32
    var before: Int32


@fieldwise_init
struct _Top(ImplicitlyCopyable, TrivialRegisterPassable):
    """A `_MaxTree` node: the greatest (priority, leaf) under it; an absent leaf is the lowest priority and leaf -1."""

    var pri: Float64
    var leaf: Int

    @inline(.always)
    def __ge__(self, other: Self) -> Bool:
        return self.pri > other.pri or (self.pri == other.pri and self.leaf >= other.leaf)


comptime _ABSENT = _Top(pri=-1e300, leaf=-1)


struct _MaxTree(Movable):
    """A segment tree over the leaves (the anchors in (query end, index) order) holding each present leaf's RMQ
    priority: a node keeps the greatest (priority, leaf), i.e. ties to the greater (query end, index)."""

    var size: Int
    var top: List[_Top]

    def __init__(out self, n: Int):
        self.size = 1
        while self.size < n:
            self.size <<= 1
        self.top = List[_Top](length=2 * self.size, fill=_ABSENT)

    def set(mut self, leaf: Int, top: _Top):
        var p = self.top.unsafe_ptr()
        var m = leaf + self.size
        p[unsafe_offset=m] = top
        m >>= 1
        while m > 0:
            var l = p[unsafe_offset=2 * m]
            var r = p[unsafe_offset=2 * m + 1]
            var best = l if l >= r else r
            if best.leaf == p[unsafe_offset=m].leaf and best.pri == p[unsafe_offset=m].pri:
                return
            p[unsafe_offset=m] = best
            m >>= 1

    def best(self, lo: Int, hi: Int) -> Int:
        """The greatest present leaf of leaves [lo, hi), -1 for none."""
        var p = self.top.unsafe_ptr()
        var best = _ABSENT
        var a = lo + self.size
        var b = hi + self.size
        while a < b:
            if a & 1:
                if p[unsafe_offset=a] >= best:
                    best = p[unsafe_offset=a]
                a += 1
            if b & 1:
                b -= 1
                if p[unsafe_offset=b] >= best:
                    best = p[unsafe_offset=b]
            a >>= 1
            b >>= 1
        return best.leaf


struct _Present(Movable):
    """Which leaves (the anchors in (query end, index) order) are present: a tree of counts, for the present leaf at
    or before a leaf without walking the absent ones between."""

    var size: Int
    var count: List[Int32]

    def __init__(out self, n: Int):
        self.size = 1
        while self.size < n:
            self.size <<= 1
        self.count = List[Int32](length=2 * self.size, fill=0)

    def set(mut self, leaf: Int, present: Bool):
        var m = leaf + self.size
        var change = Int32(1) if present else Int32(-1)
        if (self.count[m] > 0) == present:
            return
        while m > 0:
            self.count[m] += change
            m >>= 1

    def at_or_after(self, leaf: Int) -> Int:
        """The first present leaf at or after `leaf`, -1 for none."""
        if leaf < 0 or leaf >= self.size:
            return -1
        var v = leaf + self.size
        if self.count[v] > 0:
            return leaf
        while v > 1:
            if not (v & 1) and self.count[v + 1] > 0:
                v += 1
                while v < self.size:
                    v = 2 * v if self.count[2 * v] > 0 else 2 * v + 1
                return v - self.size
            v >>= 1
        return -1

    def at_or_before(self, leaf: Int) -> Int:
        """The last present leaf at or before `leaf`, -1 for none."""
        if leaf < 0:
            return -1
        var v = leaf + self.size
        if self.count[v] > 0:
            return leaf
        while v > 1:
            if v & 1 and self.count[v - 1] > 0:
                v -= 1
                while v < self.size:
                    v = 2 * v + 1 if self.count[2 * v + 1] > 0 else 2 * v
                return v - self.size
            v >>= 1
        return -1


@inline(.always)
def _mg_log2(x: Float32) -> Float32:
    """minimap2's mg_log2: log2 from the float's exponent and a polynomial in its mantissa (not exact; x >= 2)."""
    var i = bitcast[DType.uint32](x)
    var log_2 = Float32(Int((i >> 23) & 255) - 128)
    i &= ~(UInt32(255) << 23)
    i += UInt32(127) << 23
    var z = bitcast[DType.float32](i)
    return log_2 + (Float32(-0.34484843) * z + Float32(2.02466578)) * z - Float32(0.67487759)


@inline(.always)
def _penalty(dd: Int, pen_gap: Float32) -> Int:
    """comput_sc_simple's penalty for a gap difference `dd` (0 for none)."""
    if dd == 0:
        return 0
    return Int(pen_gap * Float32(dd) + Float32(0.5) * _mg_log2(Float32(dd + 1)))


@inline(.never)
def _penalty_far(dd: Int, pen_gap: Float32) -> Int:
    """`_penalty` past the table, out of line so the inner loop's branch stays a table read."""
    return _penalty(dd, pen_gap)


@inline(.always)
def _score(ai: Seed, aj: Seed, pen_gap: Float32, mut exact: Bool, mut width: Int) -> Int:
    """comput_sc_simple (no skip penalty), in floats as minimap2 computes it."""
    var dq = ai.q - aj.q
    var dr = ai.x - aj.x
    var dd = dr - dq if dr > dq else dq - dr
    width = dd
    var dg = min(dr, dq)
    var sc = min(aj.span, dg)
    exact = dd == 0 and dg <= aj.span
    if dd != 0 or dq > aj.span:
        sc -= _penalty(dd, pen_gap)
    return sc


@inline(.never)
def _walk(
    leaves: MutPointer[_Leaf, _],
    penalty: ImmPointer[Int32, _],
    start: Int,
    i: Int,
    xi: Int,
    yi: Int,
    bw: Int,
    max_dist_inner: Int,
    pen_gap: Float32,
    var max_f: Int,
    var max_j: Int,
) -> Tuple[Int, Int]:
    """mg_lchain_rmq's inner loop for anchor `i`: the inner window's anchors from leaf `start`, query end descending
    (minimap2's inner tree, walked back), each scored by `_score` with its penalty read from the table; the best
    score and predecessor, as given if none beats them. Out of line, the loop keeps its values in registers."""
    var n_skip = 0
    var floor = yi - max_dist_inner
    var m = start
    while m >= 0:
        ref here = leaves[unsafe_offset=m]
        m = Int(here.before)
        var yj = Int(here.q)
        if yj < floor:
            break
        var dq = yi - yj
        var dr = xi - Int(here.x)
        var dd = dr - dq if dr > dq else dq - dr
        if dd > bw:
            continue
        var s = Int(here.f) + min(Int(here.span), min(dr, dq))
        if dd < PENALTIES:
            s -= Int(penalty[unsafe_offset=dd])
        else:
            s -= _penalty_far(dd, pen_gap)
        if s > max_f:
            max_f = s
            max_j = Int(here.j)
            if n_skip > 0:
                n_skip -= 1
        elif Int(here.t) == i:
            n_skip += 1
            if n_skip > MAX_CHAIN_SKIP:
                break
        var prior = Int(here.p)
        if prior >= 0:
            leaves[unsafe_offset=prior].t = Int32(i)
    return (max_f, max_j)


def rmq_chains(a: List[Seed], bw: Int) -> List[List[Seed]]:
    """mg_lchain_rmq and mg_chain_backtrack: the chains, the best first, each its seeds in order."""
    var n = len(a)
    var chains = List[List[Seed]]()
    if n == 0:
        return chains^
    var pen_gap = Float32(CHAIN_GAP_SCALE * 0.01 * Float64(K))
    var max_dist = max(MAX_GAP, bw)
    var max_dist_inner = min(RMQ_INNER_DIST, max_dist)
    var penalty = List[Int32](capacity=PENALTIES)
    for dd in range(PENALTIES):
        penalty.append(Int32(_penalty(dd, pen_gap)))
    var f = List[Int](length=n, fill=0)
    var p = List[Int](length=n, fill=-1)
    var t = List[Int](length=n, fill=0)
    var keys = List[UInt64](capacity=n)
    for j in range(n):
        keys.append((UInt64(a[j].q) << 32) | UInt64(j))
    radix_sort(keys)
    var leaf = List[Int](length=n, fill=0)
    var leaves = List[_Leaf](capacity=n)
    var first = List[Int](length=n, fill=0)
    var near = List[Int](length=n, fill=0)
    var run = 0
    var lo = 0
    for m in range(n):
        var j = Int(keys[m] & 0xFFFFFFFF)
        leaf[j] = m
        leaves.append(
            _Leaf(x=Int32(a[j].x), q=Int32(a[j].q), span=Int32(a[j].span), f=0, j=Int32(j), p=-1, t=0, before=-1)
        )
        if a[j].q != Int(keys[run] >> 32):
            run = m
        while Int(keys[lo] >> 32) <= a[j].q - max_dist:
            lo += 1
        first[m] = run
        near[m] = lo
    var outer = _MaxTree(n)
    var inner = _Present(n)
    var after = List[Int](length=n, fill=-1)
    var i0 = 0
    var st = 0
    var st_inner = 0
    var leaves_p = leaves.unsafe_ptr()
    var penalty_p = penalty.unsafe_ptr()
    for i in range(n):
        var max_j = -1
        var max_f = a[i].span
        if i0 < i and a[i0].x != a[i].x:
            for j in range(i0, i):
                var m = leaf[j]
                outer.set(m, _Top(pri=Float64(f[j]) + 0.5 * Float64(pen_gap) * Float64(a[j].x + a[j].q), leaf=m))
                var prior = inner.at_or_before(m - 1)
                var next = after[prior] if prior >= 0 else inner.at_or_after(m + 1)
                leaves_p[unsafe_offset=m].before = Int32(prior)
                after[m] = next
                if prior >= 0:
                    after[prior] = m
                if next >= 0:
                    leaves_p[unsafe_offset=next].before = Int32(m)
                inner.set(m, True)
            i0 = i
        while st < i and a[i].x > a[st].x + max_dist:
            if st < i0:
                outer.set(leaf[st], _ABSENT)
            st += 1
        while st_inner < i and a[i].x > a[st_inner].x + max_dist_inner:
            var m = leaf[st_inner]
            var prior = Int(leaves_p[unsafe_offset=m].before)
            var next = after[m]
            if prior >= 0:
                after[prior] = next
            if next >= 0:
                leaves_p[unsafe_offset=next].before = Int32(prior)
            leaves_p[unsafe_offset=m].before = -1
            after[m] = -1
            inner.set(m, False)
            st_inner += 1
        var yi = a[i].q
        var mi = leaf[i]
        var best = outer.best(near[mi], first[mi])
        if best >= 0:
            var j = Int(leaves_p[unsafe_offset=best].j)
            var exact = False
            var width = 0
            var sc = f[j] + _score(a[i], a[j], pen_gap, exact, width)
            if width <= bw and sc > max_f:
                max_f = sc
                max_j = j
            if not exact and yi > 0:
                var walked = _walk(
                    leaves_p,
                    penalty_p,
                    inner.at_or_before(first[mi] - 1),
                    i,
                    a[i].x,
                    yi,
                    bw,
                    max_dist_inner,
                    pen_gap,
                    max_f,
                    max_j,
                )
                max_f = walked[0]
                max_j = walked[1]
        f[i] = max_f
        p[i] = max_j
        leaves_p[unsafe_offset=mi].f = Int32(max_f)
        leaves_p[unsafe_offset=mi].p = Int32(leaf[max_j]) if max_j >= 0 else -1
    var z = List[UInt64]()
    for i in range(n):
        if f[i] >= MIN_CHAIN_SCORE:
            z.append((UInt64(f[i]) << 32) | UInt64(i))
    radix_sort(z)
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
        var i = zi
        var end_i: Int
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
            var at = len(scores)
            while at > 0 and scores[at - 1] < sc:
                at -= 1
            scores.insert(at, sc)
            chains.insert(at, chain^)
    return chains^


def primary_chain(a: List[Seed], query_len: Int) -> List[Seed]:
    """The primary chain of `a`, as map.c gets it: re-chained at `BW_LONG` over every chain's seeds when there are
    several."""
    var chains = rmq_chains(a, BW)
    if len(chains) > 1:
        var first = 0
        for c in range(1, len(chains)):
            if chains[c][0].x < chains[first][0].x:
                first = c
        var st = chains[first][0].q
        var en = chains[first][len(chains[first]) - 1].q
        if query_len - (en - st) > RMQ_RESCUE_SIZE or Float64(en - st) > Float64(query_len) * RMQ_RESCUE_RATIO:
            var pooled = List[Seed]()
            for c in chains:
                for s in c:
                    pooled.append(s)
            chains = rmq_chains(sorted_seeds(pooled), BW_LONG)
    if len(chains) == 0:
        return List[Seed]()
    return chains[0].copy()


def primary_chain(index: ReferenceIndex, query: List[UInt8]) -> List[Seed]:
    """The primary chain of `query` against pre-indexed `index`."""
    return primary_chain(collect_anchors(index, query), len(query))


def primary_chain(reference: List[UInt8], query: List[UInt8]) -> List[Seed]:
    """The primary chain of `query` against `reference`."""
    return primary_chain(collect_anchors(reference, query), len(query))
