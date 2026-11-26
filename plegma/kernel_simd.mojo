"""SIMD lane primitives and predecessor/band helpers for `plegma/kernel.mojo`."""

from std.math import iota
from std.memory import bitcast
from std.sys import simd_width_of, size_of

from .graph import Graph


@inline(.always)
def _neg[dt: DType]() -> Scalar[dt]:
    """The score of a cell nothing reaches: far below every real score, far enough above the type's floor that
    a few gap costs taken from it cannot wrap (the kernel also clamps at it)."""
    comptime if dt == DType.int16:
        return Scalar[dt](-16000)
    else:
        return Scalar[dt](-(1 << 29))


@inline(.always)
def _fill[dt: DType](p: MutPointer[Scalar[dt], _], count: Int):
    comptime V = simd_width_of[dt]()
    var i = 0
    while i + V <= count:
        p.unsafe_store(i, SIMD[dt, V](_neg[dt]()))
        i += V
    while i < count:
        p[unsafe_offset=i] = _neg[dt]()
        i += 1


@inline(.always)
def _zero[dt: DType](p: MutPointer[Scalar[dt], _], count: Int):
    comptime V = simd_width_of[dt]()
    var i = 0
    while i + V <= count:
        p.unsafe_store(i, SIMD[dt, V](0))
        i += V
    while i < count:
        p[unsafe_offset=i] = 0
        i += 1


@inline(.always)
def _predecessor(
    pred_rows: List[Int], pred_at: List[Int], choices: List[UInt8], choice_at: List[Int], r: Int, cell: Int, kind: Int
) -> Int:
    """The row a move at `cell` of row `r` came from (`kind`: 0 the diagonal, 1 and 2 the deletions)."""
    var which = 0 if choice_at[r] < 0 else Int(choices[choice_at[r] + 3 * cell + kind])
    return pred_rows[pred_at[r] + which]


def _heaviest_first(graph: Graph, node: Int, left: Int, row_of: List[Int], mut rows: List[Int], mut keys: List[UInt64]):
    """Appends the rows of `node`'s predecessors in the alignment, heaviest edge first (the earlier edge first
    among equals): row 0 for `left` (or for none at all), so of equally good moves the alignment takes the one
    most sequences took. `keys` is scratch."""
    var e = Int(graph.first_in[node])
    if graph.in_degree[node] == 1:
        var p = Int(graph.edge_from[e])
        rows.append(0 if p == left or row_of[p] < 0 else row_of[p])
        return
    keys.clear()
    while e >= 0:
        var p = Int(graph.edge_from[e])
        if p == left or row_of[p] >= 0:
            keys.append((UInt64((1 << 31) - Int(graph.edge_weight[e])) << 32) | UInt64(e))
        e = Int(graph.next_in[e])
    if len(keys) == 0:
        rows.append(0)
        return
    # Equal weights keep the edges' order: edge ids grow as the edges are made. The keys are few and distinct: an
    # insertion sort, the order `sort` gives.
    for a in range(1, len(keys)):
        var key = keys[a]
        var b = a
        while b > 0 and keys[b - 1] > key:
            keys[b] = keys[b - 1]
            b -= 1
        keys[b] = key
    for key in keys:
        var p = Int(graph.edge_from[Int(key & 0xFFFFFFFF)])
        rows.append(0 if p == left else row_of[p])


@inline(.always)
def _reach(band_b: Int, band_f: Float64, band_max: Int, m: Int) -> Int:
    """The band's reach either side for a query of `m` letters: `band_b + band_f m`, at most `band_max` when that is
    not negative."""
    var reach = band_b + Int(band_f * Float64(m))
    return min(reach, band_max) if band_max >= 0 else reach


@inline(.always)
def _drop(mut rows: List[Int], row: Int):
    """`row` out of `rows` (order not kept)."""
    for k in range(len(rows)):
        if rows[k] == row:
            rows[k] = rows[len(rows) - 1]
            _ = rows.pop()
            return


@inline(.always)
def _beats(score: Int, base: Int, best: Int, best_ok: Bool, neg: Int) -> Bool:
    """Whether a cell storing `score` in a row on `base` beats the best end so far: a reached cell (above half the
    sentinel, `neg`) beats an unreached one, then the higher absolute score."""
    var ok = score >= neg // 2
    if ok != best_ok:
        return ok
    return score + base > best


@inline(.always)
def _moved[dt: DType](source: MutPointer[Scalar[dt], _], into: MutPointer[Scalar[dt], _], count: Int, delta: Int):
    """`count` scores from `source` into `into`, each plus `delta` (a lower base's to a higher one's: delta <= 0),
    cells at the sentinel kept at it. A delta past half the lanes' range (two predecessors' bases tens of thousands
    apart in 16 bits) is added in full width and clamped at the sentinel: narrowed, it wrapped and the lower row's
    scores came in far too high."""
    comptime V = simd_width_of[dt]()
    comptime Vec = SIMD[dt, V]
    if delta < Int(_neg[dt]()) or delta > -Int(_neg[dt]()):
        var floor = Int(_neg[dt]())
        for c in range(count):
            var v = source[unsafe_offset=c]
            into[unsafe_offset=c] = _neg[dt]() if v <= _neg[dt]() else Scalar[dt](
                min(max(Int(v) + delta, floor), -2 * floor)
            )
        return
    var shift = Scalar[dt](delta)
    var neg = Vec(_neg[dt]())
    var i = 0
    while i + V <= count:
        var v = source.unsafe_load[width=V](i)
        into.unsafe_store(i, v.le(neg).select(neg, _floor(v + Vec(shift))))
        i += V
    while i < count:
        var v = source[unsafe_offset=i]
        into[unsafe_offset=i] = _neg[dt]() if v <= _neg[dt]() else _floor(v + shift)
        i += 1


@inline(.always)
def _rebased[dt: DType](row: MutPointer[Scalar[dt], _], span: Int, shift: Scalar[dt]):
    """Cells 0 ..< span of `row` less `shift` (the row moving to a base `shift` higher), cells at the sentinel kept
    at it."""
    comptime V = simd_width_of[dt]()
    comptime Vec = SIMD[dt, V]
    var neg = Vec(_neg[dt]())
    var i = 0
    while i + V <= span:
        var v = row.unsafe_load[width=V](i)
        row.unsafe_store(i, v.le(neg).select(neg, _floor(v - Vec(shift))))
        i += V
    while i < span:
        var v = row[unsafe_offset=i]
        row[unsafe_offset=i] = _neg[dt]() if v <= _neg[dt]() else _floor(v - shift)
        i += 1


@inline(.always)
def _first_equal[dt: DType](h: MutPointer[Scalar[dt], _], span: Int, value: Scalar[dt]) -> Int:
    """The first of cells 0 ..< span of `h` holding `value` (one does), a vector at a time."""
    comptime V = simd_width_of[dt]()
    var lane = iota[dt, V]()
    var at = 0
    while at + V <= span:
        var hit = h.unsafe_load[width=V](at).eq(SIMD[dt, V](value))
        if hit.reduce_or():
            return at + Int(hit.select(lane, SIMD[dt, V](Scalar[dt](V))).reduce_min())
        at += V
    while h[unsafe_offset=at] != value:
        at += 1
    return at


@inline(.always)
def _last_equal[dt: DType](h: MutPointer[Scalar[dt], _], span: Int, value: Scalar[dt]) -> Int:
    """The last of cells 0 ..< span of `h` holding `value` (one does), a vector at a time from the right."""
    comptime V = simd_width_of[dt]()
    var lane = iota[dt, V]()
    var at = span - V
    while at >= 0:
        var hit = h.unsafe_load[width=V](at).eq(SIMD[dt, V](value))
        if hit.reduce_or():
            return at + Int(hit.select(lane, SIMD[dt, V](Scalar[dt](-1))).reduce_max())
        at -= V
    var k = min(span, at + V) - 1
    while h[unsafe_offset=k] != value:
        k -= 1
    return k


@inline(.always)
def _prefix_max[dt: DType, width: Int](v: SIMD[dt, width]) -> SIMD[dt, width]:
    """Each lane's maximum over itself and the lanes before it."""
    var floor = SIMD[dt, width](_neg[dt]())
    var out = v
    comptime for step in range(6):
        comptime shift = 1 << step
        comptime if shift < width:
            out = max(out, _up[shift](out, floor))
    return out


@inline(.always)
def _up[k: Int, dt: DType, width: Int](v: SIMD[dt, width], fill: SIMD[dt, width]) -> SIMD[dt, width]:
    """`v` moved `k` lanes up (lane i + k holds lane i), `fill`'s value (in every lane of it) in the `k` below. 16-bit
    lanes move as pairs in 32-bit ones and shifts within them: a 16-bit shuffle across a vector takes two or three
    micro-ops on the shuffle port (Skylake-X), a 32-bit one one, and the shifts run elsewhere."""
    comptime if size_of[dt]() == 2 and width % 2 == 0:
        comptime half = width // 2
        var pairs = bitcast[DType.uint32, half](v)
        var filled = bitcast[DType.uint32, half](fill)
        comptime if k % 2 == 0:
            return bitcast[dt, width](filled.join(pairs).slice[half, offset=half - k // 2]())
        else:
            var high = filled.join(pairs).slice[half, offset=half - (k + 1) // 2]()
            var low = filled.join(pairs).slice[half, offset=half - k // 2]()
            return bitcast[dt, width]((low << 16) | (high >> 16))
    else:
        return fill.join(v).slice[width, offset=width - k]()


@inline(.always)
def _last[dt: DType, width: Int](v: SIMD[dt, width]) -> SIMD[dt, width]:
    """`v`'s last lane in every lane; 16-bit lanes by a 32-bit broadcast, as `_up` moves them."""
    comptime if size_of[dt]() == 2 and width % 2 == 0:
        comptime half = width // 2
        var top = SIMD[DType.uint32, half](bitcast[DType.uint32, half](v)[half - 1])
        return bitcast[dt, width]((top >> 16) | (top & 0xFFFF0000))
    else:
        return SIMD[dt, width](v[width - 1])


@inline(.always)
def _floor[dt: DType, width: Int](v: SIMD[dt, width]) -> SIMD[dt, width]:
    """`v` kept at or above the sentinel in 16-bit lanes, where a few gap costs taken from it would wrap."""
    comptime if dt == DType.int16:
        return max(v, SIMD[dt, width](_neg[dt]()))
    else:
        return v


@inline(.always)
def _floor[dt: DType](v: Scalar[dt]) -> Scalar[dt]:
    """The scalar `_floor`."""
    comptime if dt == DType.int16:
        return max(v, _neg[dt]())
    else:
        return v
