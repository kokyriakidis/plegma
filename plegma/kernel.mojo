"""Aligning a sequence to the whole graph: banded dynamic programming over the graph's nodes in topological order.

Scores are maximized, as abPOA's: a match earns `match`, a mismatch costs `mismatch`, a gap of `k` letters costs
the less of `open1 + extend1 k` and `open2 + extend2 k` (convex; `open2` 0 makes it affine, `open1` 0 linear).
Each node's row covers only its band of query positions: abPOA's adaptive band, `b + f m` letters either side of
the predecessors' best cells and of the diagonal the node's remaining path (along the heaviest edges) leaves.

Memory follows the band and the graph's width, not the graph times the query. A row's scores live only until its
last successor is computed; what the traceback needs is one byte per cell (where each cell's score came from), and
at nodes with several predecessors three more (which predecessor each move came from).
"""

from std.sys import simd_width_of

from .alphabet import encode_base
from .graph import Graph
from .poa import SINK, SOURCE


comptime GLOBAL = 0
"""The whole query against a path from a graph source to a graph sink."""
comptime LOCAL = 1
"""The best-scoring part of the query against any part of a path."""
comptime EXTENSION = 2
"""From the graph's sources and the query's start, as far as it pays."""

comptime _NEG = Int32(-(1 << 29))
comptime _V = simd_width_of[DType.int32]()
comptime _I = SIMD[DType.int32, _V]

# A cell's traceback byte: where its score came from (the low three bits), and whether each gap state extended
# a gap (rather than opening one).
comptime _DIAGONAL = UInt8(0)
comptime _DELETION1 = UInt8(1)
comptime _DELETION2 = UInt8(2)
comptime _INSERTION1 = UInt8(3)
comptime _INSERTION2 = UInt8(4)
comptime _START = UInt8(5)
comptime _SOURCE = UInt8(7)
comptime _D1_EXTENDS = UInt8(8)
comptime _D2_EXTENDS = UInt8(16)
comptime _F1_EXTENDS = UInt8(32)
comptime _F2_EXTENDS = UInt8(64)


@inline(.always)
def _fill(p: MutPointer[Int32, _], count: Int):
    var i = 0
    while i + _V <= count:
        p.unsafe_store(i, SIMD[DType.int32, _V](_NEG))
        i += _V
    while i < count:
        p[unsafe_offset=i] = _NEG
        i += 1


struct Scoring(ImplicitlyCopyable):
    """abPOA's scoring: `match`, `mismatch`, and gaps of `k` letters at the less of `open1 + extend1 k` and
    `open2 + extend2 k` (convex), `open1 + extend1 k` when `open2` is 0 (affine), `extend1 k` when `open1` is 0
    (linear)."""

    var match_score: Int32
    var mismatch: Int32
    var open1: Int32
    var extend1: Int32
    var open2: Int32
    var extend2: Int32

    def __init__(
        out self,
        match_score: Int = 2,
        mismatch: Int = 4,
        open1: Int = 4,
        extend1: Int = 2,
        open2: Int = 24,
        extend2: Int = 1,
    ):
        self.match_score = Int32(match_score)
        self.mismatch = Int32(mismatch)
        self.open1 = Int32(open1)
        self.extend1 = Int32(extend1)
        self.open2 = Int32(open2)
        self.extend2 = Int32(extend2)

    def convex(self) -> Bool:
        return self.open1 > 0 and self.open2 > 0

    def gap(self, k: Int) -> Int32:
        """What a gap of `k` letters costs."""
        if k == 0:
            return 0
        var first = self.open1 + self.extend1 * Int32(k)
        return min(first, self.open2 + self.extend2 * Int32(k)) if self.convex() else first


struct GraphAlignment(Movable):
    """For each query letter, the node it is aligned to (-1: inserted, or outside a local alignment), the
    alignment's score, the DP cells computed and the most bytes the DP held at once."""

    var nodes: List[Int]
    var score: Int
    var cells: Int
    var peak_bytes: Int

    def __init__(out self, var nodes: List[Int], score: Int, cells: Int, peak_bytes: Int):
        self.nodes = nodes^
        self.score = score
        self.cells = cells
        self.peak_bytes = peak_bytes


struct Topology(Movable):
    """The graph's nodes in topological order, each node's rank, and the length of its remaining path along
    the heaviest out-edges (abPOA's band diagonal)."""

    var order: List[Int]
    var rank: List[Int]
    var remaining: List[Int]

    def __init__(out self, graph: Graph) raises:
        self.rank = graph.node_order()
        var n = len(self.rank)
        self.order = List[Int](length=n, fill=0)
        for node in range(n):
            self.order[self.rank[node]] = node
        self.remaining = List[Int](length=n, fill=0)
        for r in reversed(range(n)):
            var node = self.order[r]
            var best_weight = -1
            for i in range(len(graph.successors[node])):
                var s = graph.successors[node][i]
                var weight = graph.weights[node][i]
                if weight > best_weight or (weight == best_weight and self.remaining[s] + 1 > self.remaining[node]):
                    best_weight = weight
                    self.remaining[node] = self.remaining[s] + 1


def align_to_graph(
    graph: Graph,
    topology: Topology,
    query: Span[UInt8, _],
    scoring: Scoring,
    mode: Int = GLOBAL,
    band_b: Int = 10,
    band_f: Float64 = 0.01,
) raises -> GraphAlignment:
    """The best alignment of `query` to `graph` under `scoring` in `mode`, inside abPOA's adaptive band
    (`band_b` < 0: no band; a band the global end falls outside of is doubled and tried again). Of equal
    alignments, the one taking the diagonal first in the traceback, so gaps come first, and of equal
    predecessors the one whose edge most sequences took."""
    var row_of = List[Int](length=len(graph), fill=-1)
    return align_between(graph, topology, topology.order, SOURCE, SINK, query, scoring, row_of, mode, band_b, band_f)


def align_between(
    graph: Graph,
    topology: Topology,
    nodes: List[Int],
    left: Int,
    right: Int,
    query: Span[UInt8, _],
    scoring: Scoring,
    mut row_of: List[Int],
    mode: Int = GLOBAL,
    band_b: Int = 10,
    band_f: Float64 = 0.01,
) raises -> GraphAlignment:
    """`align_to_graph` over `nodes` alone (in topological order): those strictly between node `left` and node
    `right` (`SOURCE` and `SINK` for the graph's ends), the alignment leaving `left` and, globally, reaching
    `right`. `row_of` is scratch, as long as the graph, every entry -1, and left so."""
    var n = len(nodes)
    for r in range(n):
        row_of[nodes[r]] = r + 1
    var m = len(query)
    var width = m + 1
    var banded = band_b >= 0 and mode == GLOBAL
    var w = band_b + Int(band_f * Float64(m)) if banded else width
    var convex = scoring.convex()
    var gap1 = scoring.open1 + scoring.extend1
    var gap2 = scoring.open2 + scoring.extend2
    var e1 = scoring.extend1
    var e2 = scoring.extend2
    var local = mode == LOCAL

    # The query profile: what letter j - 1 scores against A, C, G, T or N, at [code * width + j].
    var profile = List[Int32](length=5 * width + _V, fill=-scoring.mismatch)
    for code in range(4):
        for j in range(1, width):
            if Int(encode_base(query[j - 1])) == code:
                profile[code * width + j] = scoring.match_score

    # Row 0 is the left end before every source; row r is node order[r - 1].
    var lo = List[Int](length=n + 1, fill=0)
    var hi = List[Int](length=n + 1, fill=0)
    var trace_at = List[Int](length=n + 1, fill=0)
    var choice_at = List[Int](length=n + 1, fill=-1)
    var left_best = List[Int](length=n + 1, fill=0)
    var right_best = List[Int](length=n + 1, fill=0)
    var trace = List[UInt8]()
    var choices = List[UInt8]()
    # Live score rows: H, then the two deletion states, each over the row's band.
    var pool = List[List[Int32]]()
    var slot_of = List[Int](length=n + 1, fill=-1)
    var waiting = List[Int](length=n + 1, fill=0)
    var free = List[Int]()
    var aux = List[Int32]()
    var live_bytes = 0
    var peak = 0

    hi[0] = min(m, w) if banded else m
    var base_remaining = topology.remaining[right] + 1 if right >= 0 else 0
    # Each row's predecessor rows, heaviest edge first; row 0, the left end, waits for those it feeds.
    var pred_rows = List[Int]()
    var pred_at = List[Int](length=n + 2, fill=0)
    for r in range(1, n + 1):
        pred_at[r] = len(pred_rows)
        _heaviest_first(graph, nodes[r - 1], left, row_of, pred_rows)
        for k in range(pred_at[r], len(pred_rows)):
            if pred_rows[k] == 0:
                waiting[0] += 1
    pred_at[n + 1] = len(pred_rows)
    # The left end may reach the right one directly, every letter inserted between them: keep its row.
    var direct = n == 0
    if left >= 0 and right >= 0:
        for s in graph.successors[left]:
            direct = direct or s == right
    if direct:
        waiting[0] += 1
    pool.append(List[Int32](length=3 * (hi[0] + 1), fill=_NEG))
    slot_of[0] = 0
    for j in range(hi[0] + 1):
        pool[0][j] = Int32(0) if local else -scoring.gap(j)
    live_bytes += 12 * (hi[0] + 1)
    trace.resize(hi[0] + 1, _SOURCE)

    var preds = List[Int]()
    var best = Int32(0) if local else _NEG
    var best_row = 0
    var best_col = 0
    var cells = hi[0] + 1
    for r in range(1, n + 1):
        var node = nodes[r - 1]
        preds.clear()
        for k in range(pred_at[r], pred_at[r + 1]):
            preds.append(pred_rows[k])
        var ends_here = False
        waiting[r] = 0
        for s in graph.successors[node]:
            if row_of[s] >= 0:
                waiting[r] += 1
            ends_here = ends_here or s == right
        if right == SINK:
            ends_here = len(graph.successors[node]) == 0
        var low = width
        var high = -1
        var left = width
        var right = 0
        for p in preds:
            low = min(low, lo[p])
            high = max(high, hi[p] + 1)
            left = min(left, left_best[p] + 1)
            right = max(right, right_best[p] + 1)
        high = min(high, m)
        if banded:
            var diagonal = m - max(0, topology.remaining[node] - base_remaining)
            low = max(low, min(left, diagonal) - w)
            high = min(high, max(right, diagonal) + w)
        low = min(low, high)
        lo[r] = low
        hi[r] = high
        var span = high - low + 1
        cells += span

        var slot: Int
        if len(free) > 0:
            slot = free.pop()
        else:
            slot = len(pool)
            pool.append(List[Int32]())
        live_bytes -= 4 * pool[slot].capacity()
        pool[slot].resize(unsafe_uninit_length=3 * span)
        live_bytes += 4 * pool[slot].capacity()
        slot_of[r] = slot
        var row = pool[slot].unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        _fill(row, 3 * span)
        var multi = len(preds) > 1
        aux.resize(unsafe_uninit_length=5 * span)
        var flags = aux.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        _zero(flags, 5 * span)
        var h = row
        var d1 = row.unsafe_offset(span)
        var d2 = row.unsafe_offset(2 * span)
        var x1 = flags
        var x2 = flags.unsafe_offset(span)
        var pick_h = flags.unsafe_offset(2 * span)
        var pick_d1 = flags.unsafe_offset(3 * span)
        var pick_d2 = flags.unsafe_offset(4 * span)
        var letters = profile.unsafe_ptr().unsafe_offset(Int(encode_base(graph.bases[node])) * width)
        var one = _I(1)
        var nothing = _I(0)
        for k in range(len(preds)):
            var p = preds[k]
            var pspan = hi[p] - lo[p] + 1
            var prow = pool[slot_of[p]].unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
            var ph = prow
            var pd1 = prow.unsafe_offset(pspan)
            var pd2 = prow.unsafe_offset(2 * pspan)
            var which = _I(Int32(k))
            # Deletions: this node skipped against query position j.
            var j = max(low, lo[p])
            var last = min(high, hi[p])
            while j <= last:
                var at = j - low
                var pat = j - lo[p]
                if j + _V <= last + 1:
                    var from_h = ph.unsafe_load[width=_V](pat)
                    var open1 = from_h - gap1
                    var extend1 = pd1.unsafe_load[width=_V](pat) - e1
                    var cand1 = max(open1, extend1)
                    var cur1 = d1.unsafe_load[width=_V](at)
                    var better1 = cand1.gt(cur1)
                    d1.unsafe_store(at, better1.select(cand1, cur1))
                    x1.unsafe_store(
                        at, better1.select(extend1.gt(open1).select(one, nothing), x1.unsafe_load[width=_V](at))
                    )
                    if multi:
                        pick_d1.unsafe_store(at, better1.select(which, pick_d1.unsafe_load[width=_V](at)))
                    if convex:
                        var open2 = from_h - gap2
                        var extend2 = pd2.unsafe_load[width=_V](pat) - e2
                        var cand2 = max(open2, extend2)
                        var cur2 = d2.unsafe_load[width=_V](at)
                        var better2 = cand2.gt(cur2)
                        d2.unsafe_store(at, better2.select(cand2, cur2))
                        x2.unsafe_store(
                            at, better2.select(extend2.gt(open2).select(one, nothing), x2.unsafe_load[width=_V](at))
                        )
                        if multi:
                            pick_d2.unsafe_store(at, better2.select(which, pick_d2.unsafe_load[width=_V](at)))
                    j += _V
                    continue
                var from_h = ph[unsafe_offset=pat]
                var open1 = from_h - gap1
                var extend1 = pd1[unsafe_offset=pat] - e1
                var cand1 = max(open1, extend1)
                if cand1 > d1[unsafe_offset=at]:
                    d1[unsafe_offset=at] = cand1
                    x1[unsafe_offset=at] = Int32(1) if extend1 > open1 else Int32(0)
                    pick_d1[unsafe_offset=at] = Int32(k)
                if convex:
                    var open2 = from_h - gap2
                    var extend2 = pd2[unsafe_offset=pat] - e2
                    var cand2 = max(open2, extend2)
                    if cand2 > d2[unsafe_offset=at]:
                        d2[unsafe_offset=at] = cand2
                        x2[unsafe_offset=at] = Int32(1) if extend2 > open2 else Int32(0)
                        pick_d2[unsafe_offset=at] = Int32(k)
                j += 1
            # Matches and mismatches: query letter j - 1 on this node.
            j = max(low, lo[p] + 1, 1)
            last = min(high, hi[p] + 1)
            while j <= last:
                var at = j - low
                var pat = j - 1 - lo[p]
                if j + _V <= last + 1:
                    var cand = ph.unsafe_load[width=_V](pat) + letters.unsafe_load[width=_V](j)
                    var cur = h.unsafe_load[width=_V](at)
                    var better = cand.gt(cur)
                    h.unsafe_store(at, better.select(cand, cur))
                    if multi:
                        pick_h.unsafe_store(at, better.select(which, pick_h.unsafe_load[width=_V](at)))
                    j += _V
                    continue
                var cand = ph[unsafe_offset=pat] + letters[unsafe_offset=j]
                if cand > h[unsafe_offset=at]:
                    h[unsafe_offset=at] = cand
                    pick_h[unsafe_offset=at] = Int32(k)
                j += 1

        # Gaps along the query, left to right: each cell's score, its traceback byte, the row's best cells.
        var tbase = len(trace)
        trace_at[r] = tbase
        trace.resize(unsafe_uninit_length=tbase + span)
        var bytes = trace.unsafe_ptr().unsafe_offset(tbase)
        var row_best = _NEG
        var row_left = low
        var row_right = low
        var run1 = _NEG
        var run2 = _NEG
        var previous = _NEG
        for at in range(span):
            var value = h[unsafe_offset=at]
            var source = _DIAGONAL
            var deleted1 = d1[unsafe_offset=at]
            if deleted1 > value:
                value = deleted1
                source = _DELETION1
            if convex:
                var deleted2 = d2[unsafe_offset=at]
                if deleted2 > value:
                    value = deleted2
                    source = _DELETION2
            var byte = UInt8(0)
            if x1[unsafe_offset=at] != 0:
                byte |= _D1_EXTENDS
            if x2[unsafe_offset=at] != 0:
                byte |= _D2_EXTENDS
            if at > 0:
                var open1 = previous - gap1
                var extend1 = run1 - e1
                run1 = max(open1, extend1)
                if extend1 > open1:
                    byte |= _F1_EXTENDS
                if run1 > value:
                    value = run1
                    source = _INSERTION1
                if convex:
                    var open2 = previous - gap2
                    var extend2 = run2 - e2
                    run2 = max(open2, extend2)
                    if extend2 > open2:
                        byte |= _F2_EXTENDS
                    if run2 > value:
                        value = run2
                        source = _INSERTION2
            if local and value < 0:
                value = 0
                source = _START
            h[unsafe_offset=at] = value
            bytes[unsafe_offset=at] = byte | source
            previous = value
            if value > row_best:
                row_best = value
                row_left = low + at
                row_right = low + at
            elif value == row_best:
                row_right = low + at
        left_best[r] = row_left
        right_best[r] = row_right
        if multi:
            var cbase = len(choices)
            choice_at[r] = cbase
            choices.resize(unsafe_uninit_length=cbase + 3 * span)
            var c = choices.unsafe_ptr().unsafe_offset(cbase)
            for at in range(span):
                c[unsafe_offset=3 * at] = UInt8(pick_h[unsafe_offset=at])
                c[unsafe_offset=3 * at + 1] = UInt8(pick_d1[unsafe_offset=at])
                c[unsafe_offset=3 * at + 2] = UInt8(pick_d2[unsafe_offset=at])
        if mode != GLOBAL:
            if row_best > best:
                best = row_best
                best_row = r
                best_col = row_left
        elif ends_here and high == m and h[unsafe_offset=m - low] > best:
            best = h[unsafe_offset=m - low]
            best_row = r
            best_col = m
        peak = max(peak, live_bytes + trace.capacity() + choices.capacity() + 4 * aux.capacity())
        # Rows no successor still needs give their memory back.
        for p in preds:
            waiting[p] -= 1
            if waiting[p] == 0:
                free.append(slot_of[p])
                slot_of[p] = -1
        if waiting[r] == 0:
            free.append(slot)
            slot_of[r] = -1

    if mode == GLOBAL:
        if direct and slot_of[0] >= 0 and hi[0] == m and pool[slot_of[0]][m] > best:
            best = pool[slot_of[0]][m]
            best_row = 0
            best_col = m
        if best == _NEG:
            for r in range(n):
                row_of[nodes[r]] = -1
            if banded and w < width:
                return align_between(
                    graph, topology, nodes, left, right, query, scoring, row_of, mode, 2 * band_b + 10, 2 * band_f
                )
            raise Error("no alignment reaches the graph's end")

    var aligned = List[Int](length=m, fill=-1)
    var r = best_row
    var j = best_col
    var state = UInt8(0)  # 0: H; otherwise the gap state, as the source codes number them
    while r > 0 and j > 0:
        var at = j - lo[r]
        var byte = trace[trace_at[r] + at]
        var node = nodes[r - 1]
        if state == 0:
            var source = byte & 7
            if source == _START:
                break
            if source == _DIAGONAL:
                aligned[j - 1] = node
                r = _predecessor(pred_rows, pred_at, choices, choice_at, r, at, 0)
                j -= 1
            else:
                state = source
        elif state == _DELETION1 or state == _DELETION2:
            var first = state == _DELETION1
            var extends = (byte & (_D1_EXTENDS if first else _D2_EXTENDS)) != 0
            r = _predecessor(pred_rows, pred_at, choices, choice_at, r, at, 1 if first else 2)
            if not extends:
                state = 0
        else:
            var extends = (byte & (_F1_EXTENDS if state == _INSERTION1 else _F2_EXTENDS)) != 0
            j -= 1
            if not extends:
                state = 0
    for k in range(n):
        row_of[nodes[k]] = -1
    return GraphAlignment(aligned^, Int(best), cells, peak)


@inline(.always)
def _predecessor(
    pred_rows: List[Int], pred_at: List[Int], choices: List[UInt8], choice_at: List[Int], r: Int, cell: Int, kind: Int
) -> Int:
    """The row a move at `cell` of row `r` came from (`kind`: 0 the diagonal, 1 and 2 the deletions)."""
    var which = 0 if choice_at[r] < 0 else Int(choices[choice_at[r] + 3 * cell + kind])
    return pred_rows[pred_at[r] + which]


def _heaviest_first(graph: Graph, node: Int, left: Int, row_of: List[Int], mut rows: List[Int]):
    """Appends the rows of `node`'s predecessors in the alignment, heaviest edge first (the earlier edge first
    among equals): row 0 for `left` (or for none at all), so of equally good moves the alignment takes the one
    most sequences took."""
    ref preds = graph.predecessors[node]
    var keys = List[UInt64](capacity=len(preds))
    for k in range(len(preds)):
        var p = preds[k]
        if p != left and row_of[p] < 0:
            continue
        var weight = 0
        for i in range(len(graph.successors[p])):
            if graph.successors[p][i] == node:
                weight = graph.weights[p][i]
        keys.append((UInt64((1 << 31) - weight) << 32) | UInt64(k))
    if len(keys) == 0:
        rows.append(0)
        return
    sort(keys)
    for key in keys:
        var p = preds[Int(key & 0xFFFFFFFF)]
        rows.append(0 if p == left else row_of[p])


@inline(.always)
def _zero(p: MutPointer[Int32, _], count: Int):
    var i = 0
    while i + _V <= count:
        p.unsafe_store(i, _I(0))
        i += _V
    while i < count:
        p[unsafe_offset=i] = 0
        i += 1
