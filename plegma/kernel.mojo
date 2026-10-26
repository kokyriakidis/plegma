"""Aligning a sequence to the whole graph: banded dynamic programming over the graph's nodes in topological order.

Scores are maximized, as abPOA's: a match earns `match`, a mismatch costs `mismatch`, a gap of `k` letters costs
the less of `open1 + extend1 k` and `open2 + extend2 k` (convex; `open2` 0 makes it affine, `open1` 0 linear).
Each node's row covers only its band of query positions: abPOA's adaptive band, `b + f m` letters either side of
the predecessors' best cells and of the diagonal the node's remaining path (along the heaviest edges) leaves.

Memory follows the band and the graph's width, not the graph times the query. A row's scores live only until its
last successor is computed; what the traceback needs is one byte per cell (where each cell's score came from), and
at nodes with several predecessors three more (which predecessor each move came from).
"""

from std.math import iota
from std.sys import simd_width_of, size_of

from .alphabet import encode_base
from .graph import Graph


comptime SOURCE = -1
"""The left end before every path's first node."""
comptime SINK = -2
"""The right end after every path's last node."""
comptime GLOBAL = 0
"""The whole query against a path from a graph source to a graph sink."""
comptime LOCAL = 1
"""The best-scoring part of the query against any part of a path."""
comptime EXTENSION = 2
"""From the graph's sources and the query's start, as far as it pays."""
comptime PREFIX = 3
"""A read anchored at the graph's start: from a source and the query's start, every query letter, ending at any node;
past a sink, the query's remaining letters are free (an overhang the graph does not hold yet)."""
comptime SUFFIX = 4
"""A read anchored at the graph's end: from any node (the query's letters before a source free: an overhang), every
query letter, ending at a sink."""

comptime BOTH = 0
"""A read spanning both anchors of a local assembly (aligned `GLOBAL`)."""
comptime LEFT = 1
"""A read fixed at the left anchor only (aligned `PREFIX`)."""
comptime RIGHT = 2
"""A read fixed at the right anchor only (aligned `SUFFIX`)."""


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
def _neg[dt: DType]() -> Scalar[dt]:
    """The score of a cell nothing reaches: far below every real score, far enough above the type's floor that
    a few gap costs taken from it cannot wrap (the kernel also clamps at it)."""
    comptime if dt == DType.int16:
        return Scalar[dt](-16000)
    else:
        return Scalar[dt](-(1 << 29))


@inline(.always)
def _banded(mode: Int) -> Bool:
    """Whether `mode` runs in the adaptive band: those that use every query letter."""
    return mode == GLOBAL or mode == PREFIX or mode == SUFFIX


def _fits_16_bits(scoring: Scoring, m: Int, n: Int) -> Bool:
    """Whether every score an alignment of `m` letters to `n` nodes can reach fits 16-bit lanes: a cell is
    never below deleting every node and inserting every letter, never above matching every letter, and the gap
    scan adds at most the extension cost of the whole query to that. The floor stays above half the sentinel,
    below which a score counts as unreached."""
    var low = -Int(scoring.gap(n)) - Int(scoring.gap(m))
    var high = Int(scoring.match_score) * m + Int(max(scoring.extend1, scoring.extend2)) * (m + 1)
    return low > Int(_neg[DType.int16]()) // 2 + 500 and high < 30000


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
    the heaviest out-edges (abPOA's band diagonal); or the same over a part of the graph (`allowed`)."""

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
            var best_weight = Int32(-1)
            var e = Int(graph.first_out[node])
            while e >= 0:
                var s = Int(graph.edge_to[e])
                var weight = graph.edge_weight[e]
                if weight > best_weight or (weight == best_weight and self.remaining[s] + 1 > self.remaining[node]):
                    best_weight = weight
                    self.remaining[node] = self.remaining[s] + 1
                e = Int(graph.next_out[e])

    def __init__(out self, graph: Graph, allowed: List[Bool]) raises:
        """Over the nodes `allowed` alone: they in topological order (`order`), and each one's remaining path
        through them alone; every node's rank in the whole graph."""
        self.rank = graph.node_order()
        var n = len(self.rank)
        var full = List[Int](length=n, fill=0)
        for node in range(n):
            full[self.rank[node]] = node
        self.order = List[Int]()
        for node in full:
            if allowed[node]:
                self.order.append(node)
        self.remaining = List[Int](length=n, fill=0)
        for r in reversed(range(len(self.order))):
            var node = self.order[r]
            var best_weight = Int32(-1)
            var e = Int(graph.first_out[node])
            while e >= 0:
                var s = Int(graph.edge_to[e])
                var weight = graph.edge_weight[e]
                if allowed[s] and (
                    weight > best_weight or (weight == best_weight and self.remaining[s] + 1 > self.remaining[node])
                ):
                    best_weight = weight
                    self.remaining[node] = self.remaining[s] + 1
                e = Int(graph.next_out[e])


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


struct Workspace(Movable):
    """What an alignment needs beside the graph, kept from one alignment to the next so that a run of them
    allocates once: per-row bands and traceback offsets, the live score rows, the traceback itself."""

    var profile: List[UInt8]
    var lo: List[Int]
    var hi: List[Int]
    var trace_at: List[Int]
    var choice_at: List[Int]
    var left_best: List[Int]
    var right_best: List[Int]
    var trace: List[UInt8]
    var choices: List[UInt8]
    var pool: List[List[UInt8]]
    var slot_of: List[Int]
    var waiting: List[Int]
    var free: List[Int]
    var aux: List[UInt8]
    var pred_rows: List[Int]
    var pred_at: List[Int]
    var preds: List[Int]
    var keys: List[UInt64]

    def __init__(out self):
        self.profile = List[UInt8]()
        self.lo = List[Int]()
        self.hi = List[Int]()
        self.trace_at = List[Int]()
        self.choice_at = List[Int]()
        self.left_best = List[Int]()
        self.right_best = List[Int]()
        self.trace = List[UInt8]()
        self.choices = List[UInt8]()
        self.pool = List[List[UInt8]]()
        self.slot_of = List[Int]()
        self.waiting = List[Int]()
        self.free = List[Int]()
        self.aux = List[UInt8]()
        self.pred_rows = List[Int]()
        self.pred_at = List[Int]()
        self.preds = List[Int]()
        self.keys = List[UInt64]()

    def align(
        mut self,
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
        `right`. `row_of` is scratch, as long as the graph, every entry -1, and left so. In 16-bit lanes (twice
        as many a vector) when every score the alignment can reach fits them (on shasta2's local assemblies the
        same alignments as 32-bit lanes and 2% faster, narrow rows too)."""
        var m = len(query)
        var reach = band_b + Int(band_f * Float64(m)) if band_b >= 0 and _banded(mode) else m + 1
        if _fits_16_bits(scoring, m, len(nodes)):
            return self._align[DType.int16](
                graph, topology, nodes, left, right, query, scoring, row_of, mode, band_b, band_f
            )
        return self._align[DType.int32](
            graph, topology, nodes, left, right, query, scoring, row_of, mode, band_b, band_f
        )

    def _align[
        dt: DType
    ](
        mut self,
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
        comptime V = simd_width_of[dt]()
        comptime S = Scalar[dt]
        comptime Vec = SIMD[dt, V]
        comptime NEG = _neg[dt]()
        comptime BYTES = size_of[S]()
        var n = len(nodes)
        for r in range(n):
            row_of[nodes[r]] = r + 1
        var m = len(query)
        var width = m + 1
        var banded = band_b >= 0 and _banded(mode)
        var w = band_b + Int(band_f * Float64(m)) if banded else width
        var convex = scoring.convex()
        var gap1 = S(Int(scoring.open1 + scoring.extend1))
        var gap2 = S(Int(scoring.open2 + scoring.extend2))
        var e1 = S(Int(scoring.extend1))
        var e2 = S(Int(scoring.extend2))
        var local = mode == LOCAL
        var free_start = mode == SUFFIX

        # The workspace's arrays, sized to this alignment; whatever an earlier one left is overwritten before use.
        self.profile.resize(BYTES * (5 * width + V), 0)
        var profile = self.profile.unsafe_ptr().unsafe_bitcast[S]().unsafe_origin_cast[MutUntrackedOrigin]()
        var k = 0
        while k + V <= 5 * width + V:
            profile.unsafe_store(k, Vec(S(-Int(scoring.mismatch))))
            k += V
        while k < 5 * width + V:
            profile[unsafe_offset=k] = S(-Int(scoring.mismatch))
            k += 1
        self.lo.resize(n + 1, 0)
        self.hi.resize(n + 1, 0)
        self.trace_at.resize(n + 1, 0)
        self.choice_at.resize(n + 1, -1)
        self.left_best.resize(n + 1, 0)
        self.right_best.resize(n + 1, 0)
        self.left_best[0] = 0
        self.right_best[0] = 0
        self.lo[0] = 0
        self.slot_of.resize(n + 1, -1)
        self.waiting.resize(n + 1, 0)
        self.waiting[0] = 0
        self.trace.clear()
        self.choices.clear()
        self.pred_rows.clear()
        self.pred_at.resize(n + 2, 0)
        self.free.clear()
        for slot in reversed(range(len(self.pool))):
            self.free.append(slot)

        # The query self.profile: what letter j - 1 scores against A, C, G, T or N, at [code * width + j].
        for j in range(1, width):
            var code = Int(encode_base(query[j - 1]))
            if code < 4:
                profile[unsafe_offset=code * width + j] = S(Int(scoring.match_score))

        # Row 0 is the left end before every source; row r is node order[r - 1].
        # Live score rows: H, then the two deletion states, each over the row's band.
        var live_bytes = 0
        for slot in range(len(self.pool)):
            live_bytes += self.pool[slot].capacity()
        var peak = 0

        # Row 0, before every source: a read anchored at the graph's end alone may begin with any number of letters
        # the graph does not hold (free), so all of it there; a single row.
        self.hi[0] = min(m, w) if banded and not free_start else m
        var base_remaining = topology.remaining[right] + 1 if right >= 0 else 0
        # Each row's predecessor rows, heaviest edge first; row 0, the left end, waits for those it feeds.
        for r in range(1, n + 1):
            self.pred_at[r] = len(self.pred_rows)
            _heaviest_first(graph, nodes[r - 1], left, row_of, self.pred_rows, self.keys)
            for k in range(self.pred_at[r], len(self.pred_rows)):
                if self.pred_rows[k] == 0:
                    self.waiting[0] += 1
        self.pred_at[n + 1] = len(self.pred_rows)
        # The left end may reach the right one directly, every letter inserted between them: keep its row.
        var direct = n == 0
        if left >= 0 and right >= 0:
            direct = direct or graph.edge(left, right) >= 0
        if direct:
            self.waiting[0] += 1
        var first_slot: Int
        if len(self.free) > 0:
            first_slot = self.free.pop()
        else:
            first_slot = len(self.pool)
            self.pool.append(List[UInt8]())
        var first_stride = self.hi[0] + 1 + 2 * V
        live_bytes -= self.pool[first_slot].capacity()
        self.pool[first_slot].resize(unsafe_uninit_length=BYTES * 3 * first_stride)
        live_bytes += self.pool[first_slot].capacity()
        self.slot_of[0] = first_slot
        var first_row = self.pool[first_slot].unsafe_ptr().unsafe_bitcast[S]().unsafe_origin_cast[MutUntrackedOrigin]()
        _fill[dt](first_row, 3 * first_stride)
        for j in range(self.hi[0] + 1):
            first_row[unsafe_offset=V + j] = S(0) if local or free_start else max(S(-Int(scoring.gap(j))), NEG)
        self.trace.resize(self.hi[0] + 1, _SOURCE)

        var best = S(0) if local else NEG
        var best_row = 0
        var best_col = 0
        var cells = self.hi[0] + 1
        for r in range(1, n + 1):
            var node = nodes[r - 1]
            self.choice_at[r] = -1
            self.preds.clear()
            for k in range(self.pred_at[r], self.pred_at[r + 1]):
                self.preds.append(self.pred_rows[k])
            var ends_here = False
            self.waiting[r] = 0
            var e = Int(graph.first_out[node])
            while e >= 0:
                var s = Int(graph.edge_to[e])
                if row_of[s] >= 0:
                    self.waiting[r] += 1
                ends_here = ends_here or s == right
                e = Int(graph.next_out[e])
            if right == SINK:
                # A sink of the nodes aligned to: none of its successors among them (the graph's own sinks when
                # aligning to all of it).
                ends_here = self.waiting[r] == 0
            var low = width
            var high = -1
            var left = width
            var right = 0
            for p in self.preds:
                low = min(low, self.lo[p])
                high = max(high, self.hi[p] + 1)
                left = min(left, self.left_best[p] + 1)
                right = max(right, self.right_best[p] + 1)
            high = min(high, m)
            # Only cells some predecessor reaches: [low, high] as they stand.
            var reach_low = low
            var reach_high = high
            if banded and mode == PREFIX:
                # No end to reach: the band follows the best cells.
                low = max(low, left - w)
                high = min(high, right + w)
            elif banded and mode == SUFFIX:
                # Anchored at the end alone: the diagonal the node's remaining path leaves, nothing else. With a free
                # start the best cells are chance matches of the read's beginning anywhere before its place, and
                # following them spread the band over the whole read (a 14 kb read against 37,000 nodes: 290 million
                # cells, not 11 million).
                var diagonal = m - max(0, topology.remaining[node] - base_remaining)
                low = max(low, diagonal - w)
                high = min(high, diagonal + w)
            elif banded:
                var diagonal = m - max(0, topology.remaining[node] - base_remaining)
                low = max(low, min(left, diagonal) - w)
                high = min(high, max(right, diagonal) + w)
            low = min(max(low, reach_low), reach_high)
            high = max(min(high, reach_high), low)
            self.lo[r] = low
            self.hi[r] = high
            var span = high - low + 1
            cells += span

            var slot: Int
            if len(self.free) > 0:
                slot = self.free.pop()
            else:
                slot = len(self.pool)
                self.pool.append(List[UInt8]())
            # A row is H, D1, D2, each over the band with `V` sentinel cells either side, so loads may run a
            # cell past a predecessor's band.
            var stride = span + 2 * V
            live_bytes -= self.pool[slot].capacity()
            self.pool[slot].resize(unsafe_uninit_length=BYTES * 3 * stride)
            live_bytes += self.pool[slot].capacity()
            self.slot_of[r] = slot
            var row = self.pool[slot].unsafe_ptr().unsafe_bitcast[S]().unsafe_origin_cast[MutUntrackedOrigin]()
            var multi = len(self.preds) > 1
            if multi:
                _fill[dt](row, 3 * stride)
            else:
                for third in range(3):
                    _fill[dt](row.unsafe_offset(third * stride), V)
                    _fill[dt](row.unsafe_offset(third * stride + V + span), V)
            var lanes = span + V  # a vector of slack: the last vector of a row runs whole
            self.aux.resize(unsafe_uninit_length=BYTES * 6 * lanes)
            var flags = self.aux.unsafe_ptr().unsafe_bitcast[S]().unsafe_origin_cast[MutUntrackedOrigin]()
            var h = row.unsafe_offset(V)
            var d1 = row.unsafe_offset(stride + V)
            var d2 = row.unsafe_offset(2 * stride + V)
            var x1 = flags
            var x2 = flags.unsafe_offset(lanes)
            var pick_h = flags.unsafe_offset(2 * lanes)
            var pick_d1 = flags.unsafe_offset(3 * lanes)
            var pick_d2 = flags.unsafe_offset(4 * lanes)
            var source_of = flags.unsafe_offset(5 * lanes)
            var letters = profile.unsafe_offset(Int(encode_base(graph.bases[node])) * width)
            var one = Vec(1)
            var nothing = Vec(0)
            if multi:
                _zero[dt](flags, 5 * lanes)
                for k in range(len(self.preds)):
                    var p = self.preds[k]
                    var pstride = self.hi[p] - self.lo[p] + 1 + 2 * V
                    var prow = (
                        self.pool[self.slot_of[p]]
                        .unsafe_ptr()
                        .unsafe_bitcast[S]()
                        .unsafe_origin_cast[MutUntrackedOrigin]()
                    )
                    var ph = prow.unsafe_offset(V)
                    var pd1 = prow.unsafe_offset(pstride + V)
                    var pd2 = prow.unsafe_offset(2 * pstride + V)
                    var which = Vec(S(k))
                    # Deletions: this node skipped against query position j.
                    var j = max(low, self.lo[p])
                    var last = min(high, self.hi[p])
                    while j <= last:
                        var at = j - low
                        var pat = j - self.lo[p]
                        if j + V <= last + 1:
                            var from_h = ph.unsafe_load[width=V](pat)
                            var open1 = from_h - gap1
                            var extend1 = pd1.unsafe_load[width=V](pat) - e1
                            var cand1 = _floor(max(open1, extend1))
                            var cur1 = d1.unsafe_load[width=V](at)
                            var better1 = cand1.gt(cur1)
                            d1.unsafe_store(at, better1.select(cand1, cur1))
                            x1.unsafe_store(
                                at, better1.select(extend1.gt(open1).select(one, nothing), x1.unsafe_load[width=V](at))
                            )
                            if multi:
                                pick_d1.unsafe_store(at, better1.select(which, pick_d1.unsafe_load[width=V](at)))
                            if convex:
                                var open2 = from_h - gap2
                                var extend2 = pd2.unsafe_load[width=V](pat) - e2
                                var cand2 = _floor(max(open2, extend2))
                                var cur2 = d2.unsafe_load[width=V](at)
                                var better2 = cand2.gt(cur2)
                                d2.unsafe_store(at, better2.select(cand2, cur2))
                                x2.unsafe_store(
                                    at,
                                    better2.select(extend2.gt(open2).select(one, nothing), x2.unsafe_load[width=V](at)),
                                )
                                if multi:
                                    pick_d2.unsafe_store(at, better2.select(which, pick_d2.unsafe_load[width=V](at)))
                            j += V
                            continue
                        var from_h = ph[unsafe_offset=pat]
                        var open1 = from_h - gap1
                        var extend1 = pd1[unsafe_offset=pat] - e1
                        var cand1 = _floor(max(open1, extend1))
                        if cand1 > d1[unsafe_offset=at]:
                            d1[unsafe_offset=at] = cand1
                            x1[unsafe_offset=at] = S(1) if extend1 > open1 else S(0)
                            pick_d1[unsafe_offset=at] = S(k)
                        if convex:
                            var open2 = from_h - gap2
                            var extend2 = pd2[unsafe_offset=pat] - e2
                            var cand2 = _floor(max(open2, extend2))
                            if cand2 > d2[unsafe_offset=at]:
                                d2[unsafe_offset=at] = cand2
                                x2[unsafe_offset=at] = S(1) if extend2 > open2 else S(0)
                                pick_d2[unsafe_offset=at] = S(k)
                        j += 1
                    # Matches and mismatches: query letter j - 1 on this node.
                    j = max(low, self.lo[p] + 1, 1)
                    last = min(high, self.hi[p] + 1)
                    while j <= last:
                        var at = j - low
                        var pat = j - 1 - self.lo[p]
                        if j + V <= last + 1:
                            var cand = _floor(ph.unsafe_load[width=V](pat) + letters.unsafe_load[width=V](j))
                            var cur = h.unsafe_load[width=V](at)
                            var better = cand.gt(cur)
                            h.unsafe_store(at, better.select(cand, cur))
                            if multi:
                                pick_h.unsafe_store(at, better.select(which, pick_h.unsafe_load[width=V](at)))
                            j += V
                            continue
                        var cand = _floor(ph[unsafe_offset=pat] + letters[unsafe_offset=j])
                        if cand > h[unsafe_offset=at]:
                            h[unsafe_offset=at] = cand
                            pick_h[unsafe_offset=at] = S(k)
                        j += 1

                # The best of the diagonal and the deletions, and where it came from.
                var at = 0
                while at < span:
                    var value = h.unsafe_load[width=V](at)
                    var source = Vec(0)
                    var deleted1 = d1.unsafe_load[width=V](at)
                    var take1 = deleted1.gt(value)
                    value = take1.select(deleted1, value)
                    source = take1.select(Vec(S(_DELETION1)), source)
                    if convex:
                        var deleted2 = d2.unsafe_load[width=V](at)
                        var take2 = deleted2.gt(value)
                        value = take2.select(deleted2, value)
                        source = take2.select(Vec(S(_DELETION2)), source)
                    h.unsafe_store(at, value)
                    source_of.unsafe_store(
                        at, source | (x1.unsafe_load[width=V](at) << 3) | (x2.unsafe_load[width=V](at) << 4)
                    )
                    at += V

            else:
                # One predecessor: the diagonal, both deletions and the better of them in one pass.
                var p = self.preds[0]
                var pstride = self.hi[p] - self.lo[p] + 1 + 2 * V
                var prow = (
                    self.pool[self.slot_of[p]].unsafe_ptr().unsafe_bitcast[S]().unsafe_origin_cast[MutUntrackedOrigin]()
                )
                var ph = prow.unsafe_offset(V)
                var pd1 = prow.unsafe_offset(pstride + V)
                var pd2 = prow.unsafe_offset(2 * pstride + V)
                var shift = low - self.lo[p]
                var at = 0
                while at < span:
                    var pat = at + shift
                    var from_h = ph.unsafe_load[width=V](pat)
                    var open1 = from_h - gap1
                    var extend1 = pd1.unsafe_load[width=V](pat) - e1
                    var deleted1 = _floor(max(open1, extend1))
                    var value = _floor(ph.unsafe_load[width=V](pat - 1) + letters.unsafe_load[width=V](low + at))
                    var take1 = deleted1.gt(value)
                    var source = take1.select(Vec(S(_DELETION1)), Vec(0))
                    value = take1.select(deleted1, value)
                    d1.unsafe_store(at, deleted1)
                    # The deletions' extend flags packed with the move, as the traceback byte holds them.
                    var flags = extend1.gt(open1).select(one, nothing) << 3
                    if convex:
                        var open2 = from_h - gap2
                        var extend2 = pd2.unsafe_load[width=V](pat) - e2
                        var deleted2 = _floor(max(open2, extend2))
                        var take2 = deleted2.gt(value)
                        source = take2.select(Vec(S(_DELETION2)), source)
                        value = take2.select(deleted2, value)
                        d2.unsafe_store(at, deleted2)
                        flags = flags | (extend2.gt(open2).select(one, nothing) << 4)
                    else:
                        d2.unsafe_store(at, Vec(NEG))
                    h.unsafe_store(at, value)
                    source_of.unsafe_store(at, source | flags)
                    at += V
            if free_start and low == 0 and h[unsafe_offset=0] < S(0):
                # The alignment may start at this node: the graph before it skipped free.
                h[unsafe_offset=0] = S(0)
                source_of[unsafe_offset=0] = (source_of[unsafe_offset=0] & S(~7)) | S(_START)
            # One pass over the row: both gap pieces' insertion scores, then each cell's score, its traceback byte
            # and the row's best. Gaps along the query, left to right, from the scores so far (`h`, before
            # insertions): a gap opened after another never beats one longer gap (with E1 >= E2), so the insertion
            # states need not see each other. As running maxima: F[k] = max over t < k of (H[t] + e t) - (o + e) -
            # e (k - 1), a prefix maximum in lane shifts, carried from block to block; a cell extends when its left
            # neighbour's gap, less e, beats opening from its left neighbour. The two pieces in one loop run side
            # by side (their chains are independent) and share the loads and the left neighbours.
            var tbase = len(self.trace)
            self.trace_at[r] = tbase
            self.trace.resize(unsafe_uninit_length=tbase + span + V)
            var bytes = self.trace.unsafe_ptr().unsafe_offset(tbase)
            var best_lanes = Vec(NEG)
            var lane = iota[dt, V]()
            var first_lane = lane.eq(0)
            var ramp1 = Vec(e1) * lane
            var ramp2 = Vec(e2) * lane
            var opened1 = gap1 - e1
            var opened2 = gap2 - e2
            var carry1 = NEG  # the largest h[t] + e1 t of the cells before the block
            var carry2 = NEG
            var last_f1 = NEG
            var last_f2 = NEG
            var last_h = NEG
            var at = 0
            while at < span:
                var here = h.unsafe_load[width=V](at)
                var left_h = first_lane.select(Vec(last_h), here.shift_right[1]())
                var e_index1 = ramp1 + Vec(e1 * S(at))
                var reach1 = _prefix_max(here + e_index1)
                var before1 = max(first_lane.select(Vec(carry1), reach1.shift_right[1]()), Vec(carry1))
                var inserted1 = _floor(before1 - opened1 - e_index1)
                if at == 0:
                    inserted1 = first_lane.select(Vec(NEG), inserted1)
                var left_f1 = first_lane.select(Vec(last_f1), inserted1.shift_right[1]())
                var extends1 = (left_f1 - e1).gt(left_h - gap1).select(one, nothing)
                carry1 = max(carry1, reach1[V - 1])
                last_f1 = inserted1[V - 1]
                var inserted2 = Vec(NEG)
                var extends2 = nothing
                if convex:
                    var e_index2 = ramp2 + Vec(e2 * S(at))
                    var reach2 = _prefix_max(here + e_index2)
                    var before2 = max(first_lane.select(Vec(carry2), reach2.shift_right[1]()), Vec(carry2))
                    inserted2 = _floor(before2 - opened2 - e_index2)
                    if at == 0:
                        inserted2 = first_lane.select(Vec(NEG), inserted2)
                    var left_f2 = first_lane.select(Vec(last_f2), inserted2.shift_right[1]())
                    extends2 = (left_f2 - e2).gt(left_h - gap2).select(one, nothing)
                    carry2 = max(carry2, reach2[V - 1])
                    last_f2 = inserted2[V - 1]
                last_h = here[V - 1]
                var value = here
                var packed = source_of.unsafe_load[width=V](at)
                var source = packed & Vec(S(7))
                var take1 = inserted1.gt(value)
                value = take1.select(inserted1, value)
                source = take1.select(Vec(S(_INSERTION1)), source)
                var take2 = inserted2.gt(value)
                value = take2.select(inserted2, value)
                source = take2.select(Vec(S(_INSERTION2)), source)
                if local:
                    var below = value.lt(Vec(0))
                    value = below.select(Vec(0), value)
                    source = below.select(Vec(S(_START)), source)
                var byte = source | (packed & Vec(S(~7))) | (extends1 << 5) | (extends2 << 6)
                h.unsafe_store(at, value)
                bytes.unsafe_store(at, byte.cast[DType.uint8]())
                best_lanes = max(best_lanes, (lane + S(at)).lt(S(span)).select(value, Vec(NEG)))
                at += V
            var row_best = best_lanes.reduce_max()
            # The last vector ran into the right padding: back to sentinels for the successors.
            _fill[dt](h.unsafe_offset(span), V)
            _fill[dt](d1.unsafe_offset(span), V)
            _fill[dt](d2.unsafe_offset(span), V)
            self.trace.resize(unsafe_uninit_length=tbase + span)
            var row_left = low + _first_equal[dt](h, span, row_best)
            var row_right = low + _last_equal[dt](h, span, row_best)
            self.left_best[r] = row_left
            self.right_best[r] = row_right
            if multi:
                var cbase = len(self.choices)
                self.choice_at[r] = cbase
                self.choices.resize(unsafe_uninit_length=cbase + 3 * span)
                var c = self.choices.unsafe_ptr().unsafe_offset(cbase)
                for at in range(span):
                    c[unsafe_offset=3 * at] = UInt8(pick_h[unsafe_offset=at])
                    c[unsafe_offset=3 * at + 1] = UInt8(pick_d1[unsafe_offset=at])
                    c[unsafe_offset=3 * at + 2] = UInt8(pick_d2[unsafe_offset=at])
            if mode == LOCAL or mode == EXTENSION:
                if row_best > best:
                    best = row_best
                    best_row = r
                    best_col = row_left
            elif mode == PREFIX:
                # Every query letter used, at any node; at a sink, the letters left free.
                if high == m and h[unsafe_offset=m - low] > best:
                    best = h[unsafe_offset=m - low]
                    best_row = r
                    best_col = m
                if ends_here and row_best > best:
                    best = row_best
                    best_row = r
                    best_col = row_right
            elif ends_here and high == m and h[unsafe_offset=m - low] > best:
                best = h[unsafe_offset=m - low]
                best_row = r
                best_col = m
            peak = max(peak, live_bytes + self.trace.capacity() + self.choices.capacity() + self.aux.capacity())
            # Rows no successor still needs give their memory back.
            for p in self.preds:
                self.waiting[p] -= 1
                if self.waiting[p] == 0:
                    self.free.append(self.slot_of[p])
                    self.slot_of[p] = -1
            if self.waiting[r] == 0:
                self.free.append(slot)
                self.slot_of[r] = -1

        if mode == GLOBAL or mode == SUFFIX:
            var start_row = self.pool[self.slot_of[0] if self.slot_of[0] >= 0 else 0].unsafe_ptr().unsafe_bitcast[S]()
            if direct and self.slot_of[0] >= 0 and self.hi[0] == m and start_row[unsafe_offset=V + m] > best:
                best = start_row[unsafe_offset=V + m]
                best_row = 0
                best_col = m
            # Scores near the sentinel are cells nothing reaches.
            if best < NEG // 2:
                for r in range(n):
                    row_of[nodes[r]] = -1
                if banded and w < width:
                    return self.align(
                        graph, topology, nodes, left, right, query, scoring, row_of, mode, 2 * band_b + 10, 2 * band_f
                    )
                raise Error("no alignment reaches the graph's end")

        var aligned = List[Int](length=m, fill=-1)
        var r = best_row
        var j = best_col
        var state = UInt8(0)  # 0: H; otherwise the gap state, as the source codes number them
        while r > 0 and j > 0:
            var at = j - self.lo[r]
            var byte = self.trace[self.trace_at[r] + at]
            var node = nodes[r - 1]
            if state == 0:
                var source = byte & 7
                if source == _START:
                    break
                if source == _DIAGONAL:
                    aligned[j - 1] = node
                    r = _predecessor(self.pred_rows, self.pred_at, self.choices, self.choice_at, r, at, 0)
                    j -= 1
                else:
                    state = source
            elif state == _DELETION1 or state == _DELETION2:
                var first = state == _DELETION1
                var extends = (byte & (_D1_EXTENDS if first else _D2_EXTENDS)) != 0
                r = _predecessor(self.pred_rows, self.pred_at, self.choices, self.choice_at, r, at, 1 if first else 2)
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
    """`Workspace.align` through a workspace of its own."""
    var workspace = Workspace()
    return workspace.align(graph, topology, nodes, left, right, query, scoring, row_of, mode, band_b, band_f)


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
    if len(keys) > 1:
        # Equal weights keep the edges' order: edge ids grow as the edges are made.
        sort(keys)
    for key in keys:
        var p = Int(graph.edge_from[Int(key & 0xFFFFFFFF)])
        rows.append(0 if p == left else row_of[p])


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
    var lane = iota[dt, width]()
    var out = v
    comptime for step in range(6):
        comptime shift = 1 << step
        comptime if shift < width:
            out = max(out, lane.lt(Scalar[dt](shift)).select(SIMD[dt, width](_neg[dt]()), out.shift_right[shift]()))
    return out


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
