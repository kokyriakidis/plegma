"""A sequence aligned to the graph through pins: minimap2's pins (`pins.pins`) between the sequence and the graph's
consensus path, each pin's letters matched to its path node, and the kernel's DP run only in the windows between two
pinned nodes (the subgraph a path from one to the other can visit). The DP's cells then follow the windows, not the
whole graph times the sequence; where the chain finds no pins, the one window is the whole graph."""

from std.atomic import Atomic

from karan.common import spread, thread_count

from .alphabet import encode_base
from .graph import Graph
from .kernel import GLOBAL, PREFIX, SINK, SOURCE, SUFFIX, GraphAlignment, Scoring, Topology, Workspace
from .pins import Pins, pins_with_index
from .seeds import POSITION_BITS, ReferenceIndex
from .trace import PINS, Trace


struct PinAnchors(Copyable, Movable):
    """Precomputed one-letter pin anchors between a consensus path and a query sequence."""

    var anchor_r: List[Int]
    var anchor_q: List[Int]
    var raw_count: Int
    var valid: Bool

    def __init__(out self):
        self.anchor_r = List[Int]()
        self.anchor_q = List[Int]()
        self.raw_count = 0
        self.valid = False

    def __init__(out self, var anchor_r: List[Int], var anchor_q: List[Int], raw_count: Int):
        self.anchor_r = anchor_r^
        self.anchor_q = anchor_q^
        self.raw_count = raw_count
        self.valid = True


struct _PinnedWorker(Movable):
    """Per-worker reusable state for parallel `align_pinned` window DP and snapshot batch alignment."""

    var row_of: List[Int]
    var marks: List[UInt8]
    var touched: List[Int]
    var stack: List[Int]
    var keys: List[Int]
    var workspace: Workspace

    def __init__(out self):
        self.row_of = List[Int]()
        self.marks = List[UInt8]()
        self.touched = List[Int]()
        self.stack = List[Int]()
        self.keys = List[Int]()
        self.workspace = Workspace()


struct WindowScratch(Movable):
    """Reusable marks for `window_nodes`, consensus-path `ReferenceIndex` cache, speculative precomputed pin
    anchors, heaviest-path buffers, and parallel worker pool for `align_pinned`."""

    var marks: List[UInt8]
    var touched: List[Int]
    var stack: List[Int]
    var keys: List[Int]
    var hp_starts: List[Int]
    var hp_ends: List[Int]
    var hp_score: List[Int]
    var hp_next: List[Int]
    var hp_path: List[Int]
    var hp_letters: List[UInt8]
    var cached_letters: List[UInt8]
    var cached_index: ReferenceIndex
    var has_cached_index: Bool
    var precomp_letters: List[UInt8]
    var precomp_anchors: List[PinAnchors]
    var workers: List[_PinnedWorker]

    def __init__(out self):
        self.marks = List[UInt8]()
        self.touched = List[Int]()
        self.stack = List[Int]()
        self.keys = List[Int]()
        self.hp_starts = List[Int]()
        self.hp_ends = List[Int]()
        self.hp_score = List[Int]()
        self.hp_next = List[Int]()
        self.hp_path = List[Int]()
        self.hp_letters = List[UInt8]()
        self.cached_letters = List[UInt8]()
        self.cached_index = ReferenceIndex()
        self.has_cached_index = False
        self.precomp_letters = List[UInt8]()
        self.precomp_anchors = List[PinAnchors]()
        self.workers = List[_PinnedWorker]()

    def clear_precomp(mut self):
        self.hp_letters.clear()
        self.precomp_letters.clear()
        self.precomp_anchors.clear()


def _window_nodes_into(
    graph: Graph,
    rank: List[Int],
    order: List[Int],
    left: Int,
    right: Int,
    mut marks: List[UInt8],
    mut touched: List[Int],
    mut stack: List[Int],
    mut keys: List[Int],
) -> List[Int]:
    var count = len(graph)
    if len(marks) < count:
        marks.resize(count, 0)
    var low = rank[left] if left >= 0 else -1
    var high = rank[right] if right >= 0 else count
    comptime FORWARD = UInt8(1)
    comptime BOTH = UInt8(2)
    stack.clear()
    if left >= 0:
        var e = Int(graph.first_out[left])
        while e >= 0:
            stack.append(Int(graph.edge_to[e]))
            e = Int(graph.next_out[e])
    else:
        for node in range(count):
            if graph.in_degree[node] == 0:
                stack.append(node)
    while len(stack) > 0:
        var node = stack.pop()
        if rank[node] >= high or marks[node] != 0:
            continue
        marks[node] = FORWARD
        touched.append(node)
        var e = Int(graph.first_out[node])
        while e >= 0:
            stack.append(Int(graph.edge_to[e]))
            e = Int(graph.next_out[e])
    if right >= 0:
        var e = Int(graph.first_in[right])
        while e >= 0:
            stack.append(Int(graph.edge_from[e]))
            e = Int(graph.next_in[e])
    else:
        for node in range(count):
            if graph.out_degree[node] == 0:
                stack.append(node)
    keys.clear()
    while len(stack) > 0:
        var node = stack.pop()
        if rank[node] <= low or marks[node] != FORWARD:
            continue
        marks[node] = BOTH
        keys.append(rank[node])
        var e = Int(graph.first_in[node])
        while e >= 0:
            stack.append(Int(graph.edge_from[e]))
            e = Int(graph.next_in[e])
    for node in touched:
        marks[node] = 0
    touched.clear()
    sort(keys)
    var out = List[Int](capacity=len(keys))
    for r in keys:
        out.append(order[r])
    return out^


def window_nodes(
    graph: Graph, rank: List[Int], order: List[Int], left: Int, right: Int, mut scratch: WindowScratch
) -> List[Int]:
    """The nodes strictly between `left` and `right` (or `SOURCE`, `SINK`) that a path from one to the other can
    visit, in topological order (`rank` each node's place in `order`), touching only the nodes between the two ends
    (and, at the graph's ends, its sources or sinks)."""
    return _window_nodes_into(
        graph, rank, order, left, right, scratch.marks, scratch.touched, scratch.stack, scratch.keys
    )


comptime UNIQUE_K = 19
"""A pin is kept where the k-mer centred on it is unique in both sequences (minimap2's seed length)."""
comptime _BLOOM_WORDS = 64
comptime _BLOOM_BITS = _BLOOM_WORDS * 64
comptime _BLOOM_MASK = UInt64(_BLOOM_BITS - 1)


def _kmer(seq: List[UInt8], start: Int) -> UInt64:
    """The 2-bit code of `seq[start : start + UNIQUE_K]`, past 2 * UNIQUE_K bits for a letter not ACGT (no k-mer of
    ACGT letters has those bits)."""
    var code = UInt64(0)
    for i in range(start, start + UNIQUE_K):
        var b = encode_base(seq[i])
        if b > 3:
            return UInt64(1) << UInt64(2 * UNIQUE_K)
        code = (code << 2) | UInt64(b)
    return code


@inline(.always)
def _bloom_bits(key: UInt64) -> Tuple[Int, Int]:
    var h = key ^ (key >> 19)
    return (Int(h & _BLOOM_MASK), Int((h >> 12) & _BLOOM_MASK))


def _count(seq: List[UInt8], mut wanted: Dict[UInt64, Int], bloom: Array[UInt64, _BLOOM_WORDS]):
    """`wanted`'s k-mers' occurrences in `seq` (each count reset first), pre-filtered by `bloom`."""
    for key in wanted.keys():
        wanted[key] = 0
    var mask = (UInt64(1) << UInt64(2 * UNIQUE_K)) - 1
    var code = UInt64(0)
    var valid = 0
    for i in range(len(seq)):
        var b = encode_base(seq[i])
        valid = valid + 1 if b < 4 else 0
        code = ((code << 2) | UInt64(b & 3)) & mask
        if valid >= UNIQUE_K:
            var bits = _bloom_bits(code)
            var b1 = bits[0]
            var b2 = bits[1]
            if (bloom[b1 >> 6] & (UInt64(1) << UInt64(b1 & 63))) == 0:
                continue
            if (bloom[b2 >> 6] & (UInt64(1) << UInt64(b2 & 63))) == 0:
                continue
            var seen = wanted.get(code)
            if seen:
                wanted[code] = seen.value() + 1


def unique_pins(reference: List[UInt8], query: List[UInt8], found: Pins) -> List[Bool]:
    """For each pin, whether the `UNIQUE_K`-mer centred on it occurs once in `reference` and once in `query`: in a
    tandem array of near-identical copies the chain may pair a seed with the wrong copy, and a pin there would force
    the alignment to; there the k-mers recur, so the array is left to the DP whole."""
    comptime half = UNIQUE_K // 2
    var keep = List[Bool](length=len(found.r), fill=False)
    var wanted = Dict[UInt64, Int]()
    var key_of = List[UInt64](length=len(found.r), fill=0)
    var bloom = Array[UInt64, _BLOOM_WORDS](fill=0)
    for m in range(len(found.r)):
        var r = found.r[m]
        var q = found.q[m]
        if r < half or q < half or r + half >= len(reference) or q + half >= len(query):
            continue
        var key = _kmer(reference, r - half)
        if key >> UInt64(2 * UNIQUE_K) != 0 or key != _kmer(query, q - half):
            continue
        key_of[m] = key
        keep[m] = True
        wanted[key] = 0
        var bits = _bloom_bits(key)
        bloom[bits[0] >> 6] |= UInt64(1) << UInt64(bits[0] & 63)
        bloom[bits[1] >> 6] |= UInt64(1) << UInt64(bits[1] & 63)
    if len(wanted) == 0:
        return keep^
    # Counted in each sequence: occurrences of the pins' k-mers alone.
    for side in range(2):
        _count(reference if side == 0 else query, wanted, bloom)
        for m in range(len(found.r)):
            if keep[m] and wanted.get(key_of[m], 0) != 1:
                keep[m] = False
    return keep^


def compute_pin_anchors(letters: List[UInt8], index: ReferenceIndex, query: List[UInt8]) -> PinAnchors:
    """Computes the strictly increasing one-letter unique pin anchors of `query` against pre-indexed `letters`."""
    if len(letters) >= 1 << POSITION_BITS or len(query) >= 1 << POSITION_BITS:
        return PinAnchors(List[Int](), List[Int](), 0)
    var found = pins_with_index(index, query)
    var unique = unique_pins(letters, query, found)
    var anchor_r = List[Int]()
    var anchor_q = List[Int]()
    for m in range(len(found.r)):
        var r = found.r[m]
        var q = found.q[m]
        if not unique[m] or letters[r] != query[q]:
            continue
        if len(anchor_r) > 0 and (r <= anchor_r[len(anchor_r) - 1] or q <= anchor_q[len(anchor_q) - 1]):
            continue
        anchor_r.append(r)
        anchor_q.append(q)
    return PinAnchors(anchor_r^, anchor_q^, len(found.r))


def prefetch_pin_anchors(
    mut scratch: WindowScratch,
    seqs: List[List[UInt8]],
    read_ids: List[Int],
    start_idx: Int,
    threads: Int,
):
    """Speculatively computes `PinAnchors` for `read_ids[start_idx:]` in parallel against `scratch.cached_letters`
    and `scratch.cached_index`, storing them in `scratch.precomp_anchors` indexed by read id."""
    var rem = len(read_ids) - start_idx
    if not scratch.has_cached_index or rem < 2 or len(scratch.cached_letters) >= 1 << POSITION_BITS:
        return
    var workers = thread_count(threads, rem)
    if workers <= 1:
        return
    scratch.precomp_letters = scratch.cached_letters.copy()
    scratch.precomp_anchors = List[PinAnchors](length=len(seqs), fill=PinAnchors())
    var out_p = scratch.precomp_anchors.unsafe_ptr()
    ref letters = scratch.precomp_letters
    ref index = scratch.cached_index
    var taken = Atomic[Int64](0)

    def seed_worker(
        slot: Int,
    ) {mut taken, imm start_idx, imm rem, imm read_ids, imm seqs, imm letters, imm index, imm out_p}:
        while True:
            var offset = Int(taken.fetch_add(1))
            if offset >= rem:
                return
            var r = read_ids[start_idx + offset]
            if len(seqs[r]) < 1 << POSITION_BITS:
                out_p[unsafe_offset=r] = compute_pin_anchors(letters, index, seqs[r])

    spread(seed_worker, workers, workers)


def _align_windows_serial(
    graph: Graph,
    topology: Topology,
    hp_path: List[Int],
    anchor_r: List[Int],
    anchor_q: List[Int],
    query: List[UInt8],
    scoring: Scoring,
    window_b: Int,
    mode: Int,
    mut row_of: List[Int],
    mut marks: List[UInt8],
    mut touched: List[Int],
    mut stack: List[Int],
    mut keys: List[Int],
    mut workspace: Workspace,
) raises -> GraphAlignment:
    """Runs window DP serially between precomputed `anchor_r`/`anchor_q` using caller-owned scratch buffers."""
    var path_len = len(hp_path)
    var aligned = List[Int](length=len(query), fill=-1)
    var score = 0
    var cells = 0
    var peak = 0
    var num_windows = len(anchor_r) + 1
    var r0 = -1
    var q0 = -1
    for m in range(num_windows):
        var r1 = anchor_r[m] if m < len(anchor_r) else path_len
        var q1 = anchor_q[m] if m < len(anchor_r) else len(query)
        if q1 - q0 > 1:
            var left = hp_path[r0] if r0 >= 0 else SOURCE
            var right = hp_path[r1] if r1 < path_len else SINK
            var nodes = _window_nodes_into(
                graph, topology.rank, topology.order, left, right, marks, touched, stack, keys
            )
            if len(nodes) == 0:
                if not ((mode == SUFFIX and r0 < 0) or (mode == PREFIX and r1 >= path_len)):
                    score -= Int(scoring.gap(q1 - q0 - 1))
            else:
                var window = Span(query)[q0 + 1 : q1]
                var window_mode = GLOBAL
                if mode == SUFFIX and r0 < 0:
                    window_mode = SUFFIX
                elif mode == PREFIX and r1 >= path_len:
                    window_mode = PREFIX
                var got = workspace.align(
                    graph, topology, nodes, left, right, window, scoring, row_of, window_mode, window_b, 0.0
                )
                for j in range(len(got.nodes)):
                    aligned[q0 + 1 + j] = got.nodes[j]
                score += got.score
                cells += got.cells
                peak = max(peak, got.peak_bytes)
        elif q1 == q0 + 1 and r1 > r0 + 1:
            var free = (r0 < 0 and mode == SUFFIX) or (r1 >= path_len and mode == PREFIX)
            var joined = r0 >= 0 and r1 < path_len and graph.edge(hp_path[r0], hp_path[r1]) >= 0
            if not free and not joined:
                score -= Int(scoring.gap(r1 - r0 - 1))
        if m < len(anchor_r):
            aligned[q1] = hp_path[r1]
            score += Int(scoring.match_score)
        r0 = r1
        q0 = q1
    return GraphAlignment(aligned^, score, cells, peak)


def align_pinned_snapshot(
    graph: Graph,
    topology: Topology,
    hp_path: List[Int],
    hp_letters: List[UInt8],
    index: ReferenceIndex,
    query: List[UInt8],
    scoring: Scoring,
    band_b: Int,
    band_f: Float64,
    mut row_of: List[Int],
    mut marks: List[UInt8],
    mut touched: List[Int],
    mut stack: List[Int],
    mut keys: List[Int],
    mut workspace: Workspace,
    mode: Int = GLOBAL,
) raises -> GraphAlignment:
    """Aligns `query` against a read-only `(graph, topology, hp_path, hp_letters, index)` snapshot using strictly
    worker-owned buffers (thread-safe for parallel mini-batch POA)."""
    var anchors = compute_pin_anchors(hp_letters, index, query)
    var window_b = band_b + Int(band_f * Float64(len(query))) if band_b >= 0 else band_b
    return _align_windows_serial(
        graph,
        topology,
        hp_path,
        anchors.anchor_r,
        anchors.anchor_q,
        query,
        scoring,
        window_b,
        mode,
        row_of,
        marks,
        touched,
        stack,
        keys,
        workspace,
    )


def align_pinned(
    graph: Graph,
    topology: Topology,
    query: List[UInt8],
    scoring: Scoring,
    band_b: Int,
    band_f: Float64,
    mut row_of: List[Int],
    mut scratch: WindowScratch,
    mut workspace: Workspace,
    mode: Int = GLOBAL,
    threads: Int = 1,
    read_id: Int = -1,
) raises -> GraphAlignment:
    """`query` aligned to the graph through the pins of its alignment to the consensus path (the heaviest bundle): a
    pin whose letters agree anchors that query letter to its path node, and between consecutive anchors the kernel
    aligns the query's letters to the nodes between their nodes. `mode` (`GLOBAL`, `PREFIX` or `SUFFIX`) is the
    alignment's as a whole: the first window starts as it does, the last ends as it does, those between are global.
    `row_of` holds a slot per node."""
    graph.heaviest_path_into(
        topology.order,
        scratch.hp_starts,
        scratch.hp_ends,
        scratch.hp_score,
        scratch.hp_next,
        scratch.hp_path,
    )
    scratch.hp_letters.clear()
    scratch.hp_letters.reserve(len(scratch.hp_path))
    for node in scratch.hp_path:
        scratch.hp_letters.append(graph.bases[node])

    var anchor_r: List[Int]
    var anchor_q: List[Int]
    var raw_count: Int
    var used_precomp = False
    if (
        read_id >= 0
        and read_id < len(scratch.precomp_anchors)
        and scratch.precomp_anchors[read_id].valid
        and scratch.precomp_letters == scratch.hp_letters
    ):
        anchor_r = scratch.precomp_anchors[read_id].anchor_r.copy()
        anchor_q = scratch.precomp_anchors[read_id].anchor_q.copy()
        raw_count = scratch.precomp_anchors[read_id].raw_count
        used_precomp = True
    else:
        if len(scratch.hp_letters) >= 1 << POSITION_BITS or len(query) >= 1 << POSITION_BITS:
            var trace_limit = Trace()
            trace_limit.line(
                PINS,
                String("no pins: ", len(scratch.hp_letters), " x ", len(query), " letters, past 2^", POSITION_BITS),
            )
            anchor_r = List[Int]()
            anchor_q = List[Int]()
            raw_count = 0
        else:
            if not scratch.has_cached_index or scratch.cached_letters != scratch.hp_letters:
                scratch.cached_letters = scratch.hp_letters.copy()
                scratch.cached_index = ReferenceIndex(scratch.cached_letters)
                scratch.has_cached_index = True
            var anchors = compute_pin_anchors(scratch.hp_letters, scratch.cached_index, query)
            raw_count = anchors.raw_count
            anchor_r = List[Int]()
            anchor_q = List[Int]()
            swap(anchor_r, anchors.anchor_r)
            swap(anchor_q, anchors.anchor_q)

    var trace = Trace()
    if trace.on(PINS):
        if used_precomp:
            pass
        trace.line(
            PINS,
            String(
                "pinned alignment: ",
                len(query),
                " letters to a path of ",
                len(scratch.hp_path),
                "; ",
                raw_count,
                " pins, ",
                len(anchor_r),
                " kept (19-mer unique in both, letters agreeing), ",
                len(anchor_r) + 1,
                " windows",
            ),
        )
    # Every window gets the band the whole sequence would (b + f * its length, not the window's): a window of a few
    # hundred letters would otherwise get a dozen either side, too few for a copy-number change inside it.
    var window_b = band_b + Int(band_f * Float64(len(query))) if band_b >= 0 else band_b
    var aligned = List[Int](length=len(query), fill=-1)
    var score = 0
    var cells = 0
    var peak = 0
    var num_windows = len(anchor_r) + 1
    var dp_windows = List[Int]()
    var dp_bases = 0
    var path_len = len(scratch.hp_path)
    var r0 = -1
    var q0 = -1
    for m in range(num_windows):
        var r1 = anchor_r[m] if m < len(anchor_r) else path_len
        var q1 = anchor_q[m] if m < len(anchor_r) else len(query)
        if q1 - q0 > 1:
            dp_windows.append(m)
            dp_bases += q1 - q0 - 1
        elif q1 == q0 + 1 and r1 > r0 + 1:
            var free = (r0 < 0 and mode == SUFFIX) or (r1 >= path_len and mode == PREFIX)
            var joined = r0 >= 0 and r1 < path_len and graph.edge(scratch.hp_path[r0], scratch.hp_path[r1]) >= 0
            if not free and not joined:
                score -= Int(scoring.gap(r1 - r0 - 1))
        if m < len(anchor_r):
            aligned[q1] = scratch.hp_path[r1]
            score += Int(scoring.match_score)
        r0 = r1
        q0 = q1

    var workers = thread_count(threads, len(dp_windows))
    if workers > 1 and len(dp_windows) >= 2 and dp_bases >= 512:
        # Align independent pin windows concurrently, longest first; each writes to its own disjoint query slice
        # and reuses its own pooled worker workspace across reads and calls.
        while len(scratch.workers) < workers:
            scratch.workers.append(_PinnedWorker())
        var order = List[UInt64](capacity=len(dp_windows))
        for m in dp_windows:
            var q_prev = anchor_q[m - 1] if m > 0 else -1
            var q_next = anchor_q[m] if m < len(anchor_r) else len(query)
            order.append((UInt64(q_next - q_prev) << 32) | UInt64(m))
        sort(order)
        order.reverse()
        var w_scores = List[Int](length=num_windows, fill=0)
        var w_cells = List[Int](length=num_windows, fill=0)
        var w_peaks = List[Int](length=num_windows, fill=0)
        var w_failed = List[Bool](length=num_windows, fill=False)
        var path_p = scratch.hp_path.unsafe_ptr()
        var pool_p = scratch.workers.unsafe_ptr()
        var aligned_p = aligned.unsafe_ptr()
        var scores_p = w_scores.unsafe_ptr()
        var cells_p = w_cells.unsafe_ptr()
        var peaks_p = w_peaks.unsafe_ptr()
        var failed_p = w_failed.unsafe_ptr()
        var taken = Atomic[Int64](0)
        var graph_len = len(graph)

        def run_pin_windows(
            slot: Int,
        ) {
            mut taken,
            imm order,
            imm anchor_r,
            imm anchor_q,
            imm path_p,
            imm path_len,
            imm graph,
            imm topology,
            imm query,
            imm scoring,
            imm window_b,
            imm mode,
            imm graph_len,
            imm pool_p,
            imm aligned_p,
            imm scores_p,
            imm cells_p,
            imm peaks_p,
            imm failed_p,
        }:
            ref w = pool_p[unsafe_offset=slot]
            w.row_of.resize(graph_len, -1)
            while True:
                var idx = Int(taken.fetch_add(1))
                if idx >= len(order):
                    return
                var m = Int(order[idx] & 0xFFFFFFFF)
                var wr0 = anchor_r[m - 1] if m > 0 else -1
                var wq0 = anchor_q[m - 1] if m > 0 else -1
                var wr1 = anchor_r[m] if m < len(anchor_r) else path_len
                var wq1 = anchor_q[m] if m < len(anchor_r) else len(query)
                var left = path_p[unsafe_offset=wr0] if wr0 >= 0 else SOURCE
                var right = path_p[unsafe_offset=wr1] if wr1 < path_len else SINK
                var nodes = _window_nodes_into(
                    graph, topology.rank, topology.order, left, right, w.marks, w.touched, w.stack, w.keys
                )
                if len(nodes) == 0:
                    if not ((mode == SUFFIX and wr0 < 0) or (mode == PREFIX and wr1 >= path_len)):
                        scores_p[unsafe_offset=m] = -Int(scoring.gap(wq1 - wq0 - 1))
                    continue
                var window = Span(query)[wq0 + 1 : wq1]
                var window_mode = GLOBAL
                if mode == SUFFIX and wr0 < 0:
                    window_mode = SUFFIX
                elif mode == PREFIX and wr1 >= path_len:
                    window_mode = PREFIX
                try:
                    var got = w.workspace.align(
                        graph, topology, nodes, left, right, window, scoring, w.row_of, window_mode, window_b, 0.0
                    )
                    for j in range(len(got.nodes)):
                        aligned_p[unsafe_offset=wq0 + 1 + j] = got.nodes[j]
                    scores_p[unsafe_offset=m] = got.score
                    cells_p[unsafe_offset=m] = got.cells
                    peaks_p[unsafe_offset=m] = got.peak_bytes
                except:
                    for i in range(len(w.row_of)):
                        w.row_of[i] = -1
                    failed_p[unsafe_offset=m] = True

        spread(run_pin_windows, workers, workers)
        for m in dp_windows:
            if w_failed[m]:
                var wr0 = anchor_r[m - 1] if m > 0 else -1
                var wq0 = anchor_q[m - 1] if m > 0 else -1
                var wr1 = anchor_r[m] if m < len(anchor_r) else path_len
                var wq1 = anchor_q[m] if m < len(anchor_r) else len(query)
                var left = scratch.hp_path[wr0] if wr0 >= 0 else SOURCE
                var right = scratch.hp_path[wr1] if wr1 < path_len else SINK
                var nodes = window_nodes(graph, topology.rank, topology.order, left, right, scratch)
                var window = Span(query)[wq0 + 1 : wq1]
                var window_mode = GLOBAL
                if mode == SUFFIX and wr0 < 0:
                    window_mode = SUFFIX
                elif mode == PREFIX and wr1 >= path_len:
                    window_mode = PREFIX
                var got = workspace.align(
                    graph, topology, nodes, left, right, window, scoring, row_of, window_mode, window_b, 0.0
                )
                for j in range(len(got.nodes)):
                    aligned[wq0 + 1 + j] = got.nodes[j]
                score += got.score
                cells += got.cells
                peak = max(peak, got.peak_bytes)
            else:
                score += w_scores[m]
                cells += w_cells[m]
                peak = max(peak, w_peaks[m])
        return GraphAlignment(aligned^, score, cells, peak)

    for m in dp_windows:
        var wr0 = anchor_r[m - 1] if m > 0 else -1
        var wq0 = anchor_q[m - 1] if m > 0 else -1
        var wr1 = anchor_r[m] if m < len(anchor_r) else path_len
        var wq1 = anchor_q[m] if m < len(anchor_r) else len(query)
        var left = scratch.hp_path[wr0] if wr0 >= 0 else SOURCE
        var right = scratch.hp_path[wr1] if wr1 < path_len else SINK
        var nodes = window_nodes(graph, topology.rank, topology.order, left, right, scratch)
        if len(nodes) == 0:
            if not ((mode == SUFFIX and wr0 < 0) or (mode == PREFIX and wr1 >= path_len)):
                score -= Int(scoring.gap(wq1 - wq0 - 1))
        else:
            var window = Span(query)[wq0 + 1 : wq1]
            var window_mode = GLOBAL
            if mode == SUFFIX and wr0 < 0:
                window_mode = SUFFIX
            elif mode == PREFIX and wr1 >= path_len:
                window_mode = PREFIX
            var got = workspace.align(
                graph, topology, nodes, left, right, window, scoring, row_of, window_mode, window_b, 0.0
            )
            for j in range(len(got.nodes)):
                aligned[wq0 + 1 + j] = got.nodes[j]
            score += got.score
            cells += got.cells
            peak = max(peak, got.peak_bytes)
    return GraphAlignment(aligned^, score, cells, peak)
