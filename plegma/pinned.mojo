"""A sequence aligned to the graph through pins: minimap2's pins (`pins.pins`) between the sequence and the graph's
consensus path, each pin's letters matched to its path node, and the kernel's DP run only in the windows between two
pinned nodes (the subgraph a path from one to the other can visit). The DP's cells then follow the windows, not the
whole graph times the sequence; where the chain finds no pins, the one window is the whole graph."""

from .alphabet import encode_base
from .graph import Graph
from .kernel import GLOBAL, PREFIX, SINK, SOURCE, SUFFIX, GraphAlignment, Scoring, Topology, Workspace
from .pins import Pins, pins


struct WindowScratch(Movable):
    """Reusable marks for `window_nodes`: one byte a node, reset after each call for the nodes it touched."""

    var marks: List[UInt8]
    var touched: List[Int]

    def __init__(out self):
        self.marks = List[UInt8]()
        self.touched = List[Int]()


def window_nodes(
    graph: Graph, rank: List[Int], order: List[Int], left: Int, right: Int, mut scratch: WindowScratch
) -> List[Int]:
    """The nodes strictly between `left` and `right` (or `SOURCE`, `SINK`) that a path from one to the other can
    visit, in topological order (`rank` each node's place in `order`), touching only the nodes between the two ends
    (and, at the graph's ends, its sources or sinks)."""
    var count = len(graph)
    if len(scratch.marks) < count:
        scratch.marks.resize(count, 0)
    var low = rank[left] if left >= 0 else -1
    var high = rank[right] if right >= 0 else count
    comptime FORWARD = UInt8(1)
    comptime BOTH = UInt8(2)
    var stack = List[Int]()
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
        if rank[node] >= high or scratch.marks[node] != 0:
            continue
        scratch.marks[node] = FORWARD
        scratch.touched.append(node)
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
    var keys = List[Int]()
    while len(stack) > 0:
        var node = stack.pop()
        if rank[node] <= low or scratch.marks[node] != FORWARD:
            continue
        scratch.marks[node] = BOTH
        keys.append(rank[node])
        var e = Int(graph.first_in[node])
        while e >= 0:
            stack.append(Int(graph.edge_from[e]))
            e = Int(graph.next_in[e])
    for node in scratch.touched:
        scratch.marks[node] = 0
    scratch.touched.clear()
    sort(keys)
    var out = List[Int](capacity=len(keys))
    for r in keys:
        out.append(order[r])
    return out^


comptime UNIQUE_K = 19
"""A pin is kept where the k-mer centred on it is unique in both sequences (minimap2's seed length)."""


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


def _count(seq: List[UInt8], mut wanted: Dict[UInt64, Int]):
    """`wanted`'s k-mers' occurrences in `seq` (each count reset first)."""
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
    if len(wanted) == 0:
        return keep^
    # Counted in each sequence: occurrences of the pins' k-mers alone.
    for side in range(2):
        _count(reference if side == 0 else query, wanted)
        for m in range(len(found.r)):
            if keep[m] and wanted.get(key_of[m], 0) != 1:
                keep[m] = False
    return keep^


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
) raises -> GraphAlignment:
    """`query` aligned to the graph through the pins of its alignment to the consensus path (the heaviest bundle): a
    pin whose letters agree anchors that query letter to its path node, and between consecutive anchors the kernel
    aligns the query's letters to the nodes between their nodes. `mode` (`GLOBAL`, `PREFIX` or `SUFFIX`) is the
    alignment's as a whole: the first window starts as it does, the last ends as it does, those between are global.
    `row_of` holds a slot per node."""
    var path = graph.heaviest_path(topology.order)
    var letters = List[UInt8](capacity=len(path))
    for node in path:
        letters.append(graph.bases[node])
    var found = pins(letters, query)
    var unique = unique_pins(letters, query, found)
    # One-letter anchors, each past the last on both: the windows then meet at path nodes, so the alignment's path
    # runs forward through the graph (two windows meeting between nodes could join nodes into a cycle).
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
    # Every window gets the band the whole sequence would (b + f * its length, not the window's): a window of a few
    # hundred letters would otherwise get a dozen either side, too few for a copy-number change inside it.
    var window_b = band_b + Int(band_f * Float64(len(query))) if band_b >= 0 else band_b
    var aligned = List[Int](length=len(query), fill=-1)
    var score = 0
    var cells = 0
    var peak = 0
    var r0 = -1
    var q0 = -1
    for m in range(len(anchor_r) + 1):
        var r1 = anchor_r[m] if m < len(anchor_r) else len(path)
        var q1 = anchor_q[m] if m < len(anchor_r) else len(query)
        if q1 - q0 > 1:
            var left = path[r0] if r0 >= 0 else SOURCE
            var right = path[r1] if r1 < len(path) else SINK
            var nodes = window_nodes(graph, topology.rank, topology.order, left, right, scratch)
            # No nodes between the anchors: the query's letters there are inserted.
            if len(nodes) > 0:
                var window = Span(query)[q0 + 1 : q1]
                var window_mode = GLOBAL
                if mode == SUFFIX and r0 < 0:
                    window_mode = SUFFIX
                elif mode == PREFIX and r1 >= len(path):
                    window_mode = PREFIX
                var got = workspace.align(
                    graph, topology, nodes, left, right, window, scoring, row_of, window_mode, window_b, 0.0
                )
                for j in range(len(got.nodes)):
                    aligned[q0 + 1 + j] = got.nodes[j]
                score += got.score
                cells += got.cells
                peak = max(peak, got.peak_bytes)
        if m < len(anchor_r):
            aligned[q1] = path[r1]
            score += Int(scoring.match_score)
        r0 = r1
        q0 = q1
    return GraphAlignment(aligned^, score, cells, peak)
