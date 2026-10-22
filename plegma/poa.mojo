"""Aligning a stretch of sequence to the graph between two nodes.

Gap-affine dynamic programming over the subgraph between the nodes: one row per node in topological order, one
column per letter, a row taking the best of its predecessors' rows (abPOA's and SPOA's layout; scalar for now).
The ends are fixed: the alignment leaves the left node and reaches the right one, so a window between two
anchors is aligned globally, never with a free end.
"""

from std.collections import Dict

from .graph import Graph


comptime SOURCE = -1
"""The left end before every path's first node."""
comptime SINK = -2
"""The right end after every path's last node."""
comptime _INF = Int32(1 << 29)


struct GapAffine(ImplicitlyCopyable):
    """A substitution costs `mismatch`; a gap of `k` letters `open + extend k`, as karan's `Costs.affine`."""

    var mismatch: Int32
    var open: Int32
    var extend: Int32

    def __init__(out self, mismatch: Int, open: Int, extend: Int):
        self.mismatch = Int32(mismatch)
        self.open = Int32(open)
        self.extend = Int32(extend)


struct WindowAlignment(Copyable, Movable):
    """For each letter of a window, the node it lands on (-1: inserted), and the alignment's cost."""

    var nodes: List[Int]
    var cost: Int

    def __init__(out self, var nodes: List[Int], cost: Int):
        self.nodes = nodes^
        self.cost = cost


def window_nodes(graph: Graph, rank: List[Int], order: List[Int], left: Int, right: Int) -> List[Int]:
    """The nodes strictly between `left` and `right` (or `SOURCE`, `SINK`) that a path from one to the other
    can visit, in topological order (`rank` each node's place in `order`)."""
    var count = len(graph)
    var low = rank[left] if left >= 0 else -1
    var high = rank[right] if right >= 0 else count
    var forward = List[Bool](length=count, fill=False)
    var stack = List[Int]()
    if left >= 0:
        for s in graph.successors[left]:
            stack.append(s)
    else:
        for node in range(count):
            if len(graph.predecessors[node]) == 0:
                stack.append(node)
    while len(stack) > 0:
        var node = stack.pop()
        if rank[node] >= high or forward[node]:
            continue
        forward[node] = True
        for s in graph.successors[node]:
            stack.append(s)
    var both = List[Bool](length=count, fill=False)
    if right >= 0:
        for p in graph.predecessors[right]:
            stack.append(p)
    else:
        for node in range(count):
            if len(graph.successors[node]) == 0:
                stack.append(node)
    while len(stack) > 0:
        var node = stack.pop()
        if rank[node] <= low or both[node] or not forward[node]:
            continue
        both[node] = True
        for p in graph.predecessors[node]:
            stack.append(p)
    var out = List[Int]()
    for r in range(low + 1, high):
        if both[order[r]]:
            out.append(order[r])
    return out^


def align_window(
    graph: Graph,
    rank: List[Int],
    order: List[Int],
    left: Int,
    right: Int,
    query: Span[UInt8, _],
    costs: GapAffine,
    max_cells: Int,
) raises -> Optional[WindowAlignment]:
    """An optimal alignment of `query` from `left` to `right`: the node each letter lands on (-1: inserted);
    None when the subgraph times the letters passes `max_cells`. Of equal alignments, the one whose gaps come
    first."""
    var nodes = window_nodes(graph, rank, order, left, right)
    var n = len(nodes)
    var m = len(query)
    var width = m + 1
    if (n + 1) * width > max_cells:
        return None
    var row_of = Dict[Int, Int]()
    for r in range(n):
        row_of[nodes[r]] = r + 1
    # Predecessor rows of each row: row 0 is the left end.
    var preds = List[List[Int]](length=n + 1, fill=List[Int]())
    for r in range(n):
        var node = nodes[r]
        if len(graph.predecessors[node]) == 0 and left == SOURCE:
            preds[r + 1].append(0)
        for p in graph.predecessors[node]:
            if p == left:
                preds[r + 1].append(0)
            elif p in row_of:
                preds[r + 1].append(row_of[p])
    var ends = List[Int]()
    for r in range(n):
        var node = nodes[r]
        if right == SINK:
            if len(graph.successors[node]) == 0:
                ends.append(r + 1)
        else:
            for s in graph.successors[node]:
                if s == right:
                    ends.append(r + 1)
    var direct = False
    if left >= 0 and right >= 0:
        for s in graph.successors[left]:
            direct = direct or s == right
    elif left >= 0:
        direct = len(graph.successors[left]) == 0
    elif right >= 0:
        direct = len(graph.predecessors[right]) == 0
    else:
        direct = n == 0
    if direct:
        ends.append(0)
    if len(ends) == 0:
        raise Error(t"no path from node {left} to node {right}")

    var cells = (n + 1) * width
    var h = List[Int32](length=cells, fill=_INF)
    var e = List[Int32](length=cells, fill=_INF)
    var f = List[Int32](length=cells, fill=_INF)
    var gap = costs.open + costs.extend
    h[0] = 0
    for j in range(1, width):
        f[j] = costs.open + costs.extend * Int32(j)
        h[j] = f[j]
    for r in range(1, n + 1):
        var base = graph.bases[nodes[r - 1]]
        var row = r * width
        for p in preds[r]:
            var prow = p * width
            for j in range(width):
                e[row + j] = min(e[row + j], h[prow + j] + gap, e[prow + j] + costs.extend)
            for j in range(1, width):
                var sub = Int32(0) if query[j - 1] == base else costs.mismatch
                h[row + j] = min(h[row + j], h[prow + j - 1] + sub)
        h[row] = min(h[row], e[row])
        for j in range(1, width):
            f[row + j] = min(h[row + j - 1] + gap, f[row + j - 1] + costs.extend)
            h[row + j] = min(h[row + j], e[row + j], f[row + j])

    var best = _INF
    var at = -1
    for r in ends:
        if h[r * width + m] < best:
            best = h[r * width + m]
            at = r
    var aligned = List[Int](length=m, fill=-1)
    var j = m
    var state = 0  # 0: H, 1: E (a node deleted), 2: F (a letter inserted)
    while at > 0 or j > 0:
        var cell = at * width + j
        if state == 0:
            if at > 0 and j > 0:
                var sub = Int32(0) if query[j - 1] == graph.bases[nodes[at - 1]] else costs.mismatch
                var moved = False
                for p in preds[at]:
                    if h[p * width + j - 1] + sub == h[cell]:
                        aligned[j - 1] = nodes[at - 1]
                        at = p
                        j -= 1
                        moved = True
                        break
                if moved:
                    continue
            if at > 0 and h[cell] == e[cell]:
                state = 1
            elif j > 0 and h[cell] == f[cell]:
                state = 2
            else:
                raise Error(t"traceback lost at row {at}, letter {j}")
        elif state == 1:
            var moved = False
            for p in preds[at]:
                if h[p * width + j] + gap == e[cell]:
                    at = p
                    state = 0
                    moved = True
                    break
                if e[p * width + j] + costs.extend == e[cell]:
                    at = p
                    moved = True
                    break
            if not moved:
                raise Error(t"traceback lost in a deletion at row {at}, letter {j}")
        else:
            state = 0 if h[cell - 1] + gap == f[cell] else 2
            j -= 1
    return WindowAlignment(aligned^, Int(best))
