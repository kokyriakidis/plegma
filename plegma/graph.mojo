"""The partial order graph: one node per letter, nodes aligned to each other grouped into columns, and every
sequence's path through it kept.

A sequence joins by its letters' nodes (`fuse_nodes`, the nodes a graph alignment pairs them with): a letter lands
on its node, or on the node of that letter in the node's column (made if missing); an inserted letter on a node of its
own column. An MSA's rows make a graph directly (`from_msa`).
"""


comptime MATCH = UInt8(ord("="))
comptime MISMATCH = UInt8(ord("X"))
comptime INSERTION = UInt8(ord("I"))
comptime DELETION = UInt8(ord("D"))
comptime GAP = UInt8(ord("-"))


struct Graph(Movable, Sized):
    """Nodes (`bases`, `columns`), weighted edges and the paths of the sequences added so far, in flat arrays:
    a node costs a few 32-bit fields, an edge five.

    Walk a node's out-edges with `var e = graph.first_out[node]`, `while e >= 0`, `e = graph.next_out[e]`
    (in insertion order; `edge_to[e]`, `edge_weight[e]`), its in-edges likewise through `first_in`, `next_in`
    and `edge_from`; a column's nodes through `column_first` and `column_next`."""

    var bases: List[UInt8]
    """Each node's letter, as the sequence spelled it."""
    var columns: List[Int32]
    """Each node's column."""
    var column_first: List[Int32]
    var column_next: List[Int32]
    """The next node in the node's column (-1: the last)."""
    var first_out: List[Int32]
    var last_out: List[Int32]
    var first_in: List[Int32]
    var last_in: List[Int32]
    var out_degree: List[Int32]
    var in_degree: List[Int32]
    var edge_from: List[Int32]
    var edge_to: List[Int32]
    var edge_weight: List[Int32]
    """How many paths take each edge."""
    var next_out: List[Int32]
    var next_in: List[Int32]
    var paths: List[List[Int32]]
    var names: List[String]

    def __init__(out self):
        self.bases = List[UInt8]()
        self.columns = List[Int32]()
        self.column_first = List[Int32]()
        self.column_next = List[Int32]()
        self.first_out = List[Int32]()
        self.last_out = List[Int32]()
        self.first_in = List[Int32]()
        self.last_in = List[Int32]()
        self.out_degree = List[Int32]()
        self.in_degree = List[Int32]()
        self.edge_from = List[Int32]()
        self.edge_to = List[Int32]()
        self.edge_weight = List[Int32]()
        self.next_out = List[Int32]()
        self.next_in = List[Int32]()
        self.paths = List[List[Int32]]()
        self.names = List[String]()

    def __len__(self) -> Int:
        """The number of nodes."""
        return len(self.bases)

    def column_count(self) -> Int:
        return len(self.column_first)

    def successors(self, node: Int) -> List[Int]:
        """`node`'s successors, in the order their edges were made (a copy: walk `first_out` in hot loops)."""
        var out = List[Int](capacity=Int(self.out_degree[node]))
        var e = Int(self.first_out[node])
        while e >= 0:
            out.append(Int(self.edge_to[e]))
            e = Int(self.next_out[e])
        return out^

    def predecessors(self, node: Int) -> List[Int]:
        """`node`'s predecessors, in the order their edges were made (a copy)."""
        var out = List[Int](capacity=Int(self.in_degree[node]))
        var e = Int(self.first_in[node])
        while e >= 0:
            out.append(Int(self.edge_from[e]))
            e = Int(self.next_in[e])
        return out^

    def _node(mut self, base: UInt8, column: Int) -> Int:
        var node = len(self.bases)
        self.bases.append(base)
        var at = column
        if at < 0:
            at = len(self.column_first)
            self.column_first.append(-1)
        self.columns.append(Int32(at))
        # Appended at the column's head: a column holds at most one node per letter.
        self.column_next.append(self.column_first[at])
        self.column_first[at] = Int32(node)
        self.first_out.append(-1)
        self.last_out.append(-1)
        self.first_in.append(-1)
        self.last_in.append(-1)
        self.out_degree.append(0)
        self.in_degree.append(0)
        return node

    def edge(self, before: Int, after: Int) -> Int:
        """The edge from `before` to `after`, or -1."""
        var e = Int(self.first_out[before])
        while e >= 0:
            if Int(self.edge_to[e]) == after:
                return e
            e = Int(self.next_out[e])
        return -1

    def _link(mut self, before: Int, after: Int):
        var found = self.edge(before, after)
        if found >= 0:
            self.edge_weight[found] += 1
            return
        var e = Int32(len(self.edge_to))
        self.edge_from.append(Int32(before))
        self.edge_to.append(Int32(after))
        self.edge_weight.append(1)
        self.next_out.append(-1)
        self.next_in.append(-1)
        if self.last_out[before] >= 0:
            self.next_out[Int(self.last_out[before])] = e
        else:
            self.first_out[before] = e
        self.last_out[before] = e
        if self.last_in[after] >= 0:
            self.next_in[Int(self.last_in[after])] = e
        else:
            self.first_in[after] = e
        self.last_in[after] = e
        self.out_degree[before] += 1
        self.in_degree[after] += 1

    def _on_path(mut self, mut path: List[Int32], node: Int):
        if len(path) > 0:
            self._link(Int(path[len(path) - 1]), node)
        path.append(Int32(node))

    def add(mut self, name: String, sequence: Span[UInt8, _]) -> Int:
        """Adds `sequence` as a path of new nodes, each its own column; returns the path's index."""
        var path = List[Int32](capacity=len(sequence))
        for c in sequence:
            self._on_path(path, self._node(c, -1))
        self.paths.append(path^)
        self.names.append(name)
        return len(self.paths) - 1

    def _lands(mut self, node: Int, base: UInt8) -> Int:
        """The node a letter aligned to `node` lands on: `node` itself, or the node of that letter in its
        column, made if missing."""
        if self.bases[node] == base:
            return node
        var column = Int(self.columns[node])
        var m = Int(self.column_first[column])
        while m >= 0:
            if self.bases[m] == base:
                return m
            m = Int(self.column_next[m])
        return self._node(base, column)

    def fuse_nodes(mut self, name: String, sequence: Span[UInt8, _], aligned: List[Int]) -> Int:
        """Adds `sequence`, its letter `j` aligned to node `aligned[j]` (-1: inserted); returns the new path's
        index."""
        var path = List[Int32](capacity=len(sequence))
        for j in range(len(sequence)):
            var node = aligned[j]
            self._on_path(path, self._lands(node, sequence[j]) if node >= 0 else self._node(sequence[j], -1))
        self.paths.append(path^)
        self.names.append(name)
        return len(self.paths) - 1

    def node_order(self) raises -> List[Int]:
        """Each node's rank in a topological order of the nodes (Kahn's, lowest node first among the ready)."""
        var count = len(self.bases)
        var incoming = List[Int32](capacity=count)
        var ready = List[Int]()
        for node in range(count):
            incoming.append(self.in_degree[node])
            if self.in_degree[node] == 0:
                ready.append(node)
        var rank = List[Int](length=count, fill=-1)
        var next = 0
        var head = 0
        while head < len(ready):
            var node = ready[head]
            head += 1
            rank[node] = next
            next += 1
            var e = Int(self.first_out[node])
            while e >= 0:
                var s = Int(self.edge_to[e])
                incoming[s] -= 1
                if incoming[s] == 0:
                    ready.append(s)
                e = Int(self.next_out[e])
        if next != count:
            raise Error(t"the nodes form a cycle: {count - next} of {count} unplaced")
        return rank^

    @staticmethod
    def from_msa(names: List[String], rows: List[String], skip: Int = -1) raises -> Graph:
        """The graph of an MSA's rows (all but row `skip`): in each column, one node per letter there, the rows'
        paths through them. Raises on rows of different lengths, or fewer names than rows."""
        var graph = Graph()
        if len(rows) == 0:
            return graph^
        if len(names) < len(rows):
            raise Error(t"an MSA of {len(rows)} rows and {len(names)} names")
        var width = rows[0].byte_length()
        for r in range(len(rows)):
            if rows[r].byte_length() != width:
                raise Error(t"MSA row {names[r]} has {rows[r].byte_length()} columns, the first {width}")
        var at = List[Int](length=256, fill=-1)
        for c in range(width):
            var column = graph.column_count()
            graph.column_first.append(-1)
            for b in range(256):
                at[b] = -1
            for r in range(len(rows)):
                if r == skip:
                    continue
                var letter = rows[r].as_bytes()[c]
                if letter != GAP and at[Int(letter)] < 0:
                    at[Int(letter)] = graph._node(letter, column)
        for r in range(len(rows)):
            if r == skip:
                continue
            var path = List[Int32]()
            var bytes = rows[r].as_bytes()
            for c in range(width):
                if bytes[c] == GAP:
                    continue
                var landed = Int(graph.column_first[c])
                while landed >= 0 and graph.bases[landed] != bytes[c]:
                    landed = Int(graph.column_next[landed])
                graph._on_path(path, landed)
            graph.paths.append(path^)
            graph.names.append(names[r])
        return graph^

    def heaviest_path(self, order: List[Int]) -> List[Int]:
        """abPOA's heaviest bundle over every path: from a virtual source to a virtual sink, each node followed by its
        heaviest out-edge (an edge weighing the paths that take it), the sink's edge weighing the paths that end there
        (so the consensus stops where most sequences do), ties to the heavier rest of the path. `order` holds the
        nodes in topological order."""
        var count = len(self.bases)
        var starts = List[Int](length=count, fill=0)
        var ends = List[Int](length=count, fill=0)
        for path in self.paths:
            if len(path) > 0:
                starts[Int(path[0])] += 1
                ends[Int(path[len(path) - 1])] += 1
        return self._bundle(order, starts, ends, self.edge_weight)

    def heaviest_path(self, order: List[Int], members: List[Int], spread: Int = 0) -> List[Int]:
        """`heaviest_path` with only the paths `members` weighing the edges, starts and ends; with `spread`, the start
        voted on over that many nodes along the walk (`_spread_vote`)."""
        var count = len(self.bases)
        var starts = List[Int](length=count, fill=0)
        var ends = List[Int](length=count, fill=0)
        var weight = List[Int32](length=len(self.edge_to), fill=0)
        for p in members:
            ref path = self.paths[p]
            if len(path) == 0:
                continue
            starts[Int(path[0])] += 1
            ends[Int(path[len(path) - 1])] += 1
            for k in range(len(path) - 1):
                weight[self.edge(Int(path[k]), Int(path[k + 1]))] += 1
        return self._bundle(order, starts, ends, weight, spread=spread)

    def heaviest_path_spanning(self, order: List[Int], spanning: List[Int], bounded: Bool = True) -> List[Int]:
        """`heaviest_path` from where the paths `spanning` start to where they end, every path weighing the edges:
        the reads anchored at one end only add their votes inside without ending the consensus where they stop.
        `bounded`: never past the last of the paths `spanning` end at, however many reads go on (a read fixed at one
        anchor may run past the other: a local assembler cuts it to the gap's estimated length and a margin);
        otherwise the reads going on past their ends carry the walk on (an assembly from one anchor alone). Where the
        votes tie, the paths `spanning` break the tie: a read on both anchors is the surer of a local assembly (a read
        fixed at one may be a paralog's, sharing that anchor alone); so every edge weighs its paths times one more than
        there are paths, plus its paths of `spanning`."""
        var count = len(self.bases)
        var starts = List[Int](length=count, fill=0)
        var ends = List[Int](length=count, fill=0)
        var scale = len(self.paths) + 1
        # 64 bits: the paths squared pass 32 bits at some 46,000 paths.
        var weight = List[Int64](capacity=len(self.edge_weight))
        for w in self.edge_weight:
            weight.append(Int64(w) * Int64(scale))
        for p in spanning:
            ref path = self.paths[p]
            if len(path) > 0:
                starts[Int(path[0])] += 1
                ends[Int(path[len(path) - 1])] += scale + 1
            for k in range(len(path) - 1):
                weight[self.edge(Int(path[k]), Int(path[k + 1]))] += 1
        return self._bundle(order, starts, ends, weight, bounded)

    def heaviest_path_backward(self, order: List[Int], members: List[Int], spread: Int = 0) -> List[Int]:
        """`heaviest_path` over the paths `members` walked from the end: from a virtual sink to a virtual source,
        each node preceded by its heaviest in-edge, the source's edge weighing the members that start there (so the
        consensus starts where most of them do). For reads anchored at the graph's end alone; with `spread`, the end
        voted on over that many nodes along the walk (`_spread_vote`)."""
        var count = len(self.bases)
        var starts = List[Int](length=count, fill=0)
        var ends = List[Int](length=count, fill=0)
        var weight = List[Int32](length=len(self.edge_to), fill=0)
        for p in members:
            ref path = self.paths[p]
            if len(path) == 0:
                continue
            starts[Int(path[0])] += 1
            ends[Int(path[len(path) - 1])] += 1
            for k in range(len(path) - 1):
                weight[self.edge(Int(path[k]), Int(path[k + 1]))] += 1
        var score = List[Int](length=count, fill=0)
        var before = List[Int](length=count, fill=-1)
        for node in order:
            var best_weight = starts[node] if starts[node] > 0 else -1
            var best_score = 0
            var choice = -1
            var e = Int(self.first_in[node])
            while e >= 0:
                var p = Int(self.edge_from[e])
                var w = Int(weight[e])
                if w > 0 and (w > best_weight or (w == best_weight and score[p] >= best_score)):
                    best_weight = w
                    best_score = score[p]
                    choice = p
                e = Int(self.next_in[e])
            before[node] = choice
            score[node] = max(best_weight, 0) + (score[choice] if choice >= 0 else 0)
        var at = -1
        if spread > 0:
            at = _spread_vote(order, ends, score, before, spread, True)
        else:
            for k in reversed(range(len(order))):
                var node = order[k]
                if ends[node] > 0 and (
                    at < 0 or ends[node] > ends[at] or (ends[node] == ends[at] and score[node] > score[at])
                ):
                    at = node
        var path = List[Int]()
        while at >= 0:
            path.append(at)
            at = before[at]
        path.reverse()
        return path^

    def _bundle[
        dt: DType
    ](
        self,
        order: List[Int],
        starts: List[Int],
        ends: List[Int],
        weight: List[Scalar[dt]],
        bounded: Bool = False,
        spread: Int = 0,
    ) -> List[Int]:
        var count = len(self.bases)
        var score = List[Int](length=count, fill=0)
        var next = List[Int](length=count, fill=-1)
        # `bounded`: only nodes an end can still be reached from (all of them otherwise).
        var reaches = List[Bool](length=count, fill=not bounded)
        for k in reversed(range(len(order))):
            var node = order[k]
            var best_weight = ends[node] if ends[node] > 0 else -1
            var best_score = 0
            var choice = -1
            if ends[node] > 0:
                reaches[node] = True
            var e = Int(self.first_out[node])
            while e >= 0:
                var s = Int(self.edge_to[e])
                var w = Int(weight[e])
                if reaches[s]:
                    reaches[node] = True
                else:
                    w = 0
                if w > 0 and (w > best_weight or (w == best_weight and score[s] >= best_score)):
                    best_weight = w
                    best_score = score[s]
                    choice = s
                e = Int(self.next_out[e])
            next[node] = choice
            score[node] = max(best_weight, 0) + (score[choice] if choice >= 0 else 0)
        var at = -1
        if spread > 0:
            at = _spread_vote(order, starts, score, next, spread, False)
        else:
            for node in order:
                if starts[node] > 0 and (
                    at < 0 or starts[node] > starts[at] or (starts[node] == starts[at] and score[node] > score[at])
                ):
                    at = node
        var path = List[Int]()
        while at >= 0:
            path.append(at)
            at = next[at]
        return path^

    def column_order(self) raises -> List[Int]:
        """Each column's place left to right: a topological order of the columns (Kahn's)."""
        var count = self.column_count()
        var incoming = List[Int32](length=count, fill=0)
        for e in range(len(self.edge_to)):
            incoming[Int(self.columns[Int(self.edge_to[e])])] += 1
        var ready = List[Int]()
        for c in range(count):
            if incoming[c] == 0:
                ready.append(c)
        var place = List[Int](length=count, fill=-1)
        var next = 0
        var head = 0
        while head < len(ready):
            var c = ready[head]
            head += 1
            place[c] = next
            next += 1
            var m = Int(self.column_first[c])
            while m >= 0:
                var e = Int(self.first_out[m])
                while e >= 0:
                    var b = Int(self.columns[Int(self.edge_to[e])])
                    incoming[b] -= 1
                    if incoming[b] == 0:
                        ready.append(b)
                    e = Int(self.next_out[e])
                m = Int(self.column_next[m])
        if next != count:
            raise Error(t"the columns form a cycle: {count - next} of {count} unplaced")
        return place^

    def msa(self) raises -> List[String]:
        """Every path as a row of one length, `-` in the columns it skips."""
        var place = self.column_order()
        var rows = List[String]()
        for path in self.paths:
            var row = List[UInt8](length=self.column_count(), fill=GAP)
            for node in path:
                row[place[Int(self.columns[Int(node)])]] = self.bases[Int(node)]
            rows.append(String(unsafe_from_utf8=row))
        return rows^


struct Runs(Movable):
    """A CIGAR's runs: their letters and lengths."""

    var ops: List[UInt8]
    var counts: List[Int]

    def __init__(out self):
        self.ops = List[UInt8]()
        self.counts = List[Int]()


def cigar_runs(cigar: StringSlice) raises -> Runs:
    """`cigar`'s runs."""
    var runs = Runs()
    var length = 0
    for c in cigar.as_bytes():
        if c >= UInt8(ord("0")) and c <= UInt8(ord("9")):
            length = length * 10 + Int(c - UInt8(ord("0")))
            continue
        if length == 0:
            raise Error(t"run without a length in {cigar}")
        runs.ops.append(c)
        runs.counts.append(length)
        length = 0
    return runs^


def _spread_vote(
    order: List[Int], votes: List[Int], score: List[Int], step: List[Int], spread: Int, backward: Bool
) -> Int:
    """Where a heaviest walk starts (or, `backward`, ends) when the reads' own starts (`votes`, per node) spread over
    a few nodes: each node a read starts at heads a run of `spread` nodes along the walk (`step`: each node's next,
    or its previous `backward`), the run holding the most starts wins (of two as many, the heavier walk, `score`),
    and the walk starts at its median start. Node by node, a few reads sharing one start outvoted many whose
    anchored ends differ by a letter or a few (where they were cut, or their errors there)."""
    var best = -1
    var best_votes = 0
    var candidates = order.copy()
    if backward:
        candidates.reverse()
    for node in candidates:
        if votes[node] == 0:
            continue
        var total = 0
        var at = node
        var k = 0
        while at >= 0 and k < spread:
            total += votes[at]
            at = step[at]
            k += 1
        if best < 0 or total > best_votes or (total == best_votes and score[node] > score[best]):
            best = node
            best_votes = total
    if best < 0:
        return -1
    var seen = 0
    var at = best
    while at >= 0:
        seen += votes[at]
        if 2 * seen >= best_votes:
            return at
        at = step[at]
    return best
