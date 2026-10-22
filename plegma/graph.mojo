"""The partial order graph: one node per letter, nodes aligned to each other grouped into columns, and every
sequence's path through it kept.

A sequence joins by an alignment to one path already in the graph, as runs of `=`, `X`, `I` and `D` (karan's
CIGAR letters, the path the reference): a match lands on the path's node, a substitution on the node of that letter
in the path node's column (made if missing), an insertion on a node of its own column.
"""

from std.collections import Dict


comptime MATCH = UInt8(ord("="))
comptime MISMATCH = UInt8(ord("X"))
comptime INSERTION = UInt8(ord("I"))
comptime DELETION = UInt8(ord("D"))
comptime GAP = UInt8(ord("-"))


struct Graph(Movable, Sized):
    """Nodes (`bases`, `columns`), weighted edges and the paths of the sequences added so far."""

    var bases: List[UInt8]
    """Each node's letter, as the sequence spelled it."""
    var columns: List[Int]
    """Each node's column."""
    var members: List[List[Int]]
    """Each column's nodes, one per letter."""
    var successors: List[List[Int]]
    var predecessors: List[List[Int]]
    var weights: List[List[Int]]
    """How many paths take each edge, beside `successors`."""
    var paths: List[List[Int]]
    var names: List[String]

    def __init__(out self):
        self.bases = List[UInt8]()
        self.columns = List[Int]()
        self.members = List[List[Int]]()
        self.successors = List[List[Int]]()
        self.predecessors = List[List[Int]]()
        self.weights = List[List[Int]]()
        self.paths = List[List[Int]]()
        self.names = List[String]()

    def __len__(self) -> Int:
        """The number of nodes."""
        return len(self.bases)

    def _node(mut self, base: UInt8, column: Int) -> Int:
        var node = len(self.bases)
        self.bases.append(base)
        var at = column
        if at < 0:
            at = len(self.members)
            self.members.append(List[Int]())
        self.columns.append(at)
        self.members[at].append(node)
        self.successors.append(List[Int]())
        self.predecessors.append(List[Int]())
        self.weights.append(List[Int]())
        return node

    def _link(mut self, before: Int, after: Int):
        for i in range(len(self.successors[before])):
            if self.successors[before][i] == after:
                self.weights[before][i] += 1
                return
        self.successors[before].append(after)
        self.predecessors[after].append(before)
        self.weights[before].append(1)

    def _on_path(mut self, mut path: List[Int], node: Int):
        if len(path) > 0:
            self._link(path[len(path) - 1], node)
        path.append(node)

    def add(mut self, name: String, sequence: Span[UInt8, _]) -> Int:
        """Adds `sequence` as a path of new nodes, each its own column; returns the path's index."""
        var path = List[Int](capacity=len(sequence))
        for c in sequence:
            self._on_path(path, self._node(c, -1))
        self.paths.append(path^)
        self.names.append(name)
        return len(self.paths) - 1

    def fuse(
        mut self, name: String, sequence: Span[UInt8, _], partner: Int, ops: List[UInt8], counts: List[Int]
    ) raises -> Int:
        """Adds `sequence` aligned to path `partner` by the runs `ops[i]` x `counts[i]`, which must span both;
        returns the new path's index."""
        var path = List[Int](capacity=len(sequence))
        var along = self.paths[partner].copy()
        var i = 0
        var j = 0
        for r in range(len(ops)):
            var op = ops[r]
            for _ in range(counts[r]):
                if op == DELETION:
                    i += 1
                    continue
                if op == INSERTION:
                    self._on_path(path, self._node(sequence[j], -1))
                    j += 1
                    continue
                if i >= len(along) or j >= len(sequence):
                    raise Error(t"alignment of {name} runs past its sequences")
                self._on_path(path, self._lands(along[i], sequence[j]))
                i += 1
                j += 1
        if i != len(along) or j != len(sequence):
            raise Error(t"alignment of {name} covers {i} of {len(along)} and {j} of {len(sequence)} letters")
        self.paths.append(path^)
        self.names.append(name)
        return len(self.paths) - 1

    def _lands(mut self, node: Int, base: UInt8) -> Int:
        """The node a letter aligned to `node` lands on: `node` itself, or the node of that letter in its
        column, made if missing."""
        if self.bases[node] == base:
            return node
        var column = self.columns[node]
        for m in self.members[column]:
            if self.bases[m] == base:
                return m
        return self._node(base, column)

    def fuse_nodes(mut self, name: String, sequence: Span[UInt8, _], aligned: List[Int]) -> Int:
        """Adds `sequence`, its letter `j` aligned to node `aligned[j]` (-1: inserted); returns the new path's
        index."""
        var path = List[Int](capacity=len(sequence))
        for j in range(len(sequence)):
            var node = aligned[j]
            self._on_path(path, self._lands(node, sequence[j]) if node >= 0 else self._node(sequence[j], -1))
        self.paths.append(path^)
        self.names.append(name)
        return len(self.paths) - 1

    def node_order(self) raises -> List[Int]:
        """Each node's rank in a topological order of the nodes (Kahn's, lowest node first among the ready)."""
        var count = len(self.bases)
        var incoming = List[Int](length=count, fill=0)
        for node in range(count):
            incoming[node] = len(self.predecessors[node])
        var ready = List[Int]()
        for node in range(count):
            if incoming[node] == 0:
                ready.append(node)
        var rank = List[Int](length=count, fill=-1)
        var next = 0
        var head = 0
        while head < len(ready):
            var node = ready[head]
            head += 1
            rank[node] = next
            next += 1
            for s in self.successors[node]:
                incoming[s] -= 1
                if incoming[s] == 0:
                    ready.append(s)
        if next != count:
            raise Error(t"the nodes form a cycle: {count - next} of {count} unplaced")
        return rank^

    @staticmethod
    def from_msa(names: List[String], rows: List[String], skip: Int = -1) -> Graph:
        """The graph of an MSA's rows (all but row `skip`): in each column, one node per letter there, the rows'
        paths through them."""
        var graph = Graph()
        if len(rows) == 0:
            return graph^
        var width = rows[0].byte_length()
        var at = List[Int](length=256, fill=-1)
        for c in range(width):
            var column = len(graph.members)
            graph.members.append(List[Int]())
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
            var path = List[Int]()
            var bytes = rows[r].as_bytes()
            for c in range(width):
                if bytes[c] == GAP:
                    continue
                var landed = -1
                for m in graph.members[c]:
                    if graph.bases[m] == bytes[c]:
                        landed = m
                graph._on_path(path, landed)
            graph.paths.append(path^)
            graph.names.append(names[r])
        return graph^

    def heaviest_path(self, order: List[Int]) -> List[Int]:
        """abPOA's heaviest bundle: from a virtual source to a virtual sink, each node followed by its heaviest
        out-edge, the sink's edge weighing the paths that end there (so the consensus stops where most sequences
        do), ties to the heavier rest of the path. `order` holds the nodes in topological order."""
        var count = len(self.bases)
        var starts = List[Int](length=count, fill=0)
        var ends = List[Int](length=count, fill=0)
        for path in self.paths:
            if len(path) > 0:
                starts[path[0]] += 1
                ends[path[len(path) - 1]] += 1
        var score = List[Int](length=count, fill=0)
        var next = List[Int](length=count, fill=-1)
        for k in reversed(range(len(order))):
            var node = order[k]
            var best_weight = ends[node] if ends[node] > 0 else -1
            var best_score = 0
            var choice = -1
            for i in range(len(self.successors[node])):
                var s = self.successors[node][i]
                var weight = self.weights[node][i]
                if weight > best_weight or (weight == best_weight and score[s] >= best_score):
                    best_weight = weight
                    best_score = score[s]
                    choice = s
            next[node] = choice
            score[node] = max(best_weight, 0) + (score[choice] if choice >= 0 else 0)
        var at = -1
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
        var count = len(self.members)
        var incoming = List[Int](length=count, fill=0)
        var after = List[List[Int]](length=count, fill=List[Int]())
        for node in range(len(self.bases)):
            for s in self.successors[node]:
                var a = self.columns[node]
                var b = self.columns[s]
                after[a].append(b)
                incoming[b] += 1
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
            for b in after[c]:
                incoming[b] -= 1
                if incoming[b] == 0:
                    ready.append(b)
        if next != count:
            raise Error(t"the columns form a cycle: {count - next} of {count} unplaced")
        return place^

    def msa(self) raises -> List[String]:
        """Every path as a row of one length, `-` in the columns it skips."""
        var place = self.column_order()
        var rows = List[String]()
        for path in self.paths:
            var row = List[UInt8](length=len(self.members), fill=GAP)
            for node in path:
                row[place[self.columns[node]]] = self.bases[node]
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
