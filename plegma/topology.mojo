"""Graph topological order and heaviest remaining path lengths (`Topology`), with reusable in-place computation."""

from .graph import Graph


struct Topology(Movable):
    """The graph's nodes in topological order, each node's rank, and the length of its remaining path along
    the heaviest out-edges (abPOA's band diagonal); or the same over a part of the graph (`allowed`)."""

    var order: List[Int]
    var rank: List[Int]
    var remaining: List[Int]

    def __init__(out self):
        self.order = List[Int]()
        self.rank = List[Int]()
        self.remaining = List[Int]()

    def compute_order(mut self, graph: Graph, mut incoming: List[Int32], mut ready: List[Int]) raises:
        """Recomputes topological `order` and `rank` in-place without computing `remaining`."""
        graph.node_order_into(self.rank, incoming, ready)
        var n = len(self.rank)
        self.order.resize(n, 0)
        for node in range(n):
            self.order[self.rank[node]] = node

    def compute(mut self, graph: Graph, mut incoming: List[Int32], mut ready: List[Int]) raises:
        """Recomputes topological order, rank, and remaining path lengths in-place."""
        self.compute_order(graph, incoming, ready)
        self.recompute_remaining(graph)

    def recompute_remaining(mut self, graph: Graph):
        """Recomputes heaviest remaining path lengths in-place, reusing the existing full `self.order`."""
        var n = len(self.order)
        self.remaining.resize(n, 0)
        for r in reversed(range(n)):
            var node = self.order[r]
            var rem = 0
            var best_weight = Int32(-1)
            var e = Int(graph.first_out[node])
            while e >= 0:
                var s = Int(graph.edge_to[e])
                var weight = graph.edge_weight[e]
                if weight > best_weight or (weight == best_weight and self.remaining[s] + 1 > rem):
                    best_weight = weight
                    rem = self.remaining[s] + 1
                e = Int(graph.next_out[e])
            self.remaining[node] = rem

    def compute_allowed(
        mut self,
        graph: Graph,
        allowed: List[Bool],
        mut incoming: List[Int32],
        mut ready: List[Int],
        mut full: List[Int],
    ) raises:
        """Recomputes topological order over `allowed` nodes in-place."""
        graph.node_order_into(self.rank, incoming, ready)
        var n = len(self.rank)
        full.resize(n, 0)
        for node in range(n):
            full[self.rank[node]] = node
        self.order.clear()
        for node in full:
            if allowed[node]:
                self.order.append(node)
        self.remaining.resize(n, 0)
        for node in range(n):
            self.remaining[node] = 0
        for r in reversed(range(len(self.order))):
            var node = self.order[r]
            var rem = 0
            var best_weight = Int32(-1)
            var e = Int(graph.first_out[node])
            while e >= 0:
                var s = Int(graph.edge_to[e])
                var weight = graph.edge_weight[e]
                if allowed[s] and (weight > best_weight or (weight == best_weight and self.remaining[s] + 1 > rem)):
                    best_weight = weight
                    rem = self.remaining[s] + 1
                e = Int(graph.next_out[e])
            self.remaining[node] = rem

    def __init__(out self, graph: Graph) raises:
        self.order = List[Int]()
        self.rank = List[Int]()
        self.remaining = List[Int]()
        var incoming = List[Int32]()
        var ready = List[Int]()
        self.compute(graph, incoming, ready)

    def __init__(out self, graph: Graph, allowed: List[Bool]) raises:
        """Over the nodes `allowed` alone: they in topological order (`order`), and each one's remaining path
        through them alone; every node's rank in the whole graph."""
        self.order = List[Int]()
        self.rank = List[Int]()
        self.remaining = List[Int]()
        var incoming = List[Int32]()
        var ready = List[Int]()
        var full = List[Int]()
        self.compute_allowed(graph, allowed, incoming, ready, full)
