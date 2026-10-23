"""A pairwise alignment through a backbone: the anchors matched, each window between them aligned exactly by karan.
"""

from karan import Aligner, Costs, Mode, Ties

from .anchors import Anchor
from .backbone import Level, backbone, windows
from .chain import ChainCosts
from .graph import DELETION, Graph, INSERTION, MATCH, cigar_runs
from .poa import GapAffine, SINK, SOURCE, align_window


struct Anchored(Movable):
    """An alignment as CIGAR runs (`ops` x `counts`), its cost, and the backbone it went through."""

    var ops: List[UInt8]
    var counts: List[Int]
    var cost: Int
    var anchors: List[Anchor]

    def __init__(out self):
        self.ops = List[UInt8]()
        self.counts = List[Int]()
        self.cost = 0
        self.anchors = List[Anchor]()

    def _run(mut self, op: UInt8, count: Int):
        if count == 0:
            return
        var last = len(self.ops) - 1
        if last >= 0 and self.ops[last] == op:
            self.counts[last] += count
        else:
            self.ops.append(op)
            self.counts.append(count)


def anchored_alignment(
    mut aligner: Aligner,
    partner: ImmSpan[UInt8, _],
    query: ImmSpan[UInt8, _],
    levels: List[Level],
    chain_costs: ChainCosts,
    costs: Costs = Costs.edit(),
    budget: Int = 0,
    stay: Bool = False,
    repeat_k: Int = 0,
    repeat_distance: Int = 0,
    mask_depth: Int = Int.MAX,
) raises -> Anchored:
    """`query` aligned to `partner` end to end through the backbone `levels` find (its repeat mask as `backbone`
    takes it), every window filled with an optimal alignment under `costs`, indels placed left."""
    var out = Anchored()
    out.anchors = backbone(
        partner,
        query,
        levels,
        chain_costs,
        budget,
        stay,
        repeat_k=repeat_k,
        repeat_distance=repeat_distance,
        mask_depth=mask_depth,
    )
    var gaps = windows(out.anchors, len(partner), len(query))
    for w in range(len(gaps)):
        var window = gaps[w]
        var ref_len = window.ref_end - window.ref_start
        var query_len = window.query_end - window.query_start
        if ref_len + query_len > 0:
            var found = aligner.align(
                partner[window.ref_start : window.ref_end],
                query[window.query_start : window.query_end],
                costs,
                Mode.GLOBAL,
                ties=Ties.LEFT,
            )
            out.cost += found.cost
            var runs = cigar_runs(found.cigar)
            for r in range(len(runs.ops)):
                out._run(runs.ops[r], runs.counts[r])
        if w < len(out.anchors):
            out._run(MATCH, out.anchors[w].length)
    return out^


def graph_alignment(
    mut aligner: Aligner,
    graph: Graph,
    partner_path: Int,
    partner: ImmSpan[UInt8, _],
    query: ImmSpan[UInt8, _],
    levels: List[Level],
    chain_costs: ChainCosts,
    costs: GapAffine,
    max_cells: Int,
    stay: Bool = False,
    budget: Int = 0,
) raises -> List[Int]:
    """For each letter of `query`, the graph node it aligns to (-1: inserted): the backbone against path
    `partner_path` (whose sequence is `partner`) pins its anchors to that path's nodes, and each window between
    them is aligned to the graph between the anchors' nodes. A window over `max_cells` is aligned to the
    partner path alone, by karan."""
    var path = graph.paths[partner_path].copy()
    var rank = graph.node_order()
    var order = List[Int](length=len(rank), fill=0)
    for node in range(len(rank)):
        order[rank[node]] = node
    var anchors = backbone(partner, query, levels, chain_costs, budget, stay)
    var gaps = windows(anchors, len(partner), len(query))
    var aligned = List[Int](length=len(query), fill=-1)
    for w in range(len(gaps)):
        var window = gaps[w]
        var left = Int(path[window.ref_start - 1]) if window.ref_start > 0 else SOURCE
        var right = Int(path[window.ref_end]) if window.ref_end < len(path) else SINK
        if window.query_end > window.query_start:
            var found = align_window(
                graph, rank, order, left, right, query[window.query_start : window.query_end], costs, max_cells
            )
            if found:
                ref nodes = found.value().nodes
                for j in range(len(nodes)):
                    aligned[window.query_start + j] = nodes[j]
            else:
                var along = aligner.align(
                    partner[window.ref_start : window.ref_end],
                    query[window.query_start : window.query_end],
                    Costs.affine(Int(costs.mismatch), Int(costs.open), Int(costs.extend)),
                    Mode.GLOBAL,
                    ties=Ties.LEFT,
                )
                var runs = cigar_runs(along.cigar)
                var i = window.ref_start
                var j = window.query_start
                for r in range(len(runs.ops)):
                    for _ in range(runs.counts[r]):
                        var op = runs.ops[r]
                        if op == DELETION:
                            i += 1
                        elif op == INSERTION:
                            j += 1
                        else:
                            aligned[j] = Int(path[i])
                            i += 1
                            j += 1
        if w < len(anchors):
            var a = anchors[w]
            for k in range(a.length):
                aligned[a.query_pos + k] = Int(path[a.ref_pos + k])
    return aligned^
