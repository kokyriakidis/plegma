"""E6 start: a progressive POA of one region, every sequence anchored to and aligned with one partner path.

usage: progressive <in.fa> <out.msa.fa> <config> <order> <costs> <partner> [path|graph]

The distinct sequences join longest first (`longest`) or in file order (`given`). Each is anchored against a
partner path, aligned to it through its backbone and fused into the graph along it. The partner is the first
sequence (`first`, a star) or the path already in the graph nearest to it in edit distance (`nearest`).
Every record of the input gets its sequence's row, so duplicates share one.
- `<config>`: `none` (karan end to end, no anchors) or `mums` (MUMs of 20 letters, ChainX costs, recursing
  until no MUM is left).
- `<costs>`: karan's costs in the windows: `edit`, `affine` (4, 6, 2) or `convex` (4, 6, 2, 24, 1).
- `path` (the default) aligns each window to the partner path with karan; `graph` aligns it to the graph between
  the anchors' nodes (gap-affine 4, 6, 2), to the path only past 2e8 cells.
"""

from std.collections import Dict
from std.sys import argv
from std.time import perf_counter_ns

from karan import Aligner, Costs
from sporos import read_fasta

from plegma.align import anchored_alignment, graph_alignment
from plegma.poa import GapAffine
from plegma.backbone import Level, MUMS
from plegma.chain import ChainCosts
from plegma.graph import Graph


def main() raises:
    var args = argv()
    if len(args) < 7:
        raise Error(
            "usage: progressive <in.fa> <out.msa.fa> <none|mums> <longest|given> <edit|affine|convex> <first|nearest>"
        )
    var config = String(args[3])
    var order = String(args[4])
    var cost_name = String(args[5])
    var partner_rule = String(args[6])
    var costs = Costs.edit()
    if cost_name == "affine":
        costs = Costs.affine(4, 6, 2)
    elif cost_name == "convex":
        costs = Costs.two_piece(4, 6, 2, 24, 1)
    elif cost_name != "edit":
        raise Error(t"unknown costs {cost_name}")
    var fill = String(args[7]) if len(args) > 7 else String("path")
    if partner_rule != "first" and partner_rule != "nearest":
        raise Error(t"unknown partner rule {partner_rule}")
    var levels = List[Level]()
    if config == "mums":
        levels.append(Level(MUMS, 20))
    elif config != "none":
        raise Error(t"unknown config {config}")
    var records = read_fasta(String(args[1]))
    var distinct = Dict[String, Int]()
    var firsts = List[Int]()
    var row_of = List[Int]()
    for i in range(len(records)):
        var key = String(unsafe_from_utf8=List[UInt8](records.seq(i)))
        if key not in distinct:
            distinct[key] = len(firsts)
            firsts.append(i)
        row_of.append(distinct[key])
    var ranked = List[Int]()
    for d in range(len(firsts)):
        ranked.append(d)
    if order == "longest":
        # Insertion sort, longest first, file order among equals.
        for a in range(1, len(ranked)):
            var at = a
            while at > 0 and len(records.seq(firsts[ranked[at - 1]])) < len(records.seq(firsts[ranked[at]])):
                ranked.swap_elements(at - 1, at)
                at -= 1
    var t0 = perf_counter_ns()
    var graph = Graph()
    var aligner = Aligner()
    var path_of = List[Int](length=len(firsts), fill=-1)
    var texts = List[String]()
    for d in range(len(firsts)):
        texts.append(String(unsafe_from_utf8=List[UInt8](records.seq(firsts[d]))))
    var start_name = records.names[firsts[ranked[0]]]
    path_of[ranked[0]] = graph.add(start_name, records.seq(firsts[ranked[0]]))
    var total = 0
    for k in range(1, len(ranked)):
        var d = ranked[k]
        var i = firsts[d]
        var partner = ranked[0]
        if partner_rule == "nearest":
            var best = aligner.distance(texts[partner], texts[d])
            for e in range(1, k):
                var near = aligner.distance(texts[ranked[e]], texts[d], max_cost=best - 1)
                if near:
                    best = near.value()
                    partner = ranked[e]
        var name = records.names[i]
        if fill == "graph":
            var aligned = graph_alignment(
                aligner,
                graph,
                path_of[partner],
                records.seq(firsts[partner]),
                records.seq(i),
                levels,
                ChainCosts.chainx(),
                GapAffine(4, 6, 2),
                200_000_000,
                stay=True,
            )
            path_of[d] = graph.fuse_nodes(name, records.seq(i), aligned)
            continue
        var found = anchored_alignment(
            aligner, records.seq(firsts[partner]), records.seq(i), levels, ChainCosts.chainx(), costs, stay=True
        )
        total += found.cost
        path_of[d] = graph.fuse(name, records.seq(i), path_of[partner], found.ops, found.counts)
    var rows = graph.msa()
    var out = open(String(args[2]), "w")
    for i in range(len(records)):
        out.write(">", records.names[i], "\n", rows[path_of[row_of[i]]], "\n")
    out.close()
    var seconds = Float64(perf_counter_ns() - t0) / 1e9
    print(
        t"{args[1]}\t{config}\t{order}\t{cost_name}\t{partner_rule}\t{fill}\t{len(firsts)}\t{len(graph)}\t{rows[0].byte_length()}\t{total}\t{seconds}"
    )
