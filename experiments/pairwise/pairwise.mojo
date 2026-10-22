"""E1-E5: pairwise anchors, chains and filled backbones on one vntr-eval region.

usage: pairwise <hap32.fa> <region> <partner name prefix> <flank left> <flank right> <period> <out prefix>

Each distinct sequence other than the partner is a query. `<out prefix>.census.tsv` holds E1: per anchor kind
and minimum length, how many anchors fall in the repeat array and how many lie on some optimal alignment (unit
costs). `<out prefix>.chains.tsv` holds E2-E5: per backbone configuration, the chain, the windows it leaves and
the unit cost of filling them exactly, against the pair's optimum.
"""

from std.collections import Dict
from std.sys import argv
from std.time import perf_counter_ns

from karan import Aligner
from sporos import read_fasta

from plegma.anchors import Anchor, find_anchors, mums, rare, smems
from plegma.backbone import ALL, Level, MUMS, RARE, REFERENCE_UNIQUE, SMEMS, backbone, windows
from plegma.chain import ChainCosts, linked_cost


comptime SAMPLE = 200
"""Anchors of a kind tested against the optimum, at most, per pair and length."""


def _text(bytes: Span[UInt8, _]) -> String:
    return String(unsafe_from_utf8=List[UInt8](bytes))


struct Config(Copyable, Movable):
    var name: String
    var levels: List[Level]
    var costs: ChainCosts
    var budget: Int
    var repeat_last: Bool

    def __init__(
        out self, name: String, var levels: List[Level], costs: ChainCosts, budget: Int = 0, repeat_last: Bool = False
    ):
        self.name = name
        self.levels = levels^
        self.costs = costs
        self.budget = budget
        self.repeat_last = repeat_last


def _configs(period: Int) -> List[Config]:
    var out = List[Config]()
    var concave = ChainCosts(open1=2.0, extend1=1.0, open2=20.0, extend2=0.1)
    var periodic = ChainCosts(open1=2.0, extend1=1.0, open2=20.0, extend2=0.1, period=period, off_period=1.0)
    var weighted = ChainCosts(gain=1.0, rarity=True)
    for min_len in [20, 50, 100]:
        out.append(Config(String(t"heaviest.mum.L{min_len}"), [Level(MUMS, min_len)], ChainCosts.heaviest()))
        out.append(Config(String(t"chainx.mum.L{min_len}"), [Level(MUMS, min_len)], ChainCosts.chainx()))
        out.append(Config(String(t"chainx.mum.L{min_len}.rec"), [Level(MUMS, min_len)], ChainCosts.chainx(), 0, True))
    out.append(Config("heaviest.mum.L20.rec", [Level(MUMS, 20)], ChainCosts.heaviest(), 0, True))
    out.append(Config("chainx.mumref.L20", [Level(REFERENCE_UNIQUE, 20)], ChainCosts.chainx()))
    out.append(Config("chainx.rare8.L20", [Level(RARE, 20, 8)], ChainCosts.chainx()))
    out.append(Config("chainx.rare50.L31", [Level(RARE, 31, 50)], ChainCosts.chainx()))
    out.append(Config("chainx.smem.L20", [Level(SMEMS, 20)], ChainCosts.chainx()))
    out.append(Config("chainx.all.L31", [Level(ALL, 31)], ChainCosts.chainx()))
    out.append(Config("concave.mum.L20", [Level(MUMS, 20)], concave))
    out.append(Config("concave.rare8.L20", [Level(RARE, 20, 8)], concave))
    out.append(Config("periodic.mum.L20", [Level(MUMS, 20)], periodic))
    out.append(Config("periodic.rare8.L20", [Level(RARE, 20, 8)], periodic))
    out.append(Config("weighted.rare8.L20", [Level(RARE, 20, 8)], weighted))
    var schedule: List[Level] = [Level(MUMS, 51), Level(MUMS, 31), Level(MUMS, 21), Level(MUMS, 15)]
    for budget in [1000, 2000, 5000]:
        out.append(Config(String(t"chainx.mum.51-15.W{budget}"), schedule.copy(), ChainCosts.chainx(), budget))
    out.append(Config("chainx.flank.L50", [Level(MUMS, 50)], ChainCosts.chainx()))
    return out^


def _in_array(a: Anchor, partner_len: Int, flank_left: Int, flank_right: Int) -> Bool:
    var middle = a.ref_pos + a.length // 2
    return middle >= flank_left and middle < partner_len - flank_right


def _optimal_through(mut aligner: Aligner, partner: String, query: String, a: Anchor, best: Int) raises -> Bool:
    """Whether some optimal alignment matches `a` whole: the costs before and after it add up to the optimum."""
    var before = aligner.distance(
        String(partner[byte = 0 : a.ref_pos]), String(query[byte = 0 : a.query_pos]), max_cost=best
    )
    if not before:
        return False
    var after = aligner.distance(
        String(partner[byte = a.ref_end() :]), String(query[byte = a.query_end() :]), max_cost=best - before.value()
    )
    return Bool(after) and before.value() + after.value() == best


def main() raises:
    var args = argv()
    if len(args) != 8:
        raise Error("usage: pairwise <hap32.fa> <region> <partner prefix> <flank left> <flank right> <period> <out>")
    var region = String(args[2])
    var prefix = String(args[3])
    var flank_left = Int(args[4])
    var flank_right = Int(args[5])
    var period = Int(args[6])
    var out_prefix = String(args[7])
    var records = read_fasta(String(args[1]))
    var partner_index = -1
    var seen = Dict[String, Int]()
    var queries = List[Int]()
    for i in range(len(records)):
        if partner_index < 0 and records.names[i].startswith(prefix):
            partner_index = i
    if partner_index < 0:
        raise Error(t"no record named {prefix}* in {args[1]}")
    seen[_text(records.seq(partner_index))] = partner_index
    for i in range(len(records)):
        var key = _text(records.seq(i))
        if key not in seen:
            seen[key] = i
            queries.append(i)

    var aligner = Aligner()
    var configs = _configs(period)
    var partner_bytes = records.seq(partner_index)
    var partner = _text(partner_bytes)
    var census = open(out_prefix + ".census.tsv", "w")
    census.write("region\tquery\tquery_len\tpartner_len\toptimum\tmin_len\tkind\tanchors\tin_array\ttested")
    census.write("\ttested_in_array\toptimal\toptimal_in_array\n")
    var chains = open(out_prefix + ".chains.tsv", "w")
    chains.write("region\tquery\tquery_len\tpartner_len\toptimum\toptimum_ms\tconfig\tanchors\tanchored\tchain_cost")
    chains.write("\tfilled\twindows\tlargest_window\twindows_over_1kb\tbackbone_ms\tfill_ms\n")
    for qi in queries:
        var query_bytes = records.seq(qi)
        var query = _text(query_bytes)
        var name = records.names[qi]
        var t0 = perf_counter_ns()
        var best = aligner.distance(partner, query)
        var optimum_ms = Float64(perf_counter_ns() - t0) / 1e6
        var head = String(t"{region}\t{name}\t{len(query_bytes)}\t{len(partner_bytes)}\t{best}")

        # E1: which anchors lie on an optimal alignment.
        for min_len in [15, 21, 31, 51]:
            var every = find_anchors(partner_bytes, query_bytes, min_len, 51)
            var kinds = List[Tuple[String, List[Anchor]]]()
            kinds.append(("mum", mums(every)))
            var reference_unique = List[Anchor]()
            for a in every:
                if a.ref_copies == 1:
                    reference_unique.append(a)
            kinds.append(("mumreference", reference_unique^))
            kinds.append(("rare2", rare(every, 2)))
            kinds.append(("rare8", rare(every, 8)))
            kinds.append(("rare50", rare(every, 50)))
            kinds.append(("smem", smems(every)))
            kinds.append(("all", every.copy()))
            for kind in kinds:
                ref anchors = kind[1]
                var in_array = 0
                for a in anchors:
                    if _in_array(a, len(partner_bytes), flank_left, flank_right):
                        in_array += 1
                var step = max(1, len(anchors) // SAMPLE)
                var tested = 0
                var tested_in_array = 0
                var optimal = 0
                var optimal_in_array = 0
                for i in range(0, len(anchors), step):
                    var a = anchors[i]
                    var inside = _in_array(a, len(partner_bytes), flank_left, flank_right)
                    var on = _optimal_through(aligner, partner, query, a, best)
                    tested += 1
                    tested_in_array += Int(inside)
                    optimal += Int(on)
                    optimal_in_array += Int(on and inside)
                census.write(
                    t"{head}\t{min_len}\t{kind[0]}\t{len(anchors)}\t{in_array}\t{tested}\t{tested_in_array}",
                    t"\t{optimal}\t{optimal_in_array}\n",
                )

        # E2-E5: backbones and their exact fill.
        for config in configs:
            var t1 = perf_counter_ns()
            var picked = backbone(
                partner_bytes, query_bytes, config.levels, config.costs, config.budget, config.repeat_last
            )
            var backbone_ms = Float64(perf_counter_ns() - t1) / 1e6
            var anchored = 0
            for a in picked:
                anchored += a.length
            var chain_cost = linked_cost(picked, len(partner_bytes), len(query_bytes), ChainCosts.chainx())
            var t2 = perf_counter_ns()
            var filled = 0
            var largest = 0
            var over = 0
            var gaps = windows(picked, len(partner_bytes), len(query_bytes))
            for w in gaps:
                filled += aligner.distance(
                    String(partner[byte = w.ref_start : w.ref_end]), String(query[byte = w.query_start : w.query_end])
                )
                largest = max(largest, w.size())
                over += Int(w.size() > 1000)
            var fill_ms = Float64(perf_counter_ns() - t2) / 1e6
            chains.write(
                t"{head}\t{optimum_ms}\t{config.name}\t{len(picked)}\t{anchored}\t{chain_cost}\t{filled}",
                t"\t{len(gaps)}\t{largest}\t{over}\t{backbone_ms}\t{fill_ms}\n",
            )
    census.close()
    chains.close()
