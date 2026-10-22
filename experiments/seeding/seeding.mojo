"""Seeding strategies on pairs of sequences: how close each one's backbone, filled exactly by karan, comes to the
true alignment and to the optimum.

usage: seeding <pairs.fa> <truth.tsv or -> <out.tsv> <strategy> [<strategy> ...]

Pairs are consecutive records (a then b; b is aligned to a). The truth, when given, is a CIGAR of b against a
per pair (`simulate_vntr.py`). A strategy is `levels[;option ...]`:
- `levels`: comma-separated seed levels, each a kind and a minimum length: `mum51` (unique in both),
  `smem31` (supermaximal), `ms21` (MUMs and SMEMs), `ref20` (unique in a), `rare8_21` (at most 8 copies in
  each), `combo8_21` (rare or supermaximal, the command line's -S seeds), `all31`; `none` for no anchors.
- `stay`: seed every gap again at the same level until it yields nothing, then the next level.
- `s=F`: what the chain charges an unanchored letter (0.1); `T=N`: letters each anchor gives back to the DP at
  both ends (10); `n=N`: gaps of at most N letters are not seeded again (0); `H=N`: chain look-back (200).
The fill and the optimum use two-piece gap-affine costs (4, 6, 2, 24, 1), ties left.
"""

from std.sys import argv
from std.time import perf_counter_ns

from karan import Aligner, Costs, Mode, Ties
from sporos import read_fasta

from plegma.anchors import Anchor
from plegma.backbone import (
    ALL,
    Level,
    MUMS,
    MUMS_OR_SMEMS,
    RARE,
    RARE_OR_MAXIMAL,
    REFERENCE_UNIQUE,
    SMEMS,
    backbone,
    windows,
)
from plegma.chain import ChainCosts
from plegma.graph import DELETION, INSERTION, cigar_runs


struct Strategy(Copyable, Movable):
    var text: String
    var levels: List[Level]
    var stay: Bool
    var sigma: Float64
    var trim: Int
    var budget: Int
    var lookback: Int

    def __init__(out self, text: String) raises:
        self.text = text
        self.levels = List[Level]()
        self.stay = False
        self.sigma = 0.1
        self.trim = 10
        self.budget = 0
        self.lookback = 200
        var parts = text.split(";")
        for k in range(1, len(parts)):
            var part = String(parts[k])
            if part == "stay":
                self.stay = True
            elif part.startswith("s="):
                self.sigma = Float64(String(part[byte=2:]))
            elif part.startswith("T="):
                self.trim = Int(String(part[byte=2:]))
            elif part.startswith("n="):
                self.budget = Int(String(part[byte=2:]))
            elif part.startswith("H="):
                self.lookback = Int(String(part[byte=2:]))
            else:
                raise Error(t"unknown strategy option {part} in {text}")
        var levels = String(parts[0])
        if levels == "none":
            return
        for item in levels.split(","):
            var token = String(item)
            if token.startswith("mum"):
                self.levels.append(Level(MUMS, Int(String(token[byte=3:]))))
            elif token.startswith("smem"):
                self.levels.append(Level(SMEMS, Int(String(token[byte=4:]))))
            elif token.startswith("ms"):
                self.levels.append(Level(MUMS_OR_SMEMS, Int(String(token[byte=2:]))))
            elif token.startswith("ref"):
                self.levels.append(Level(REFERENCE_UNIQUE, Int(String(token[byte=3:]))))
            elif token.startswith("rare"):
                var rest = String(token[byte=4:])
                var halves = rest.split("_")
                self.levels.append(Level(RARE, Int(String(halves[1])), Int(String(halves[0]))))
            elif token.startswith("combo"):
                var rest = String(token[byte=5:])
                var halves = rest.split("_")
                self.levels.append(Level(RARE_OR_MAXIMAL, Int(String(halves[1])), Int(String(halves[0]))))
            elif token.startswith("all"):
                self.levels.append(Level(ALL, Int(String(token[byte=3:]))))
            else:
                raise Error(t"unknown seed level {token} in {text}")


def _text(bytes: Span[UInt8, _]) -> String:
    return String(unsafe_from_utf8=List[UInt8](bytes))


def _partners(cigar: StringSlice, a_len: Int) raises -> List[Int]:
    """For each letter of a, the letter of b a CIGAR of b against a pairs it with (-1: none)."""
    var out = List[Int](length=a_len, fill=-1)
    var runs = cigar_runs(cigar)
    var i = 0
    var j = 0
    for r in range(len(runs.ops)):
        for _ in range(runs.counts[r]):
            var op = runs.ops[r]
            if op == DELETION:
                i += 1
            elif op == INSERTION:
                j += 1
            else:
                out[i] = j
                i += 1
                j += 1
    return out^


def _f1(ours: List[Int], truth: List[Int]) -> Float64:
    """F1 of the aligned letter pairs of `ours` against those of `truth`."""
    var both = 0
    var mine = 0
    var theirs = 0
    for i in range(len(ours)):
        mine += Int(ours[i] >= 0)
        theirs += Int(truth[i] >= 0)
        both += Int(ours[i] >= 0 and ours[i] == truth[i])
    return 2.0 * Float64(both) / Float64(mine + theirs) if mine + theirs > 0 else 1.0


def main() raises:
    var args = argv()
    if len(args) < 5:
        raise Error("usage: seeding <pairs.fa> <truth.tsv or -> <out.tsv> <strategy> ...")
    var records = read_fasta(String(args[1]))
    var truths = List[String]()
    if String(args[2]) != "-":
        with open(String(args[2]), "r") as f:
            # Split on newlines alone: this nightly's splitlines also splits at tabs.
            var content = f.read()
            var lines = content.split("\n")
            for k in range(1, len(lines)):
                var line = String(lines[k])
                if line.byte_length() == 0:
                    continue
                var fields = line.split("\t")
                truths.append(String(fields[len(fields) - 1]))
    var strategies = List[Strategy]()
    for k in range(4, len(args)):
        strategies.append(Strategy(String(args[k])))
    var costs = Costs.two_piece(4, 6, 2, 24, 1)
    var aligner = Aligner()
    var out = open(String(args[3]), "w")
    out.write("pair\tstrategy\ta_len\tb_len\toptimum\tcost\tf1_truth\tf1_optimum\toptimum_f1_truth\tanchors")
    out.write("\tanchors_true\tanchored\tlargest_window\tms\n")
    for k in range(len(records) // 2):
        var a = List[UInt8](records.seq(2 * k))
        var b = List[UInt8](records.seq(2 * k + 1))
        var name = records.names[2 * k]
        var a_text = _text(Span(a))
        var b_text = _text(Span(b))
        var best = aligner.align(a_text, b_text, costs, Mode.GLOBAL, ties=Ties.LEFT)
        var optimal = _partners(best.cigar, len(a))
        var has_truth = k < len(truths)
        var truth = _partners(truths[k], len(a)) if has_truth else List[Int]()
        var optimum_f1 = _f1(optimal, truth) if has_truth else -1.0
        for strategy in strategies:
            var t0 = perf_counter_ns()
            var chained = (
                backbone(
                    Span(a),
                    Span(b),
                    strategy.levels,
                    ChainCosts(substitution=strategy.sigma),
                    strategy.budget,
                    strategy.stay,
                    lookback=strategy.lookback,
                ) if len(strategy.levels)
                > 0 else List[Anchor]()
            )
            var anchors = List[Anchor]()
            for x in chained:
                if x.length > 2 * strategy.trim:
                    anchors.append(
                        Anchor(
                            x.ref_pos + strategy.trim,
                            x.query_pos + strategy.trim,
                            x.length - 2 * strategy.trim,
                            x.ref_copies,
                            x.query_copies,
                        )
                    )
            var ours = List[Int](length=len(a), fill=-1)
            var cost = 0
            var largest = 0
            var gaps = windows(anchors, len(a), len(b))
            for w in range(len(gaps)):
                var gap = gaps[w]
                largest = max(largest, gap.size())
                if gap.ref_end > gap.ref_start or gap.query_end > gap.query_start:
                    var found = aligner.align(
                        String(a_text[byte = gap.ref_start : gap.ref_end]),
                        String(b_text[byte = gap.query_start : gap.query_end]),
                        costs,
                        Mode.GLOBAL,
                        ties=Ties.LEFT,
                    )
                    cost += found.cost
                    var inside = _partners(found.cigar, gap.ref_end - gap.ref_start)
                    for i in range(len(inside)):
                        if inside[i] >= 0:
                            ours[gap.ref_start + i] = gap.query_start + inside[i]
                if w < len(anchors):
                    var x = anchors[w]
                    for t in range(x.length):
                        ours[x.ref_pos + t] = x.query_pos + t
            var ms = Float64(perf_counter_ns() - t0) / 1e6
            var true_anchors = 0
            var anchored = 0
            for x in anchors:
                anchored += x.length
                if has_truth:
                    var right = True
                    for t in range(x.length):
                        right = right and truth[x.ref_pos + t] == x.query_pos + t
                    true_anchors += Int(right)
            var truth_f1 = _f1(ours, truth) if has_truth else -1.0
            out.write(
                name,
                "\t",
                strategy.text,
                "\t",
                len(a),
                "\t",
                len(b),
                "\t",
                best.cost,
                "\t",
                cost,
                "\t",
                truth_f1,
                "\t",
                _f1(ours, optimal),
                "\t",
                optimum_f1,
                "\t",
                len(anchors),
                "\t",
                true_anchors,
                "\t",
                anchored,
                "\t",
                largest,
                "\t",
                ms,
                "\n",
            )
    out.close()
