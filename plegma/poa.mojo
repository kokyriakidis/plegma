"""Partial order alignment of a set of reads, as `plegma poa` runs it, reads anchored at both ends or at one.

A local assembly between two anchors has reads of three kinds: reads that span both anchors (`BOTH`), and reads
fixed at one anchor only, the left (`LEFT`, starting at the left anchor and ending inside the gap) or the right
(`RIGHT`). The spanning reads build the graph first (each pinned to the consensus so far and aligned between the
pins, `pinned.align_pinned`, or to the whole graph); then the left reads, each aligned from the graph's start to any
node (`kernel.PREFIX`: whatever passes a sink is new sequence the graph grows by), and the right reads, from any node
to the graph's end (`kernel.SUFFIX`). Their consensus runs from the spanning reads' start to their end, every read
weighing the edges (`Graph.heaviest_path_spanning`). A read fixed at one anchor only may be another copy's (a
tandem repeat's or the other haplotype's, sharing the one anchor and so missing the other): with spanning reads
enough, one much further from their consensus than any of them is (`_foreign`) gets no vote; it joins the graph after
the consensus is taken, and the caller leaves it out of polishing (`Poa.foreign`).

With no spanning read each side is assembled on its own: the left reads into a graph of their own, the longest
first, its consensus followed from the left anchor while most reads go on (`Graph.heaviest_path`); the right reads
likewise from the right anchor backwards (`Graph.heaviest_path_backward`). With both sides the two consensus
sequences are joined where they overlap, at the pin of minimap2's chain between them where both sides have the
most reads (or, without one, a short overlap karan finds between their ends); the graph is then rebuilt on the
joined sequence, every read aligned to it as anchored. With one side only (a telomere) its consensus is the
assembly, one-directional. Memory follows the adaptive band in every mode: the reads' length times the band, not
times each other.
"""

from std.math import isfinite

from karan import Aligner, Alignment, Costs, Mode

from .graph import Graph
from .kernel import (
    BOTH,
    EXTENSION,
    GLOBAL,
    LOCAL,
    LEFT,
    PREFIX,
    RIGHT,
    SINK,
    SOURCE,
    SUFFIX,
    GraphAlignment,
    Scoring,
    Topology,
    Workspace,
)
from .order import guide_order
from .pinned import WindowScratch, align_pinned, unique_pins
from .pins import pins
from .polish import past_the_end


comptime INPUT_ORDER = 0
comptime GUIDE_ORDER = 1
comptime LONGEST_FIRST = 2

comptime JOIN_ENDS = 2000
"""Without pins between the two sides' consensus sequences, their last and first this many letters are aligned for an
overlap."""
comptime JOIN_LEAST = 200
"""The fewest letters an overlap joining the two sides spans: shorter ones are as often two sides' tandem copies alike
as the same sequence."""
comptime JOIN_IDENTITY = 0.75
"""The least share of an overlap's letters a karan overlap matches."""
comptime JOIN_EDIT = 19
"""What an edit costs in that overlap, a match earning 1."""

comptime FOREIGN_SPANNING = 3
"""The fewest spanning reads whose consensus a read anchored at one end is checked against (`_foreign`)."""
comptime FOREIGN_RATIO = 1.5
comptime FOREIGN_SLACK = 3
"""A read anchored at one end is another copy's when its edit distance to the spanning reads' consensus passes
`FOREIGN_RATIO` times the worst spanning read's rate over the letters it covers, plus `FOREIGN_SLACK`."""


comptime BAND_MAX = 64
"""The adaptive band's reach at most, either side. abPOA's b + f L grows with the read (310 letters at 30 kb), but
the band also follows the best cells and the diagonal of the graph's remaining path, which carry a long read's
drift. Capped at 64, shasta2's long local assemblies take 2.6 times fewer cells and its chr21 assembly is byte for
byte the same (as at 32 and 100); the alignments in an MSA may place a gap elsewhere (237 of 36,082 rows over 2,093
local assemblies). Pinned windows (`align_pinned`) keep abPOA's band."""


struct PoaOptions(Copyable, Movable):
    """How `poa` aligns: the kernel's `mode` (`GLOBAL`, or abPOA's `LOCAL` and `EXTENSION` for spanning reads),
    its `scoring` and band, the spanning reads' `order`, whether each read is pinned to the consensus (`seeded`) and
    whether a read scoring poorly is tried reverse-complemented (`strand`, spanning reads only)."""

    var mode: Int
    var scoring: Scoring
    var band_b: Int
    var band_f: Float64
    var band_max: Int
    """The band's reach at most, either side (none when negative): `BAND_MAX`."""
    var order: Int
    var seeded: Bool
    var strand: Bool
    var min_cover: Int
    """A one-directional assembly (reads fixed at one anchor, none spanning) stops where fewer reads go on."""

    def __init__(out self):
        self.mode = GLOBAL
        self.scoring = Scoring()
        self.band_b = 10
        self.band_f = 0.01
        self.band_max = BAND_MAX
        self.order = GUIDE_ORDER
        self.seeded = True
        self.strand = False
        # Where five reads at 10% error agree a one-directional consensus errs about 0.2 times a kb; past it, with
        # three, about 1 (simulated telomeres, experiments/msa/anchored.py).
        self.min_cover = 5

    def check(self) raises:
        """Raises on options `poa` cannot run with, naming the value."""
        if self.mode != GLOBAL and self.mode != LOCAL and self.mode != EXTENSION:
            raise Error(t"mode {self.mode}: 0 (global), 1 (local) or 2 (extension)")
        var s = self.scoring
        if s.match_score < 0 or s.mismatch < 0 or s.open1 < 0 or s.open2 < 0 or s.extend1 <= 0 or s.extend2 < 0:
            raise Error("scoring: match, mismatch and gap opens must be 0 or more, the first gap extension more than 0")
        if not isfinite(self.band_f) or self.band_f < 0:
            raise Error(t"band_f {self.band_f}: must be 0 or more")
        if self.order != INPUT_ORDER and self.order != GUIDE_ORDER and self.order != LONGEST_FIRST:
            raise Error(t"order {self.order}: INPUT_ORDER, GUIDE_ORDER or LONGEST_FIRST")
        if self.min_cover < 1:
            raise Error(t"min_cover {self.min_cover}: must be 1 or more")


struct Poa(Movable):
    """What `poa` made: the graph; `shown`, the paths to output (a starting MSA's rows, then every read's, in input
    order); `spanning`, the reads (input indices) a consensus of all of them may be polished against; when any read
    is anchored at one end, `consensus` (one sequence, or the two sides' when they could not be joined), otherwise
    empty: the caller takes the graph's own; and `foreign`, the reads anchored at one end that are another copy's,
    left out of the consensus (and to be left out of its polishing)."""

    var graph: Graph
    var shown: List[Int]
    var spanning: List[Int]
    var consensus: List[List[UInt8]]
    var foreign: List[Int]

    def __init__(out self, var graph: Graph, var shown: List[Int], var spanning: List[Int]):
        self.graph = graph^
        self.shown = shown^
        self.spanning = spanning^
        self.consensus = List[List[UInt8]]()
        self.foreign = List[Int]()

    def parts(
        deinit self,
        mut graph: Graph,
        mut shown: List[Int],
        mut spanning: List[Int],
        mut consensus: List[List[UInt8]],
        mut foreign: List[Int],
    ):
        """The graph, `shown`, `spanning`, `consensus` and `foreign`, taken apart into those."""
        graph = self.graph^
        shown = self.shown^
        spanning = self.spanning^
        consensus = self.consensus^
        foreign = self.foreign^


struct _Aligner(Movable):
    """One read at a time into a graph: the scratch the kernel and the pinned windows reuse."""

    var row_of: List[Int]
    var scratch: WindowScratch
    var workspace: Workspace

    def __init__(out self):
        self.row_of = List[Int]()
        self.scratch = WindowScratch()
        self.workspace = Workspace()

    def align(
        mut self, graph: Graph, seq: List[UInt8], options: PoaOptions, mode: Int, anchors: List[Int] = List[Int]()
    ) raises -> GraphAlignment:
        """`seq` aligned to `graph` in `mode` (`options.mode` for a spanning read; `PREFIX` or `SUFFIX`). With
        `anchors` (paths), a read anchored at one end is aligned to the part of the graph between their ends alone
        (`_between`); pinned (`options.seeded`) only while that part is the whole graph, as the pins' windows reach
        the graph's every source and sink (a right read would end at a left read's letters past the right anchor)."""
        self.row_of.resize(len(graph), -1)
        if len(anchors) > 0 and (mode == PREFIX or mode == SUFFIX):
            var allowed = _between(graph, anchors, mode)
            if not options.seeded or False in allowed:
                var part = Topology(graph, allowed)
                return self.workspace.align(
                    graph,
                    part,
                    part.order,
                    SOURCE,
                    SINK,
                    Span(seq),
                    options.scoring,
                    self.row_of,
                    mode,
                    options.band_b,
                    options.band_f,
                    options.band_max,
                )
        var topology = Topology(graph)
        if options.seeded and (mode == GLOBAL or mode == PREFIX or mode == SUFFIX):
            return align_pinned(
                graph,
                topology,
                seq,
                options.scoring,
                options.band_b,
                options.band_f,
                self.row_of,
                self.scratch,
                self.workspace,
                mode,
            )
        return self.workspace.align(
            graph,
            topology,
            topology.order,
            SOURCE,
            SINK,
            Span(seq),
            options.scoring,
            self.row_of,
            mode,
            options.band_b,
            options.band_f,
            options.band_max,
        )


def poa(
    mut names: List[String],
    mut seqs: List[List[UInt8]],
    anchors: List[Int],
    options: PoaOptions,
    var start: Graph,
) raises -> Poa:
    """The POA of `seqs` (`names` theirs), each anchored as `anchors` says (`BOTH`, `LEFT` or `RIGHT`), into
    `start` (a graph made from an existing MSA, or empty). With `options.strand` a spanning read that scores poorly
    may be reverse-complemented (it and its name then change in `seqs` and `names`). No reads and an empty `start`:
    an empty graph and no consensus."""
    if len(seqs) == 0 and len(start.paths) == 0:
        return Poa(start^, List[Int](), List[Int]())
    var spanning = List[Int]()
    var left = List[Int]()
    var right = List[Int]()
    for r in range(len(seqs)):
        if anchors[r] == LEFT:
            left.append(r)
        elif anchors[r] == RIGHT:
            right.append(r)
        else:
            spanning.append(r)
    var shown = List[Int]()
    var starting = len(start.paths)
    for p in range(starting):
        shown.append(p)
    var path_of = List[Int](length=len(seqs), fill=-1)
    var aligner = _Aligner()
    var graph = start^
    if len(spanning) > 0 or len(graph.paths) > 0:
        _add_spanning(graph, names, seqs, _ordered(seqs, spanning, options.order), options, path_of, aligner)
        # The starting MSA's rows and the spanning reads end the consensus.
        var members = List[Int]()
        for p in range(starting):
            members.append(p)
        for r in spanning:
            members.append(path_of[r])
        var foreign = List[Int]()
        if len(spanning) >= FOREIGN_SPANNING and len(left) + len(right) > 0:
            var core = _letters(graph, graph.heaviest_path_spanning(Topology(graph).order, members))
            foreign = _foreign(core, seqs, spanning, left, right)
        var kept_left = _without(left, foreign)
        var kept_right = _without(right, foreign)
        _add_anchored(graph, names, seqs, _longest_first(seqs, kept_left), options, PREFIX, path_of, aligner, members)
        _add_anchored(graph, names, seqs, _longest_first(seqs, kept_right), options, SUFFIX, path_of, aligner, members)
        var consensus = List[List[UInt8]]()
        if len(left) + len(right) > 0:
            consensus.append(_letters(graph, graph.heaviest_path_spanning(Topology(graph).order, members)))
        # Another copy's reads in the graph too (the MSA shows every read), once the consensus is taken.
        _add_anchored(
            graph,
            names,
            seqs,
            _longest_first(seqs, _without(left, kept_left)),
            options,
            PREFIX,
            path_of,
            aligner,
            members,
        )
        _add_anchored(
            graph,
            names,
            seqs,
            _longest_first(seqs, _without(right, kept_right)),
            options,
            SUFFIX,
            path_of,
            aligner,
            members,
        )
        for r in range(len(seqs)):
            shown.append(path_of[r])
        var out = Poa(graph^, shown^, spanning^)
        out.consensus = consensus^
        out.foreign = foreign^
        return out^
    # No read spans both anchors: each side on its own, then joined.
    var parts = List[List[UInt8]]()
    if len(left) > 0:
        var side = _one_side(names, seqs, left, options, PREFIX)
        parts.append(_letters(side[0], side[0].heaviest_path(Topology(side[0]).order, side[1])))
    if len(right) > 0:
        var side = _one_side(names, seqs, right, options, SUFFIX)
        parts.append(_letters(side[0], side[0].heaviest_path_backward(Topology(side[0]).order, side[1])))
    if len(parts) == 2:
        var left_reach = List[Int](capacity=len(left))
        for r in left:
            left_reach.append(len(seqs[r]))
        var right_reach = List[Int](capacity=len(right))
        for r in right:
            right_reach.append(len(seqs[r]))
        var joined = _join(parts[0], parts[1], left_reach, right_reach)
        if len(joined) > 0:
            parts.clear()
            parts.append(joined^)
    elif len(left) > 0:
        _trim(parts[0], seqs, left, options.min_cover, True)
    else:
        _trim(parts[0], seqs, right, options.min_cover, False)
    # The graph every read is aligned into: on the joined sequence, or on the one side's.
    _ = graph.add("plegma_backbone", Span(parts[0]))
    if len(parts) == 1:
        var backbone: List[Int] = [0]
        _add_anchored(graph, names, seqs, _longest_first(seqs, left), options, PREFIX, path_of, aligner, backbone)
        _add_anchored(graph, names, seqs, _longest_first(seqs, right), options, SUFFIX, path_of, aligner, backbone)
        # The consensus again from the graph every read voted in, between the assembled sequence's ends (its path one
        # vote among the reads'): the reads mend what joining or one side's graph got wrong.
        parts[0] = _letters(
            graph, graph.heaviest_path_spanning(Topology(graph).order, backbone, len(left) > 0 and len(right) > 0)
        )
        # One side alone: the reads going on past the cut carry the walk past it again.
        if len(right) == 0:
            _trim(parts[0], seqs, left, options.min_cover, True)
        elif len(left) == 0:
            _trim(parts[0], seqs, right, options.min_cover, False)
    else:
        # Not joined: the left reads on the left side's sequence; the right reads beside it, on their own.
        _add_anchored(graph, names, seqs, _longest_first(seqs, left), options, PREFIX, path_of, aligner)
        var order = _longest_first(seqs, right)
        if len(order) > 0:
            path_of[order[0]] = graph.add(names[order[0]], Span(seqs[order[0]]))
            # Each aligned to the right side's part of the graph alone (what reaches the longest's end), not pinned
            # to the left side's path.
            var right_side: List[Int] = [path_of[order[0]]]
            for k in range(1, len(order)):
                var r = order[k]
                var found = aligner.align(graph, seqs[r], options, SUFFIX, right_side)
                path_of[r] = graph.fuse_nodes(names[r], Span(seqs[r]), found.nodes)
    for r in range(len(seqs)):
        shown.append(path_of[r])
    var out = Poa(graph^, shown^, spanning^)
    out.consensus = parts^
    return out^


def _foreign(
    core: List[UInt8], seqs: List[List[UInt8]], spanning: List[Int], left: List[Int], right: List[Int]
) raises -> List[Int]:
    """The reads of `left` and `right` (anchored at the left or the right end) another copy's: further from `core`,
    the spanning reads' consensus, than `FOREIGN_RATIO` times the worst spanning read's rate over the letters each
    covers, plus `FOREIGN_SLACK`. Edit distances; a one-sided read's as the polisher aligns it."""
    var aligner = Aligner(warm=False)
    var worst = 0.0
    for r in spanning:
        var found = aligner.align(Span(core), Span(seqs[r]), Costs.edit(), Mode.GLOBAL)
        worst = max(worst, Float64(found.cost) / Float64(max(1, len(core))))
    var out = List[Int]()
    for r in left:
        if _far(aligner, core, seqs[r], True, worst):
            out.append(r)
    for r in right:
        if _far(aligner, core, seqs[r], False, worst):
            out.append(r)
    sort(out)
    return out^


def _far(mut aligner: Aligner, core: List[UInt8], seq: List[UInt8], left: Bool, worst: Float64) raises -> Bool:
    """Whether `seq`, anchored at `core`'s start (`left`) or end, is further from it than `_foreign` allows."""
    var found = _one_sided(aligner, core, seq, left)
    var covered = found.reference_end - found.reference_start
    return Float64(found.cost) > FOREIGN_RATIO * worst * Float64(covered) + Float64(FOREIGN_SLACK)


def _one_sided(mut aligner: Aligner, consensus: List[UInt8], seq: List[UInt8], left: Bool) raises -> Alignment:
    """`seq`, anchored at the consensus's start (`left`) or end, against it: ending inside it (the whole read against
    a prefix, or a suffix) or past it (the whole consensus against a prefix of the read, or a suffix, the rest free),
    the cheaper."""
    var found = aligner.align(Span(consensus), Span(seq), Costs.edit(), Mode.PREFIX if left else Mode.SUFFIX)
    if len(consensus) - len(seq) < found.cost:
        var other = past_the_end(aligner, consensus, seq, left)
        if other.cost < found.cost:
            found = other^
    return found^


def _without(reads: List[Int], drop: List[Int]) -> List[Int]:
    """`reads` but those in `drop` (sorted)."""
    var out = List[Int]()
    for r in reads:
        var lo = 0
        var hi = len(drop)
        while lo < hi:
            var mid = (lo + hi) // 2
            if drop[mid] < r:
                lo = mid + 1
            else:
                hi = mid
        if lo == len(drop) or drop[lo] != r:
            out.append(r)
    return out^


def _add_spanning(
    mut graph: Graph,
    mut names: List[String],
    mut seqs: List[List[UInt8]],
    order: List[Int],
    options: PoaOptions,
    mut path_of: List[Int],
    mut aligner: _Aligner,
) raises:
    """The spanning reads in `order`, each aligned to the whole graph (or between pins) and fused into it; a read
    the same as one before it fused along that one's path (`_twin`)."""
    var seen = Dict[String, Int]()
    for r in order:
        if len(graph) == 0:
            path_of[r] = graph.add(names[r], Span(seqs[r]))
            if not options.strand:
                seen[String(unsafe_from_utf8=seqs[r])] = r
            continue
        if not options.strand and _twin(graph, names, seqs, r, seen, path_of):
            continue
        var found = aligner.align(graph, seqs[r], options, options.mode)
        if options.strand:
            # abPOA's rule: a score under a third of the best a match run of the shorter length could earn.
            var floor = Int(options.scoring.match_score) * min(len(seqs[r]), len(graph.paths[0])) // 3
            if found.score < floor:
                var flipped = _reverse_complement(seqs[r])
                var other = aligner.align(graph, flipped, options, options.mode)
                if other.score > found.score:
                    seqs[r] = flipped^
                    found = other^
                    names[r] += "_reverse_complement"
        path_of[r] = graph.fuse_nodes(names[r], Span(seqs[r]), found.nodes)
        if not options.strand:
            seen[String(unsafe_from_utf8=seqs[r])] = r


def _twin(
    mut graph: Graph,
    names: List[String],
    seqs: List[List[UInt8]],
    r: Int,
    seen: Dict[String, Int],
    mut path_of: List[Int],
) raises -> Bool:
    """Read `r` fused along the path of an earlier read with the same letters (in `seen`), if there is one: aligning
    it again would find that path, every letter matched (a local assembler's reads repeat often: on shasta2's chr21
    a third of them, the shorter ones)."""
    var key = String(unsafe_from_utf8=seqs[r])
    var twin = seen.get(key, -1)
    if twin < 0:
        return False
    var nodes = List[Int](capacity=len(seqs[r]))
    for node in graph.paths[path_of[twin]]:
        nodes.append(Int(node))
    path_of[r] = graph.fuse_nodes(names[r], Span(seqs[r]), nodes)
    return True


def _add_anchored(
    mut graph: Graph,
    names: List[String],
    seqs: List[List[UInt8]],
    order: List[Int],
    options: PoaOptions,
    mode: Int,
    mut path_of: List[Int],
    mut aligner: _Aligner,
    anchors: List[Int] = List[Int](),
) raises:
    """The reads in `order`, anchored at one end (`mode` `PREFIX` for the left, `SUFFIX` for the right), each
    aligned and fused into the graph; letters past the graph's end (or before its start) join it as new nodes.
    `anchors`: the paths whose ends are the anchors (the spanning reads, or the assembled backbone)."""
    var seen = Dict[String, Int]()
    for r in order:
        if _twin(graph, names, seqs, r, seen, path_of):
            continue
        var found = aligner.align(graph, seqs[r], options, mode, anchors)
        path_of[r] = graph.fuse_nodes(names[r], Span(seqs[r]), found.nodes)
        seen[String(unsafe_from_utf8=seqs[r])] = r


def _between(graph: Graph, anchors: List[Int], mode: Int) -> List[Bool]:
    """The nodes a read anchored at one end may align to: for a right read (`SUFFIX`) those an end of the paths
    `anchors` can still be reached from, for a left read (`PREFIX`) those reachable from a start of them. A read
    fixed at one anchor may run past the other (a local assembler cuts it to the gap's estimated length and a
    margin), and its letters past the graph's end join it as new nodes after the end: aligned to the whole graph,
    a right read would have to end at their end, the anchor no longer a sink (and a left read start at the new
    nodes before the start)."""
    var n = len(graph)
    var allowed = List[Bool](length=n, fill=False)
    var work = List[Int]()
    for p in anchors:
        ref path = graph.paths[p]
        if len(path) == 0:
            continue
        var node = Int(path[len(path) - 1]) if mode == SUFFIX else Int(path[0])
        if not allowed[node]:
            allowed[node] = True
            work.append(node)
    # Backwards from the ends (a right read), forwards from the starts (a left read).
    while len(work) > 0:
        var node = work.pop()
        var e = Int(graph.first_in[node]) if mode == SUFFIX else Int(graph.first_out[node])
        while e >= 0:
            var other = Int(graph.edge_from[e]) if mode == SUFFIX else Int(graph.edge_to[e])
            if not allowed[other]:
                allowed[other] = True
                work.append(other)
            e = Int(graph.next_in[e]) if mode == SUFFIX else Int(graph.next_out[e])
    return allowed^


def _one_side(
    names: List[String], seqs: List[List[UInt8]], reads: List[Int], options: PoaOptions, mode: Int
) raises -> Tuple[Graph, List[Int]]:
    """The graph of one side's reads alone: the longest first, the rest anchored as `mode` says; with its paths."""
    var graph = Graph()
    var order = _longest_first(seqs, reads)
    var paths = List[Int]()
    var aligner = _Aligner()
    for k in range(len(order)):
        var r = order[k]
        if k == 0:
            paths.append(graph.add(names[r], Span(seqs[r])))
            continue
        var found = aligner.align(graph, seqs[r], options, mode)
        paths.append(graph.fuse_nodes(names[r], Span(seqs[r]), found.nodes))
    return (graph^, paths^)


def _trim(mut consensus: List[UInt8], seqs: List[List[UInt8]], reads: List[Int], least: Int, from_left: Bool):
    """A one-directional consensus cut where fewer than `least` reads go on (or, with fewer reads than that, where the
    first of them stops): from the left anchor (`from_left`) past the `least`-th longest read's end, else before its
    start from the right anchor. The reads' lengths stand for their reach, as they all start (or end) at the anchor."""
    if least <= 1 or len(reads) == 0:
        return
    var lengths = List[Int]()
    for r in reads:
        lengths.append(len(seqs[r]))
    sort(lengths)
    var reach = lengths[len(lengths) - min(least, len(lengths))]
    if reach >= len(consensus):
        return
    if from_left:
        consensus.resize(reach, 0)
    else:
        var kept = List[UInt8](capacity=reach)
        for i in range(len(consensus) - reach, len(consensus)):
            kept.append(consensus[i])
        consensus = kept^


def _support(reach: List[Int], at: Int) -> Int:
    """How many reads reaching `reach` letters from their anchor get past `at` letters from it."""
    var n = 0
    for length in reach:
        if length > at:
            n += 1
    return n


def _join(left: List[UInt8], right: List[UInt8], left_reach: List[Int], right_reach: List[Int]) raises -> List[UInt8]:
    """`left` and `right` (the two sides' consensus sequences) joined where `right`'s start overlaps `left`'s end: at
    a pin of minimap2's chain between them whose 19-mer is unique in both (a tandem array's copies may chain to the
    wrong copy), when they place `right` starting inside `left` and reaching past its end by an overlap of
    `JOIN_LEAST` letters or more; else at an overlap karan finds between `left`'s last and `right`'s first `JOIN_ENDS`
    letters, as long and `JOIN_IDENTITY` alike; empty when they do not overlap.

    Of the pins, the one where the side worse supported there is best supported: left's letters before it come from
    the left reads that reach it (`left_reach`, each read's length from the left anchor), right's after it from the
    right reads that reach back to it (`right_reach`, from the right anchor). Each side's sequence near its far end is
    the few longest reads' alone. Taking the middle pin, HG00512's HLA-DRB1 haplotype 1 (no read spanning, 9 left
    reads and 18 right) had a stretch of 1.1 kb from the one right read reaching that far, against 9 left reads
    there, and its 9 insertions and 2 deletions (12 of 13 reads disagreeing at each) shifted exons 3 to 6 out of
    frame; the allele was read *08:01:05 for *08:03."""
    var found = pins(left, right)
    var unique = unique_pins(left, right, found)
    var kept = List[Int]()
    for m in range(len(found.r)):
        # Right starts inside left (q letters before the pin, r of left's) and ends past it.
        var r = found.r[m]
        var q = found.q[m]
        if unique[m] and r >= q and len(right) - q > len(left) - r and len(left) - (r - q) >= JOIN_LEAST:
            kept.append(m)
    if len(kept) > 0:
        # As well supported, the pin nearer the middle of the chain, where the join was taken before.
        var best = -1
        var best_support = -1
        var best_off = 0
        for k in range(len(kept)):
            var m = kept[k]
            var support = min(_support(left_reach, found.r[m]), _support(right_reach, len(right) - found.q[m]))
            var off = abs(k - len(kept) // 2)
            if support > best_support or (support == best_support and off < best_off):
                best = m
                best_support = support
                best_off = off
        var r = found.r[best]
        var q = found.q[best]
        var out = List[UInt8](capacity=r + len(right) - q)
        for i in range(r):
            out.append(left[i])
        for j in range(q, len(right)):
            out.append(right[j])
        return out^
    var a_start = max(0, len(left) - JOIN_ENDS)
    var a = List[UInt8](capacity=len(left) - a_start)
    for i in range(a_start, len(left)):
        a.append(left[i])
    var b = List[UInt8](capacity=min(len(right), JOIN_ENDS))
    for j in range(min(len(right), JOIN_ENDS)):
        b.append(right[j])
    var mode = Mode.ends_free(reference_start=len(a), query_end=len(b)).with_match_score(1)
    var aligner = Aligner(warm=False)
    # Edits far dearer than a match earns (minimap2's 19 to 1): a longer overlap shifted by tandem copies, its letters
    # a few percent apart, loses to the true one, its letters as alike as the two consensus sequences.
    var got = aligner.align(Span(a), Span(b), Costs.linear(JOIN_EDIT, JOIN_EDIT), mode)
    var spans = got.query_end - got.query_start
    if got.reference_end != len(a) or got.query_start != 0 or spans < JOIN_LEAST:
        return List[UInt8]()
    if Float64(got.cost) > Float64(JOIN_EDIT) * (1.0 - JOIN_IDENTITY) * Float64(spans):
        return List[UInt8]()
    var out = List[UInt8](capacity=len(left) + len(right) - got.query_end)
    for i in range(len(left)):
        out.append(left[i])
    for j in range(got.query_end, len(right)):
        out.append(right[j])
    return out^


def _ordered(seqs: List[List[UInt8]], reads: List[Int], order: Int) raises -> List[Int]:
    """`reads` in input order, abPOA's guide order or longest first."""
    if order == GUIDE_ORDER and len(reads) > 1:
        var subset = List[List[UInt8]]()
        for r in reads:
            subset.append(seqs[r].copy())
        var out = List[Int]()
        for k in guide_order(subset):
            out.append(reads[k])
        return out^
    if order == LONGEST_FIRST:
        return _longest_first(seqs, reads)
    return reads.copy()


def _longest_first(seqs: List[List[UInt8]], reads: List[Int]) -> List[Int]:
    """`reads` longest first, input order among equals."""
    var keys = List[UInt64]()
    for r in reads:
        keys.append((UInt64((1 << 32) - len(seqs[r])) << 32) | UInt64(r))
    sort(keys)
    var out = List[Int](capacity=len(keys))
    for key in keys:
        out.append(Int(key & 0xFFFFFFFF))
    return out^


def _letters(graph: Graph, path: List[Int]) -> List[UInt8]:
    var out = List[UInt8](capacity=len(path))
    for node in path:
        out.append(graph.bases[node])
    return out^


def _reverse_complement(seq: List[UInt8]) -> List[UInt8]:
    var out = List[UInt8](capacity=len(seq))
    for i in reversed(range(len(seq))):
        var c = seq[i]
        var lower = c | 0x20
        if lower == UInt8(ord("a")):
            out.append(UInt8(ord("T")))
        elif lower == UInt8(ord("c")):
            out.append(UInt8(ord("G")))
        elif lower == UInt8(ord("g")):
            out.append(UInt8(ord("C")))
        elif lower == UInt8(ord("t")) or lower == UInt8(ord("u")):
            out.append(UInt8(ord("A")))
        else:
            out.append(c)
    return out^
