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

from .alphabet import encode_base
from .graph import DELETION, INSERTION, MATCH, Graph, cigar_runs
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
from .bands import past_the_end
from .trace import POA, Trace


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
comptime JOIN_EDIT_STRICT = 100
"""The same, for the overlap of the ends found again with edits dearer still (a longer overlap shifted by tandem
copies a few percent apart loses to the true one even then)."""
comptime JOIN_MAX_PERIOD = 500
comptime JOIN_PERIOD_K = 8
comptime JOIN_PERIOD_SHARE = 0.3
"""An overlap whose `JOIN_PERIOD_K`-mers recur one period on (at most `JOIN_MAX_PERIOD`) this often is a tandem
array's (`_period`)."""
comptime JOIN_SHIFTS = 64
"""The most copies either way a placement inside a tandem array is shifted by, to be scored too."""
comptime JOIN_DRIFT = 16
"""Pins whose diagonals are within this many letters of each other place right alike (the two sides' insertions and
deletions move it)."""
comptime JOIN_MARGIN = 1.5
comptime JOIN_MARGIN_FLOOR = 0.002
"""A placement on another diagonal agreeing within `JOIN_MARGIN` times the best's edits a letter, plus
`JOIN_MARGIN_FLOOR`, leaves the two sides unjoined: a tandem array whose copies are as alike as the two consensus
sequences."""

comptime FOREIGN_SPANNING = 3
"""The fewest spanning reads whose consensus a read anchored at one end is checked against (`_foreign`)."""
comptime FOREIGN_RATIO = 1.5
comptime FOREIGN_SLACK = 3
"""A read anchored at one end is another copy's when its edit distance to the spanning reads' consensus passes
`FOREIGN_RATIO` times the worst spanning read's rate over the letters it covers, plus `FOREIGN_SLACK`."""


comptime ANCHOR_SPREAD = 16
"""A one-sided assembly's anchored end is voted on over this many nodes (`Graph._spread_vote`): HG002's HLA-DRB1
under GRCh38, 12 left reads starting over 14 nodes (no more than 2 at one) and 8 short reads at one other place,
whose start the side took, its first 250 letters not the haplotype's."""
comptime END_LETTERS = 100
"""A spanning read's end checked against the reads anchored at that end (`_ends_checked`): its first or last this
many letters."""
comptime END_EDITS = 35
"""Two reads' ends agree within this many edits over `END_LETTERS`: two reads at 10% errors differ by some 20 there,
two unrelated sequences by some 52."""
comptime _STRAY = -1
"""A read given as spanning whose ends both disagree with the reads anchored there (`_ends_checked`)."""
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
    left out of the consensus (and to be left out of its polishing), and the reads given as spanning whose ends agree
    with neither anchor's reads; `reversed`, the reads `options.strand`
    reverse-complemented (on the other strand to the consensus now), ascending."""

    var graph: Graph
    var shown: List[Int]
    var spanning: List[Int]
    var consensus: List[List[UInt8]]
    var foreign: List[Int]
    var reversed: List[Int]
    var anchors: List[Int]
    """Each read's anchor as built with (`_ends_checked`), when any read is anchored at one end; otherwise empty."""

    def __init__(out self, var graph: Graph, var shown: List[Int], var spanning: List[Int]):
        self.graph = graph^
        self.shown = shown^
        self.spanning = spanning^
        self.consensus = List[List[UInt8]]()
        self.foreign = List[Int]()
        self.reversed = List[Int]()
        self.anchors = List[Int]()

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
        (`_between`). Pinned (`options.seeded`), the pins' windows reach the graph's every source and sink, so a
        right read could end on a left read's letters past the right anchor (or a left read's dead-end tip), or on
        another right read's letters past the right anchor itself (a read anchored by an estimate of where the
        anchor is, not by its alignment): its letters on nodes outside the part are the read's letters past its
        anchor, and are inserted there (`_past_anchor`). Unpinned, the read is aligned to the part alone."""
        self.row_of.resize(len(graph), -1)
        if len(anchors) > 0 and (mode == PREFIX or mode == SUFFIX):
            var allowed = _between(graph, anchors, mode)
            if options.seeded:
                var topology = Topology(graph)
                var found = align_pinned(
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
                _past_anchor(found, allowed)
                return found^
            else:
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
    var anchor_of = _ends_checked(seqs, anchors)
    var spanning = List[Int]()
    var left = List[Int]()
    var right = List[Int]()
    var strays = List[Int]()
    for r in range(len(seqs)):
        if anchor_of[r] == LEFT:
            left.append(r)
        elif anchor_of[r] == RIGHT:
            right.append(r)
        elif anchor_of[r] == _STRAY:
            strays.append(r)
        else:
            spanning.append(r)
    var trace = Trace()
    if trace.on(POA):
        for r in range(len(seqs)):
            if anchor_of[r] != anchors[r]:
                var now = "left" if anchor_of[r] == LEFT else ("right" if anchor_of[r] == RIGHT else "neither")
                trace.line(POA, String("read ", r, " given as spanning: its ends agree with the ", now, " reads'"))
    if trace.on(POA):
        var order_name = "guide" if options.order == GUIDE_ORDER else (
            "longest first" if options.order == LONGEST_FIRST else "input"
        )
        trace.line(
            POA,
            String(
                len(seqs),
                " reads: ",
                len(spanning),
                " spanning, ",
                len(left),
                " left, ",
                len(right),
                " right; ",
                len(start.paths),
                " starting rows; ",
                order_name,
                " order, ",
                "pinned" if options.seeded else "unpinned",
            ),
        )
    var shown = List[Int]()
    var starting = len(start.paths)
    for p in range(starting):
        shown.append(p)
    var path_of = List[Int](length=len(seqs), fill=-1)
    var aligner = _Aligner()
    var graph = start^
    if len(spanning) > 0 or len(graph.paths) > 0:
        var reversed = _add_spanning(
            graph, names, seqs, _ordered(seqs, spanning, options.order), options, path_of, aligner
        )
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
            if trace.on(POA) and len(foreign) > 0:
                trace.line(POA, String("another copy's (no vote): reads ", _listed(foreign)))
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
        _add_strays(graph, names, seqs, strays, options, path_of, aligner)
        for r in range(len(seqs)):
            shown.append(path_of[r])
        var out = Poa(graph^, shown^, spanning^)
        out.consensus = consensus^
        foreign.extend(strays.copy())
        sort(foreign)
        out.foreign = foreign^
        out.anchors = _as_built(anchor_of)
        sort(reversed)
        if trace.on(POA) and len(reversed) > 0:
            trace.line(POA, String("reverse-complemented: reads ", _listed(reversed)))
        out.reversed = reversed^
        return out^
    # No read spans both anchors: each side on its own, then joined.
    var parts = List[List[UInt8]]()
    # Each side's reads' reach along its consensus, from its anchor (`_reach`).
    var left_reach = List[Int]()
    var right_reach = List[Int]()
    if len(left) > 0:
        var side = _one_side(names, seqs, left, options, PREFIX)
        var path = side[0].heaviest_path(Topology(side[0]).order, side[1], ANCHOR_SPREAD)
        parts.append(_letters(side[0], path))
        left_reach = _reach(side[0], path, side[1], True)
    if len(right) > 0:
        var side = _one_side(names, seqs, right, options, SUFFIX)
        var path = side[0].heaviest_path_backward(Topology(side[0]).order, side[1], ANCHOR_SPREAD)
        parts.append(_letters(side[0], path))
        right_reach = _reach(side[0], path, side[1], False)
    if len(parts) == 2:
        var joined = _join(parts[0], parts[1], left_reach, right_reach)
        trace.line(
            POA,
            String(
                "no spanning read: left side ",
                len(parts[0]),
                " letters, right ",
                len(parts[1]),
                "; ",
                String("joined, ", len(joined), " letters") if len(joined) > 0 else String("not joined"),
            ),
        )
        if len(joined) > 0:
            parts.clear()
            parts.append(joined^)
    elif len(left) > 0:
        _trim(parts[0], left_reach, options.min_cover, True)
        trace.line(
            POA,
            String("left side alone: ", len(parts[0]), " letters, cut where under ", options.min_cover, " reads go on"),
        )
    else:
        _trim(parts[0], right_reach, options.min_cover, False)
        trace.line(
            POA,
            String(
                "right side alone: ", len(parts[0]), " letters, cut where under ", options.min_cover, " reads go on"
            ),
        )
    # The graph every read is aligned into: on the joined sequence, or on the one side's.
    _ = graph.add("plegma_backbone", Span(parts[0]))
    if len(parts) == 1:
        var backbone: List[Int] = [0]
        _add_anchored(graph, names, seqs, _longest_first(seqs, left), options, PREFIX, path_of, aligner, backbone)
        _add_anchored(graph, names, seqs, _longest_first(seqs, right), options, SUFFIX, path_of, aligner, backbone)
        # The consensus again from the graph every read voted in, between the assembled sequence's ends (its path one
        # vote among the reads'): the reads mend what joining or one side's graph got wrong.
        var path = graph.heaviest_path_spanning(Topology(graph).order, backbone, len(left) > 0 and len(right) > 0)
        parts[0] = _letters(graph, path)
        # One side alone: the reads going on past the cut carry the walk past it again.
        if len(right) == 0 or len(left) == 0:
            var paths = List[Int]()
            for r in left if len(right) == 0 else right:
                paths.append(path_of[r])
            _trim(parts[0], _reach(graph, path, paths, len(right) == 0), options.min_cover, len(right) == 0)
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
    _add_strays(graph, names, seqs, strays, options, path_of, aligner)
    for r in range(len(seqs)):
        shown.append(path_of[r])
    var out = Poa(graph^, shown^, spanning^)
    out.consensus = parts^
    out.foreign = strays^
    out.anchors = _as_built(anchor_of)
    return out^


def _ends_checked(seqs: List[List[UInt8]], anchors: List[Int]) raises -> List[Int]:
    """`anchors`, each spanning read's ends checked against the reads anchored at that end: where most of those (the
    one-sided reads anchored there and the other spanning reads, at least `END_LETTERS` long) disagree with its end
    (`_ends_agree`), it is anchored at its other end alone (`LEFT` or `RIGHT`), or at neither (`_STRAY`: another
    copy's, or another place's). Only at an end that one-sided reads are anchored at, so reads all spanning are as
    given. A read anchored by counting letters to an anchor ends where the haplotype's length is not the
    reference's: HG03009's HLA-DQA1 haplotype 1 had one read given as spanning, ending 5.3 kb before the end its 22
    right reads reach, and the consensus ended there, 4.3 kb of the right reads' letters inserted before its end."""
    var out = anchors.copy()
    var left = List[Int]()
    var right = List[Int]()
    var spanning = List[Int]()
    for r in range(len(seqs)):
        if anchors[r] == LEFT:
            left.append(r)
        elif anchors[r] == RIGHT:
            right.append(r)
        elif len(seqs[r]) >= END_LETTERS:
            spanning.append(r)
    if len(spanning) == 0 or len(left) + len(right) == 0:
        return out^
    var aligner = Aligner(warm=False)
    for s in spanning:
        var held = 0
        for at_start in [True, False]:
            var ones = left.copy() if at_start else right.copy()
            if len(ones) == 0:
                held += 1 if at_start else 2
                continue
            var agree = 0
            var disagree = 0
            for r in ones + spanning.copy():
                if r == s or len(seqs[r]) < END_LETTERS:
                    continue
                if _ends_agree(aligner, seqs[s], seqs[r], at_start):
                    agree += 1
                else:
                    disagree += 1
            if disagree <= agree:
                held += 1 if at_start else 2
        # held: 1 the start, 2 the end, 3 both (an end with no one-sided read holds).
        if held == 0:
            out[s] = _STRAY
        elif held == 1:
            out[s] = LEFT
        elif held == 2:
            out[s] = RIGHT
    return out^


def _ends_agree(mut aligner: Aligner, a: List[UInt8], b: List[UInt8], at_start: Bool) raises -> Bool:
    """Whether `a`'s and `b`'s first (`at_start`) or last `END_LETTERS` letters are within `END_EDITS` edits."""
    var x = Span(a)[0:END_LETTERS] if at_start else Span(a)[len(a) - END_LETTERS :]
    var y = Span(b)[0:END_LETTERS] if at_start else Span(b)[len(b) - END_LETTERS :]
    return aligner.align(x, y, Costs.edit(), Mode.GLOBAL).cost <= END_EDITS


def _as_built(anchor_of: List[Int]) -> List[Int]:
    """Each read's anchor as `poa` built with it, a stray as `BOTH` (it is in `Poa.foreign`, and polishes nothing)."""
    var out = anchor_of.copy()
    for r in range(len(out)):
        if out[r] == _STRAY:
            out[r] = BOTH
    return out^


def _add_strays(
    mut graph: Graph,
    names: List[String],
    seqs: List[List[UInt8]],
    strays: List[Int],
    options: PoaOptions,
    mut path_of: List[Int],
    mut aligner: _Aligner,
) raises:
    """The reads `_ends_checked` found anchored at neither end, into the graph once the consensus is taken (the MSA
    shows every read), aligned as spanning reads are."""
    for r in strays:
        var found = aligner.align(graph, seqs[r], options, options.mode)
        path_of[r] = graph.fuse_nodes(names[r], Span(seqs[r]), found.nodes)


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
) raises -> List[Int]:
    """The spanning reads in `order`, each aligned to the whole graph (or between pins) and fused into it; a read
    the same as one before it fused along that one's path (`_twin`). The reads reverse-complemented (`strand`)."""
    var seen = Dict[String, Int]()
    var reversed = List[Int]()
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
                    reversed.append(r)
        path_of[r] = graph.fuse_nodes(names[r], Span(seqs[r]), found.nodes)
        if not options.strand:
            seen[String(unsafe_from_utf8=seqs[r])] = r
    return reversed^


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


def _past_anchor(mut found: GraphAlignment, allowed: List[Bool]):
    """A one-sided read's letters `found` aligns to nodes outside `allowed` (`_between`) made inserted letters: those
    past the read's own anchor. They are one run at the read's anchored end, as no node outside the part leads into
    it (a right read's after the anchor, a left read's before it), and the letters before them keep their pinned
    places. The part's alignment of such a read alone, unpinned, lost the read's path in the band whenever its
    letters past the anchor outran the band (kilobases where a read was anchored by counting letters to the anchor):
    on HG03009's HLA-DQA1 every right read after the first that ran past it, and the consensus 4.6 kb short."""
    for i in range(len(found.nodes)):
        var node = found.nodes[i]
        if node >= 0 and not allowed[node]:
            found.nodes[i] = -1


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


def _reach(graph: Graph, consensus: List[Int], paths: List[Int], from_left: Bool) -> List[Int]:
    """How far along `consensus` (its nodes) each of `paths` reaches from the anchor it starts at (`from_left`: from
    the consensus's start, else back from its end): letters to the furthest consensus column it passes (a letter
    unlike the consensus's there is in its column; an inserted one in none). A read's own length would stand for it
    only where the consensus has the read's letters one for one."""
    var at_column = List[Int](length=graph.column_count(), fill=-1)
    for k in range(len(consensus)):
        at_column[Int(graph.columns[consensus[k]])] = k
    var out = List[Int](capacity=len(paths))
    for p in paths:
        var furthest = 0
        for node in graph.paths[p]:
            var k = at_column[Int(graph.columns[Int(node)])]
            if k >= 0:
                furthest = max(furthest, k + 1 if from_left else len(consensus) - k)
        out.append(furthest)
    return out^


def _trim(mut consensus: List[UInt8], reach: List[Int], least: Int, from_left: Bool):
    """A one-directional consensus cut where fewer than `least` reads go on (or, with fewer reads than that, where the
    first of them stops): from the left anchor (`from_left`) past the `least`-th furthest read's reach along it
    (`_reach`), else before it from the right anchor. With `least` 1 the cut is past the furthest read."""
    if least <= 0 or len(reach) == 0:
        return
    var reaches = reach.copy()
    sort(reaches)
    var cut = reaches[len(reaches) - min(least, len(reaches))]
    if cut >= len(consensus):
        return
    if from_left:
        consensus.resize(cut, 0)
    else:
        var kept = List[UInt8](capacity=cut)
        for i in range(len(consensus) - cut, len(consensus)):
            kept.append(consensus[i])
        consensus = kept^


def _support(reach: List[Int], at: Int) -> Int:
    """How many reads reaching `reach` letters from their anchor get past `at` letters from it."""
    var n = 0
    for length in reach:
        if length > at:
            n += 1
    return n


@fieldwise_init
struct _Placing(Copyable, Movable):
    """A way to join the two sides: right's letters from `cut_right` on after left's first `cut_left` (right
    starting at `diagonal` in left's letters there); how well the two agree over the overlap it implies (`edits` in
    `letters`), and where it puts each of left's letters in right (`path`)."""

    var diagonal: Int
    var cut_left: Int
    var cut_right: Int
    var edits: Int
    var letters: Int
    var how: String
    var path: List[Int]
    """The right letter each of left's overlap letters lands on, from left's letter `path_from`."""
    var path_from: Int
    var edits_after: List[Int]
    """The alignment's edits at left's letters `path_from + i` and on, by `i` (one
    more, 0, past the last), so two are compared on the left letters both align."""

    def edits_from(self, at_left: Int) -> Int:
        """This placement's edits at left's letters from `at_left` on."""
        var i = at_left - self.path_from
        if i <= 0:
            return self.edits_after[0] if len(self.edits_after) > 0 else self.edits
        return self.edits_after[min(i, len(self.edits_after) - 1)]

    def lands(self, at_left: Int) -> Int:
        """The right letter left's letter `at_left` lands on (outside the overlap: as the diagonal puts it)."""
        var k = at_left - self.path_from
        if k >= 0 and k < len(self.path):
            return self.path[k]
        return at_left - self.diagonal

    def rate(self) -> Float64:
        return Float64(self.edits) / Float64(max(1, self.letters))


def _walk(found: Alignment, left_from: Int, right_from: Int, mut path: List[Int]) raises:
    """Append to `path` the right letter each left letter of `found` lands on (left the query, starting at
    `left_from`; right the reference, at `right_from` plus the alignment's start in it; a left letter right lacks at
    the next right letter)."""
    var runs = cigar_runs(found.cigar)
    var j = right_from + found.reference_start
    for x in range(len(runs.ops)):
        var op = runs.ops[x]
        for _ in range(runs.counts[x]):
            if op == DELETION:
                j += 1
            else:
                path.append(j)
                if op != INSERTION:
                    j += 1


def _placing(
    mut aligner: Aligner, left: List[UInt8], right: List[UInt8], cut_left: Int, cut_right: Int, var how: String
) raises -> _Placing:
    """Right placed so that left's letter `cut_left` meets right's `cut_right` (a pin; or left's end meeting right's
    `cut_right`, an overlap of the ends): the two sides aligned from there back to the overlap's start and on to
    left's end, right's far ends free (the two sides' insertions and deletions move the overlap's ends); its edits,
    its letters in left, and the right letter each of them lands on (`path`, from the overlap's start)."""
    var lo = max(0, cut_left - cut_right)
    var out = _Placing(cut_left - cut_right, cut_left, cut_right, 0, len(left) - lo, how^, List[Int](), lo, List[Int]())
    # An alignment past karan's memory (a stray pin's, between unrelated stretches) raises, and the caller leaves that
    # placement out.
    var at = List[Int](length=len(left) - lo + 1, fill=0)
    if cut_left > lo:
        # Back from the meeting point: left's letters before it against a suffix of right's before it.
        var back = aligner.align(Span(right)[0:cut_right], Span(left)[lo:cut_left], Costs.edit(), Mode.SUFFIX)
        out.edits += back.cost
        _walk(back, lo, 0, out.path)
        _tally(back, 0, at)
    if cut_left < len(left):
        # On from it: left's letters after it against a prefix of right's after it.
        var on = aligner.align(Span(right)[cut_right:], Span(left)[cut_left:], Costs.edit(), Mode.PREFIX)
        out.edits += on.cost
        _walk(on, cut_left, cut_right, out.path)
        _tally(on, cut_left - lo, at)
    for k in reversed(range(len(at) - 1)):
        at[k] += at[k + 1]
    out.edits_after = at^
    return out^


def _tally(found: Alignment, first: Int, mut at: List[Int]) raises:
    """Each edit of `found` at the left letter (its query, from `at[first]`) it is at; a right letter left lacks, at
    the next left letter."""
    var runs = cigar_runs(found.cigar)
    var i = first
    for x in range(len(runs.ops)):
        var op = runs.ops[x]
        for _ in range(runs.counts[x]):
            if op == DELETION:
                at[min(i, len(at) - 1)] += 1
            else:
                if op != MATCH:
                    at[i] += 1
                i += 1


def _copy_join(
    left: List[UInt8], right: List[UInt8], hypotheses: List[_Placing], period: Int, trace: Trace
) -> List[UInt8]:
    """The two sides joined at the copy number whose placement has the fewest edits a letter, when every other copy
    number has more than `JOIN_MARGIN` times as many plus `JOIN_MARGIN_FLOOR`, on the left letters both align
    (`_versus`; a longer overlap's thin far end does not count against it); else apart: the copies there are too
    alike to count."""
    var best = _best(hypotheses)
    if best < 0:
        trace.line(POA, String("not joined: no copy number ", JOIN_IDENTITY, " alike over ", JOIN_LEAST, " letters"))
        return List[UInt8]()
    ref chosen = hypotheses[best]
    for c in range(len(hypotheses)):
        ref h = hypotheses[c]
        if c == best or h.letters < JOIN_LEAST or abs(h.diagonal - chosen.diagonal) < period // 2:
            continue
        var start = max(h.path_from, chosen.path_from)
        var e = _versus(h, chosen)
        var letters = Float64(max(1, len(left) - start))
        var mine = Float64(e[1]) / letters
        var theirs = Float64(e[0]) / letters
        if theirs <= mine * JOIN_MARGIN + JOIN_MARGIN_FLOOR:
            trace.line(
                POA,
                String(
                    "not joined: ",
                    chosen.how,
                    " has ",
                    e[1],
                    " edits, ",
                    h.how,
                    " ",
                    e[0],
                    " in the ",
                    len(left) - start,
                    " left letters both align",
                ),
            )
            return List[UInt8]()
    trace.line(POA, String("join at ", chosen.how, ": ", chosen.edits, " edits in ", chosen.letters, " letters"))
    var out = List[UInt8](capacity=chosen.cut_left + len(right) - chosen.cut_right)
    for i in range(chosen.cut_left):
        out.append(left[i])
    for j in range(chosen.cut_right, len(right)):
        out.append(right[j])
    return out^


def _best(placings: List[_Placing]) -> Int:
    """The placement with the fewest edits a letter (of two as good, the longer overlap), among those over `JOIN_LEAST`
    letters and `JOIN_IDENTITY` alike; -1 for none."""
    var best = -1
    for c in range(len(placings)):
        ref p = placings[c]
        if p.letters < JOIN_LEAST or p.rate() > 1.0 - JOIN_IDENTITY:
            continue
        if (
            best < 0
            or p.rate() < placings[best].rate()
            or (p.rate() == placings[best].rate() and p.letters > placings[best].letters)
        ):
            best = c
    return best


def _period(seq: List[UInt8], start: Int, end: Int) -> Int:
    """The period `seq[start:end]` repeats with, or 0 for none: the distance at which most of its `JOIN_PERIOD_K`-mers
    recur (each to its last occurrence, at most `JOIN_MAX_PERIOD` back), counted within a twentieth of it either side
    so copies differing by insertions and deletions (a noisy consensus's, a VNTR's) still count; at least
    `JOIN_PERIOD_SHARE` of the k-mers, and two copies. (Comparing letters a fixed period apart lost the phase at the
    first indel.)"""
    var n = end - start
    if n < 2 * JOIN_PERIOD_K:
        return 0
    var at_distance = List[Int](length=JOIN_MAX_PERIOD + 1, fill=0)
    var last = Dict[UInt64, Int]()
    var kmers = 0
    var key = UInt64(0)
    var valid = 0
    var mask = (UInt64(1) << UInt64(2 * JOIN_PERIOD_K)) - 1
    for x in range(start, end):
        var code = encode_base(seq[x])
        valid = valid + 1 if code < 4 else 0
        key = ((key << 2) | UInt64(code & 3)) & mask
        if valid < JOIN_PERIOD_K:
            continue
        kmers += 1
        var seen = last.get(key)
        if seen and x - seen.value() <= JOIN_MAX_PERIOD:
            at_distance[x - seen.value()] += 1
        last[key] = x
    var best = 0
    var best_count = 0
    for p in range(2, min(JOIN_MAX_PERIOD, n // 2) + 1):
        var width = max(1, p // 20)
        var count = 0
        for d in range(max(2, p - width), min(JOIN_MAX_PERIOD, p + width) + 1):
            count += at_distance[d]
        if count > best_count:
            best = p
            best_count = count
    if Float64(best_count) < JOIN_PERIOD_SHARE * Float64(kmers):
        return 0
    return best


def _same_phase(right: List[UInt8], at: Int, near: Int, half: Int) -> Int:
    """The letter of right within `half` of `near` whose surroundings (`JOIN_PERIOD_K` either side) are most like
    those of right's letter `at`: the same place in another copy of a tandem array, whatever indels make its copies'
    lengths differ."""
    var best = near
    var best_same = -1
    for j in range(max(0, near - half), min(len(right), near + half + 1)):
        var same = 0
        for t in range(-JOIN_PERIOD_K, JOIN_PERIOD_K):
            var a = at + t
            var b = j + t
            if a >= 0 and b >= 0 and a < len(right) and b < len(right) and right[a] == right[b]:
                same += 1
        if same > best_same or (same == best_same and abs(j - near) < abs(best - near)):
            best = j
            best_same = same
    return best


def _versus(a: _Placing, b: _Placing) -> Tuple[Int, Int]:
    """Two copy hypotheses' edits on the left letters both align (from the later start of the two): a longer
    overlap's noisiest letters, the far ends', do not count against it."""
    var start = max(a.path_from, b.path_from)
    return (a.edits_from(start), b.edits_from(start))


def _join(left: List[UInt8], right: List[UInt8], left_reach: List[Int], right_reach: List[Int]) raises -> List[UInt8]:
    """`left` and `right` (the two sides' consensus sequences) joined where `right`'s start overlaps `left`'s end, by
    `JOIN_LEAST` letters or more; empty when they do not overlap, or when where they do is not clear.

    The ways to place right: the pins of minimap2's chain between them whose 19-mer is unique in both, grouped by
    the diagonal they put right on (`JOIN_DRIFT`), each placing right so that it reaches left's end; and the
    overlaps karan finds between left's last and right's first `JOIN_ENDS` letters, a match earning 1 against edits
    of `JOIN_EDIT` and of `JOIN_EDIT_STRICT`. Each is scored by how well the two sides agree over all the overlap it
    implies: the fewest edits a letter wins, if under `1 - JOIN_IDENTITY`. One 19-mer cannot place a join inside a
    tandem array, where a copy's 19-mer may be unique in both sides and still another copy's; the whole overlap's
    agreement can, as the two sides' copies differ. Where a placement on another diagonal agrees nearly as well
    (`JOIN_MARGIN`), the copies there are too alike to tell, and the two sides are not joined.

    At a pin's diagonal, the join is at the pin where the side worse supported there is best supported: left's
    letters before it come from the left reads that reach it (`left_reach`, each read's reach along left from the
    left anchor), right's after it from the right reads that reach back to it (`right_reach`). Each side's sequence
    near its far end is the few longest reads' alone. Taking the middle pin, HG00512's HLA-DRB1 haplotype 1 (no read
    spanning, 9 left reads and 18 right) had a stretch of 1.1 kb from the one right read reaching that far, against 9
    left reads there, and its 9 insertions and 2 deletions (12 of 13 reads disagreeing at each) shifted exons 3 to 6
    out of frame; the allele was read *08:01:05 for *08:03."""
    var trace = Trace()
    var aligner = Aligner(warm=False)
    var placings = List[_Placing]()
    var found = pins(left, right)
    var unique = unique_pins(left, right, found)
    # The pins that place right reaching left's end with JOIN_LEAST letters of overlap, by diagonal.
    var keys = List[UInt64]()
    for m in range(len(found.r)):
        # In left's letters right runs from r - q (before left's start when a right read reaches past the left
        # anchor) to r - q + its length, which must reach left's end (or the join would drop left's last letters).
        var start = found.r[m] - found.q[m]
        var end = start + len(right)
        if unique[m] and end >= len(left) and min(len(left), end) - max(0, start) >= JOIN_LEAST:
            keys.append((UInt64(start + (1 << 30)) << 32) | UInt64(m))
    sort(keys)
    var k = 0
    while k < len(keys):
        # One group: pins whose diagonals are within JOIN_DRIFT of the one before.
        var group = List[Int]()
        var last = Int(keys[k] >> 32) - (1 << 30)
        while k < len(keys) and Int(keys[k] >> 32) - (1 << 30) - last <= JOIN_DRIFT:
            last = Int(keys[k] >> 32) - (1 << 30)
            group.append(Int(keys[k] & 0xFFFFFFFF))
            k += 1
        # Its join: the pin best supported on both sides; as well supported, the one nearer the group's middle.
        var best = -1
        var best_support = -1
        var best_off = 0
        for g in range(len(group)):
            var m = group[g]
            var support = min(_support(left_reach, found.r[m]), _support(right_reach, len(right) - found.q[m]))
            var off = abs(g - len(group) // 2)
            if support > best_support or (support == best_support and off < best_off):
                best = m
                best_support = support
                best_off = off
        try:
            placings.append(
                _placing(
                    aligner,
                    left,
                    right,
                    found.r[best],
                    found.q[best],
                    String(
                        "the pin left ",
                        found.r[best],
                        " / right ",
                        found.q[best],
                        " (",
                        len(group),
                        " pins; ",
                        best_support,
                        " reads either side)",
                    ),
                )
            )
        except:
            pass
    # The overlaps of the ends.
    var a_start = max(0, len(left) - JOIN_ENDS)
    var a = List[UInt8](capacity=len(left) - a_start)
    for i in range(a_start, len(left)):
        a.append(left[i])
    var b = List[UInt8](capacity=min(len(right), JOIN_ENDS))
    for j in range(min(len(right), JOIN_ENDS)):
        b.append(right[j])
    var mode = Mode.ends_free(reference_start=len(a), query_end=len(b)).with_match_score(1)
    for edit in [JOIN_EDIT, JOIN_EDIT_STRICT]:
        var got = aligner.align(Span(a), Span(b), Costs.linear(edit, edit), mode)
        if got.reference_end != len(a) or got.query_start != 0 or got.query_end - got.query_start < JOIN_LEAST:
            continue
        try:
            placings.append(
                _placing(
                    aligner,
                    left,
                    right,
                    len(left),
                    got.query_end,
                    String("an overlap of the ends, ", got.query_end, " letters (edits ", edit, ")"),
                )
            )
        except:
            pass
    var best = _best(placings)
    if best >= 0:
        # Inside a tandem array the chain follows one copy, and the placement a copy or more off may be no candidate
        # at all: with the overlap periodic, every shift of the best by whole periods is one (each copy number the
        # two sides could have between them), scored alike.
        var at = placings[best].diagonal
        var period = _period(left, max(0, at), len(left))
        if period > 0:
            # The candidates are now the copy numbers, every one scored alike (the best's own placement too): its
            # meeting point moved by whole copies along right, to the same place in that copy (`_same_phase`), and
            # aligned out from there as the best was, the far ends free (the consensus sequences' indels drift a long
            # overlap's ends by more than a copy).
            var cut_left = placings[best].cut_left
            var cut_right = placings[best].cut_right
            var from_how = placings[best].how.copy()
            var hypotheses = List[_Placing]()
            for k in range(-JOIN_SHIFTS, JOIN_SHIFTS + 1):
                var meets = cut_right if k == 0 else _same_phase(right, cut_right, cut_right - k * period, period // 2)
                var diagonal = cut_left - meets
                if meets < 0 or meets > len(right) or diagonal + len(right) < len(left):
                    continue
                if len(left) - max(0, diagonal) < JOIN_LEAST:
                    continue
                try:
                    hypotheses.append(
                        _placing(
                            aligner,
                            left,
                            right,
                            cut_left,
                            meets,
                            String("copy shift ", k, " (period ", period, ") of ", from_how),
                        )
                    )
                except:
                    pass
            if len(hypotheses) > 0:
                return _copy_join(left, right, hypotheses, period, trace)
            best = _best(placings)
    if best < 0:
        trace.line(
            POA,
            String(
                "not joined: ", len(placings), " placings, none ", JOIN_IDENTITY, " alike over ", JOIN_LEAST, " letters"
            ),
        )
        return List[UInt8]()
    ref chosen = placings[best]
    # Another placement whose meeting point lies on the chosen one's path (within JOIN_DRIFT) is the same one, its
    # diagonal moved by the two sides' insertions and deletions. (An overlap of the ends meets at left's end: its
    # last letter, one right letter before.)
    for c in range(len(placings)):
        ref p = placings[c]
        var at_end = p.cut_left >= len(left)
        var meets = chosen.lands(p.cut_left - 1 if at_end else p.cut_left)
        var same = abs(meets - (p.cut_right - 1 if at_end else p.cut_right)) <= JOIN_DRIFT
        if not same and p.letters >= JOIN_LEAST:
            if p.rate() <= chosen.rate() * JOIN_MARGIN + JOIN_MARGIN_FLOOR:
                trace.line(
                    POA,
                    String(
                        "not joined: ",
                        chosen.how,
                        " at ",
                        chosen.rate(),
                        " edits a letter, ",
                        p.how,
                        " ",
                        p.diagonal - chosen.diagonal,
                        " letters off at ",
                        p.rate(),
                    ),
                )
                return List[UInt8]()
    trace.line(POA, String("join at ", chosen.how, ": ", chosen.edits, " edits in ", chosen.letters, " letters"))
    var out = List[UInt8](capacity=chosen.cut_left + len(right) - chosen.cut_right)
    for i in range(chosen.cut_left):
        out.append(left[i])
    for j in range(chosen.cut_right, len(right)):
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


def _listed(reads: List[Int]) -> String:
    """`reads` as a comma-separated list, for a trace line."""
    var out = String()
    for k in range(len(reads)):
        out.write("," if k > 0 else "", reads[k])
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
