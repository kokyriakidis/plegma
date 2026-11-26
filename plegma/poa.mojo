"""Partial order alignment of a set of reads, as `plegma poa` runs it, reads anchored at both ends or at one.

Anchored-read quality checks and reach/trimming live in `anchored.mojo`; two-side joining lives in `join.mojo`;
this module manages `PoaOptions`, `Poa`, reusable `PoaScratch`, and graph construction (`poa`, `poa_with_scratch`).
"""

from std.math import isfinite

from karan import Aligner
from karan.common import spread

from .alphabet import reverse_complement
from .anchored import (
    FOREIGN_SPANNING,
    _STRAY,
    as_built,
    ends_checked,
    find_foreign,
    read_reach,
    trim_consensus,
    without_reads,
)
from .graph import Graph
from .join import join_sides
from .kernel import (
    BOTH,
    EXTENSION,
    GLOBAL,
    LEFT,
    LOCAL,
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
from .pinned import WindowScratch, align_pinned
from .trace import POA, Trace

comptime INPUT_ORDER = 0
comptime GUIDE_ORDER = 1
comptime LONGEST_FIRST = 2

comptime ANCHOR_SPREAD = 16
"""A one-sided assembly's anchored end is voted on over this many nodes (`Graph._spread_vote`)."""
comptime BAND_MAX = 64
"""The adaptive band's reach at most, either side."""


struct PoaOptions(Copyable, Movable):
    """How `poa` aligns: the kernel's `mode` (`GLOBAL`, `LOCAL` or `EXTENSION` for spanning reads), its `scoring`
    and band, the spanning reads' `order`, whether each read is pinned to the consensus (`seeded`), whether a read
    scoring poorly is tried reverse-complemented (`strand`, spanning reads only), `min_cover` for one-directional
    trimming, and `threads` for parallel stages (1: caller's thread only)."""

    var mode: Int
    var scoring: Scoring
    var band_b: Int
    var band_f: Float64
    var band_max: Int
    var order: Int
    var seeded: Bool
    var strand: Bool
    var min_cover: Int
    var threads: Int

    def __init__(out self):
        self.mode = GLOBAL
        self.scoring = Scoring()
        self.band_b = 10
        self.band_f = 0.01
        self.band_max = BAND_MAX
        self.order = GUIDE_ORDER
        self.seeded = True
        self.strand = False
        self.min_cover = 5
        self.threads = 1

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
        if self.threads < 1:
            raise Error(t"threads {self.threads}: must be 1 or more")


struct Poa(Movable):
    """What `poa` made: the graph, output path indices `shown`, `spanning`, `consensus`, `foreign`, `reversed`,
    and `anchors` as built."""

    var graph: Graph
    var shown: List[Int]
    var spanning: List[Int]
    var consensus: List[List[UInt8]]
    var foreign: List[Int]
    var reversed: List[Int]
    var anchors: List[Int]

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


struct PoaScratch(Movable):
    """Reusable scratch buffers for `poa`: kernel workspace, topology buffers, pinned-window scratch, and
    one-sided reachability buffers."""

    var row_of: List[Int]
    var scratch: WindowScratch
    var workspace: Workspace
    var topology: Topology
    var topo_in: List[Int32]
    var topo_ready: List[Int]
    var topo_full: List[Int]
    var allowed: List[Bool]
    var bfs_work: List[Int]

    def __init__(out self):
        self.row_of = List[Int]()
        self.scratch = WindowScratch()
        self.workspace = Workspace()
        self.topology = Topology()
        self.topo_in = List[Int32]()
        self.topo_ready = List[Int]()
        self.topo_full = List[Int]()
        self.allowed = List[Bool]()
        self.bfs_work = List[Int]()

    def order_of(mut self, graph: Graph) raises -> ref[self.topology.order] List[Int]:
        self.topology.compute(graph, self.topo_in, self.topo_ready)
        return self.topology.order

    def align(
        mut self, graph: Graph, seq: List[UInt8], options: PoaOptions, mode: Int, anchors: List[Int] = List[Int]()
    ) raises -> GraphAlignment:
        """`seq` aligned to `graph` in `mode`, reusing all internal topology and workspace buffers."""
        self.row_of.resize(len(graph), -1)
        if len(anchors) > 0 and (mode == PREFIX or mode == SUFFIX):
            _between_into(graph, anchors, mode, self.allowed, self.bfs_work)
            if options.seeded:
                self.topology.compute(graph, self.topo_in, self.topo_ready)
                var found = align_pinned(
                    graph,
                    self.topology,
                    seq,
                    options.scoring,
                    options.band_b,
                    options.band_f,
                    self.row_of,
                    self.scratch,
                    self.workspace,
                    mode,
                    options.threads,
                )
                _past_anchor(found, self.allowed)
                return found^
            else:
                self.topology.compute_allowed(graph, self.allowed, self.topo_in, self.topo_ready, self.topo_full)
                return self.workspace.align(
                    graph,
                    self.topology,
                    self.topology.order,
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
        self.topology.compute(graph, self.topo_in, self.topo_ready)
        if options.seeded and (mode == GLOBAL or mode == PREFIX or mode == SUFFIX):
            return align_pinned(
                graph,
                self.topology,
                seq,
                options.scoring,
                options.band_b,
                options.band_f,
                self.row_of,
                self.scratch,
                self.workspace,
                mode,
                options.threads,
            )
        return self.workspace.align(
            graph,
            self.topology,
            self.topology.order,
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


comptime _Aligner = PoaScratch


def poa(
    mut names: List[String],
    mut seqs: List[List[UInt8]],
    anchors: List[Int],
    options: PoaOptions,
    var start: Graph,
) raises -> Poa:
    """The POA of `seqs` (`names` theirs), each anchored as `anchors` says (`BOTH`, `LEFT` or `RIGHT`), into
    `start` (a graph made from an existing MSA, or empty)."""
    var aligner = PoaScratch()
    var karan_aligner = Aligner(warm=False)
    return poa_with_scratch(names, seqs, anchors, options, start^, aligner, karan_aligner)


def poa_with_scratch(
    mut names: List[String],
    mut seqs: List[List[UInt8]],
    anchors: List[Int],
    options: PoaOptions,
    var start: Graph,
    mut aligner: PoaScratch,
    mut karan_aligner: Aligner,
) raises -> Poa:
    """As `poa`, reusing `aligner` and `karan_aligner` across calls."""
    if len(seqs) == 0 and len(start.paths) == 0:
        return Poa(start^, List[Int](), List[Int]())
    var anchor_of = ends_checked(karan_aligner, seqs, anchors, options.threads)
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
    var graph = start^
    if len(spanning) > 0 or len(graph.paths) > 0:
        var reversed = _add_spanning(
            graph, names, seqs, _ordered(seqs, spanning, options.order, options.threads), options, path_of, aligner
        )
        var members = List[Int]()
        for p in range(starting):
            members.append(p)
        for r in spanning:
            members.append(path_of[r])
        var foreign = List[Int]()
        if len(spanning) >= FOREIGN_SPANNING and len(left) + len(right) > 0:
            var core = _letters(graph, graph.heaviest_path_spanning(aligner.order_of(graph), members))
            foreign = find_foreign(karan_aligner, core, seqs, spanning, left, right, options.threads)
            if trace.on(POA) and len(foreign) > 0:
                trace.line(POA, String("another copy's (no vote): reads ", _listed(foreign)))
        var kept_left = without_reads(left, foreign)
        var kept_right = without_reads(right, foreign)
        _add_anchored(graph, names, seqs, _longest_first(seqs, kept_left), options, PREFIX, path_of, aligner, members)
        _add_anchored(graph, names, seqs, _longest_first(seqs, kept_right), options, SUFFIX, path_of, aligner, members)
        var consensus = List[List[UInt8]]()
        if len(left) + len(right) > 0:
            consensus.append(_letters(graph, graph.heaviest_path_spanning(aligner.order_of(graph), members)))
        _add_anchored(
            graph,
            names,
            seqs,
            _longest_first(seqs, without_reads(left, kept_left)),
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
            _longest_first(seqs, without_reads(right, kept_right)),
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
        out.anchors = as_built(anchor_of)
        sort(reversed)
        if trace.on(POA) and len(reversed) > 0:
            trace.line(POA, String("reverse-complemented: reads ", _listed(reversed)))
        out.reversed = reversed^
        return out^

    # No read spans both anchors: each side on its own (in parallel when both sides exist and threads > 1), then joined.
    var parts = List[List[UInt8]]()
    var left_reach = List[Int]()
    var right_reach = List[Int]()
    if options.threads > 1 and len(left) > 0 and len(right) > 0:
        var side_parts = List[List[UInt8]](length=2, fill=List[UInt8]())
        var side_reaches = List[List[Int]](length=2, fill=List[Int]())
        var side_failed = List[Bool](length=2, fill=False)
        var sp_p = side_parts.unsafe_ptr()
        var sr_p = side_reaches.unsafe_ptr()
        var sf_p = side_failed.unsafe_ptr()
        var side_options = options.copy()
        side_options.threads = max(1, options.threads // 2)

        def build_side(
            which: Int,
        ) {imm names, imm seqs, imm left, imm right, imm side_options, imm sp_p, imm sr_p, imm sf_p}:
            var local_aligner = PoaScratch()
            try:
                if which == 0:
                    var side = _one_side(names, seqs, left, side_options, PREFIX, local_aligner)
                    var path = side[0].heaviest_path(local_aligner.order_of(side[0]), side[1], ANCHOR_SPREAD)
                    sp_p[unsafe_offset=0] = _letters(side[0], path)
                    sr_p[unsafe_offset=0] = read_reach(side[0], path, side[1], True)
                else:
                    var side = _one_side(names, seqs, right, side_options, SUFFIX, local_aligner)
                    var path = side[0].heaviest_path_backward(local_aligner.order_of(side[0]), side[1], ANCHOR_SPREAD)
                    sp_p[unsafe_offset=1] = _letters(side[0], path)
                    sr_p[unsafe_offset=1] = read_reach(side[0], path, side[1], False)
            except:
                sf_p[unsafe_offset=which] = True

        spread(build_side, 2, 2)
        if side_failed[0]:
            var side = _one_side(names, seqs, left, options, PREFIX, aligner)
            var path = side[0].heaviest_path(aligner.order_of(side[0]), side[1], ANCHOR_SPREAD)
            side_parts[0] = _letters(side[0], path)
            side_reaches[0] = read_reach(side[0], path, side[1], True)
        if side_failed[1]:
            var side = _one_side(names, seqs, right, options, SUFFIX, aligner)
            var path = side[0].heaviest_path_backward(aligner.order_of(side[0]), side[1], ANCHOR_SPREAD)
            side_parts[1] = _letters(side[0], path)
            side_reaches[1] = read_reach(side[0], path, side[1], False)
        parts.append(side_parts[0].copy())
        parts.append(side_parts[1].copy())
        left_reach = side_reaches[0].copy()
        right_reach = side_reaches[1].copy()
    else:
        if len(left) > 0:
            var side = _one_side(names, seqs, left, options, PREFIX, aligner)
            var path = side[0].heaviest_path(aligner.order_of(side[0]), side[1], ANCHOR_SPREAD)
            parts.append(_letters(side[0], path))
            left_reach = read_reach(side[0], path, side[1], True)
        if len(right) > 0:
            var side = _one_side(names, seqs, right, options, SUFFIX, aligner)
            var path = side[0].heaviest_path_backward(aligner.order_of(side[0]), side[1], ANCHOR_SPREAD)
            parts.append(_letters(side[0], path))
            right_reach = read_reach(side[0], path, side[1], False)
    if len(parts) == 2:
        var joined = join_sides(karan_aligner, parts[0], parts[1], left_reach, right_reach, options.threads)
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
        trim_consensus(parts[0], left_reach, options.min_cover, True)
        trace.line(
            POA,
            String("left side alone: ", len(parts[0]), " letters, cut where under ", options.min_cover, " reads go on"),
        )
    else:
        trim_consensus(parts[0], right_reach, options.min_cover, False)
        trace.line(
            POA,
            String(
                "right side alone: ", len(parts[0]), " letters, cut where under ", options.min_cover, " reads go on"
            ),
        )
    _ = graph.add("plegma_backbone", Span(parts[0]))
    if len(parts) == 1:
        var backbone: List[Int] = [0]
        _add_anchored(graph, names, seqs, _longest_first(seqs, left), options, PREFIX, path_of, aligner, backbone)
        _add_anchored(graph, names, seqs, _longest_first(seqs, right), options, SUFFIX, path_of, aligner, backbone)
        var path = graph.heaviest_path_spanning(aligner.order_of(graph), backbone, len(left) > 0 and len(right) > 0)
        parts[0] = _letters(graph, path)
        if len(right) == 0 or len(left) == 0:
            var paths = List[Int]()
            for r in left if len(right) == 0 else right:
                paths.append(path_of[r])
            trim_consensus(
                parts[0], read_reach(graph, path, paths, len(right) == 0), options.min_cover, len(right) == 0
            )
    else:
        _add_anchored(graph, names, seqs, _longest_first(seqs, left), options, PREFIX, path_of, aligner)
        var order = _longest_first(seqs, right)
        if len(order) > 0:
            path_of[order[0]] = graph.add(names[order[0]], Span(seqs[order[0]]))
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
    out.anchors = as_built(anchor_of)
    return out^


def _add_strays(
    mut graph: Graph,
    names: List[String],
    seqs: List[List[UInt8]],
    strays: List[Int],
    options: PoaOptions,
    mut path_of: List[Int],
    mut aligner: PoaScratch,
) raises:
    for r in strays:
        var found = aligner.align(graph, seqs[r], options, options.mode)
        path_of[r] = graph.fuse_nodes(names[r], Span(seqs[r]), found.nodes)


def _add_spanning(
    mut graph: Graph,
    mut names: List[String],
    mut seqs: List[List[UInt8]],
    order: List[Int],
    options: PoaOptions,
    mut path_of: List[Int],
    mut aligner: PoaScratch,
) raises -> List[Int]:
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
            var floor = Int(options.scoring.match_score) * min(len(seqs[r]), len(graph.paths[0])) // 3
            if found.score < floor:
                var flipped = reverse_complement(seqs[r])
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
) -> Bool:
    """Read `r` fused along the path of an earlier read with the same letters (in `seen`), if there is one."""
    var key = String(unsafe_from_utf8=seqs[r])
    var twin = seen.get(key, -1)
    if twin < 0:
        return False
    path_of[r] = graph.add_twin_path(names[r], path_of[twin])
    return True


def _add_anchored(
    mut graph: Graph,
    names: List[String],
    seqs: List[List[UInt8]],
    order: List[Int],
    options: PoaOptions,
    mode: Int,
    mut path_of: List[Int],
    mut aligner: PoaScratch,
    anchors: List[Int] = List[Int](),
) raises:
    var seen = Dict[String, Int]()
    for r in order:
        if _twin(graph, names, seqs, r, seen, path_of):
            continue
        var found = aligner.align(graph, seqs[r], options, mode, anchors)
        path_of[r] = graph.fuse_nodes(names[r], Span(seqs[r]), found.nodes)
        seen[String(unsafe_from_utf8=seqs[r])] = r


def _between_into(graph: Graph, anchors: List[Int], mode: Int, mut allowed: List[Bool], mut work: List[Int]):
    """Writes into `allowed` the nodes a read anchored at one end may align to, reusing `allowed` and `work`."""
    var n = len(graph)
    allowed.resize(n, False)
    for i in range(n):
        allowed[i] = False
    work.clear()
    for p in anchors:
        ref path = graph.paths[p]
        if len(path) == 0:
            continue
        var node = Int(path[len(path) - 1]) if mode == SUFFIX else Int(path[0])
        if not allowed[node]:
            allowed[node] = True
            work.append(node)
    while len(work) > 0:
        var node = work.pop()
        var e = Int(graph.first_in[node]) if mode == SUFFIX else Int(graph.first_out[node])
        while e >= 0:
            var other = Int(graph.edge_from[e]) if mode == SUFFIX else Int(graph.edge_to[e])
            if not allowed[other]:
                allowed[other] = True
                work.append(other)
            e = Int(graph.next_in[e]) if mode == SUFFIX else Int(graph.next_out[e])


def _past_anchor(mut found: GraphAlignment, allowed: List[Bool]):
    for i in range(len(found.nodes)):
        var node = found.nodes[i]
        if node >= 0 and not allowed[node]:
            found.nodes[i] = -1


def _one_side(
    names: List[String],
    seqs: List[List[UInt8]],
    reads: List[Int],
    options: PoaOptions,
    mode: Int,
    mut aligner: PoaScratch,
) raises -> Tuple[Graph, List[Int]]:
    """The graph of one side's reads alone: the longest first, the rest anchored as `mode` says; with its paths."""
    var graph = Graph()
    var order = _longest_first(seqs, reads)
    var paths = List[Int]()
    for k in range(len(order)):
        var r = order[k]
        if k == 0:
            paths.append(graph.add(names[r], Span(seqs[r])))
            continue
        var found = aligner.align(graph, seqs[r], options, mode)
        paths.append(graph.fuse_nodes(names[r], Span(seqs[r]), found.nodes))
    return (graph^, paths^)


def _ordered(seqs: List[List[UInt8]], reads: List[Int], order: Int, threads: Int = 1) raises -> List[Int]:
    """`reads` in input order, abPOA's guide order or longest first."""
    if order == GUIDE_ORDER and len(reads) > 1:
        var subset = List[List[UInt8]]()
        for r in reads:
            subset.append(seqs[r].copy())
        var out = List[Int]()
        for k in guide_order(subset, threads=threads):
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
    var out = String()
    for k in range(len(reads)):
        out.write("," if k > 0 else "", reads[k])
    return out^


def _letters(graph: Graph, path: List[Int]) -> List[UInt8]:
    var out = List[UInt8](capacity=len(path))
    for node in path:
        out.append(graph.bases[node])
    return out^
