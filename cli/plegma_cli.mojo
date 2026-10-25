"""The `plegma` command line: a partial order alignment of the sequences of one FASTA/FASTQ file, as abPOA runs
one (its flags, defaults and output modes), written as a consensus, an MSA or a GFA graph."""

from std.os import getenv
from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from karan import Aligner, Costs, Mode, Ties
from sporos import read_fasta

from plegma.anchors import Anchor
from plegma.backbone import Level, MUMS, backbone, split_anchors, windows
from plegma.chain import ChainCosts
from plegma.graph import DELETION, INSERTION, Graph, GAP, cigar_runs
from plegma.clusters import cluster_rows
from plegma.order import guide_order
from plegma.kernel import EXTENSION, GLOBAL, LOCAL, GraphAlignment, Scoring, Topology, Workspace
from plegma.poa import SINK, SOURCE, WindowScratch, window_nodes
from plegma.polish import polish
from plegma.progressive import on_reference, on_tree, progressive
from plegma.refine import refine


comptime USAGE = """Usage: plegma [options] <in.fa/fq> > cons.fa/msa.fa/plegma.gfa

Options:
  Alignment:
    -m INT      alignment mode [0]: 0 global, 1 local, 2 extension
    -M INT      match score [2]
    -X INT      mismatch penalty [4]
    -O INT[,INT] gap opening penalties O1,O2 [4,24]
    -E INT[,INT] gap extension penalties E1,E2 [2,1]
                a g-long gap costs min{O1+g*E1, O2+g*E2} (convex); O2 0: O1+g*E1 (affine); O1 0: g*E1 (linear)
  Adaptive banded DP:
    -b INT      extra band b [10]; < 0 disables the band
    -f FLOAT    extra band f [0.01]: b+f*L letters either side, L the sequence's length
  Seeding and partition (global mode):
    -S          anchor each sequence to the consensus path with sporos's MUMs (51, then 31, then 21 letters,
                seeded again inside every gap, none inside tandem repeats), and align only the windows
                between anchors
    -J          with -S, anchors by split_anchors: the flanks' MUMs, then inside tandem repeats the longest match
                unique in each window that leaves its length difference to one side and holds enough letters for it
    -Y          with -S, anchors where the sequences already in the graph agree: the sequence aligned optimally
                (unit costs) to the 24 nearest, each matched letter voting for the node that one spells it with, and
                every run of 21 or more letters a majority of them place on consecutive nodes held
    -n INT      gaps of at most INT letters are not seeded again [500]
    -T INT      letters each anchor gives back to the DP at both ends [10]
    -H INT      anchors each anchor's chain search looks back over, at most [200]; 0: all
    -K INT:INT  k:d, no anchor where one of its k-mers recurs within d letters of the path (a tandem
                repeat) [12:500]; 0:0 keeps them
  Order and strand:
    -L          align the sequences longest first
    -p          progressive: abPOA's guide order, most similar minimizer multisets first
    -s          ambiguous strand: a sequence scoring poorly is tried reverse-complemented, the better kept
  Refinement:
    -R INT      rounds of realignment [0]: each round starts a graph from the consensus so far and aligns every
                sequence to it again; the consensus is then the sequences' own (the backbone carries no weight)
    -I INT[,INT[,INT]] passes of sum-of-pairs refinement of the MSA [0]: each distinct sequence realigned in turn to
                the others' columns, exactly, in a band of the second number of letters [64] centred where the
                sequences' optimal pairwise alignments place it, every gap a pair opens costing the third [1];
                the fourth 1 seeds those pairwise alignments (flanks pinned by MUMs, then MUMs inside each array)
                [0]; the fifth caps the band's widening [4096] (MSA and GFA outputs)
    -G INT      the -I refinement's group moves: every split of a UPGMA tree over the rows realigned as two
                profiles within INT columns of their present merge, each pass [0; 256 with -A -N]; 0 off
    -A          the MSA instead of the graph, anchored on the first sequence (the reference): every other
                aligned to it exactly (minimap2's asm5 scores as costs, gaps left-normalized), its letters in the
                reference's columns and the runs inserted at one place aligned to each other on the longest; each
                sequence's alignment with the reference is then its own optimal one, as a caller reporting variants
                against the reference needs; with -I, then refined on the sum of pairs (-W, -I's prices) with the
                reference weighing 4 times the rows' number, so the others choose among the placements as good
                against it the one that suits each other; the graph is the MSA's
    -N          with -A, placement by votes instead (no reference): the distinct sequences, longest first, each
                placed exactly into the MSA so far in a band centred on the columns every placed sequence votes for
                its letters through their optimal alignment, at vntr-eval's affine costs (with -V, through seeds
                first: sporos's MUMs outside tandem repeats, then inside each array, then SMEMs); where every
                sequence holds one clear tandem array, a second start with the arrays aligned copy by copy; each
                refined (-I, 20 passes by default, with -G's group moves, at -W's prices) and the lower sum of pairs
                kept (lowest on vntr-eval's costs; on chr20 its graphs genotype worse, splitting copy-number changes)
    -Q          with -A, the default (kept for older command lines)
    -t INT      threads for -A's alignments to the reference [1]
    -U          with -A -N, no copy-level start
    -W INT,INT  the sum-of-pairs prices of a column of two different letters and of a letter against a gap [1,1;
                2,2 with -A] (-I's third number prices a gap opening, 1): vntr-eval's unit cost is 1,1 with no
                opening, its affine cost 4,2 with 6
    -V INT[,INT...] with -A -N, the seeds' confidence [-2,3,0,0,0]: -2 no seeds (the optimal alignments
                to the placed sequences vote, the staged seeds where a pair is too long to align); -4 the staged
                seeds first (MUMs outside tandem repeats pin the flanks, then MUMs of 1,000, 600 and 400 letters),
                then the alignments; -5 the split seeds first (backbone.split_anchors); -3 SMEMs of 600 letters or more
                first; e >= 0 sporos's MUMs and SMEMs, an
                anchor inside a repeat kept only if its diagonal lies within e letters of one end of its gap (one
                indel), -1 without that rule; the second, the costs of the optimal pairwise alignments, which also
                centre the refinement's bands: 3 vntr-eval's affine (4, 6 + 2k), 1 unit (its unit metric), 2
                gap-affine at -W's prices, 0 minimap2's two-piece; the third holds the letters that many
                of the best sources anchor to their columns in placement [0]; a fourth 1 places the sequences nearest
                first, not longest first; a fifth 1 anchors a letter where most of its votes agree (not more than
                half), each anchor weighing the votes its column won; a sixth n > 0 places them in a random order
                (after the longest) seeded by n, a restart; a seventh n > 0 places each sequence exactly (the DP
                unbanded) where that takes n million cells or fewer; an eighth n > 0 runs a refinement pass over the
                sequences placed so far after every n placements; a ninth s >= 0 keeps each placement inside its
                co-optimal corridor against its nearest placed relative (alignments within s of their optimum); a
                tenth 1 adds the subfamily start (the two clusters at the top of a tree aligned each on its own, then
                merged), an eleventh 1 the tree start (along the minimum spanning tree)
    -B          refinement starts each round's graph from the consensus polished (-P), not the heaviest path
    -P INT[,INT] rounds of consensus polishing, at most [50,20] (global mode): each sequence aligned to the
                consensus, and every one-letter or whole-repeat-unit edit kept that lowers their summed edit distance
                to it (the first), then that raises their likelihood under a pair HMM (the second); 0,0 off
  Input/Output:
    -l          the input is a list of sequence files, one per line, each aligned on its own
    -i FILE     start from the graph of an existing MSA (FASTA, rows with '-')
    -o FILE     output to FILE [stdout]
    -r INT      output [0]: 0 consensus FASTA, 1 MSA FASTA, 2 both, 3 GFA, 4 GFA with the consensus path
    -a INT      consensus [0]: 0 heaviest path, 1 most frequent letter per column
    -d INT      at most INT consensus sequences [1]; at 2 or more the sequences are clustered on the MSA's
                heterozygous columns and each cluster gets its own
    -q FLOAT    smallest share of the sequences a cluster and a column's second allele must hold [0.25]
    -v          report timings on stderr
    -h          this help
"""


def _value(args: List[String], i: Int) raises -> String:
    if i + 1 >= len(args):
        raise Error(t"{args[i]} takes a value")
    return args[i + 1]


def _pair(text: String) raises -> Tuple[Int, Int]:
    var parts = text.split(",")
    if len(parts) == 1:
        return (Int(parts[0]), -1)
    return (Int(parts[0]), Int(parts[1]))


def _consensus_most_frequent(rows: List[String]) -> String:
    var out = List[UInt8]()
    if len(rows) == 0:
        return String()
    var width = rows[0].byte_length()
    var counts = List[Int](length=256, fill=0)
    for c in range(width):
        for b in range(256):
            counts[b] = 0
        for r in rows:
            counts[Int(r.as_bytes()[c])] += 1
        var best = Int(GAP)
        for b in range(256):
            if counts[b] > counts[best]:
                best = b
        if best != Int(GAP):
            out.append(UInt8(best))
    return String(unsafe_from_utf8=out)


comptime AGREE_SHARE = 0.5
"""-Y: the share of the voters a letter's node needs."""
comptime AGREE_RUN = 21
"""-Y: the letters an agreed run needs to anchor."""
comptime AGREE_VOTERS = 24
"""-Y: the sequences in the graph that vote, the nearest by edit distance."""


def _agreed(
    graph: Graph,
    topology: Topology,
    query: Span[UInt8, _],
    share: Float64,
    least: Int,
    voters: Int,
    mut path: List[Int],
) raises -> List[Anchor]:
    """Anchors where the sequences already in the graph agree: `query` aligned optimally (unit costs) to each of
    the `voters` nearest of them (all when fewer), each matched letter voting for the node that sequence spells it
    with; a letter whose most voted node holds `share` of the votes or more takes it, and a run of `least` or more
    such letters, each node an edge after the last and the run after the one before in the graph's order, is an
    anchor. `path` gets the anchors' nodes one after another, the anchors' positions along it."""
    var aligner = Aligner()
    var members = List[Int]()
    if len(graph.paths) <= voters:
        for p in range(len(graph.paths)):
            members.append(p)
    else:
        var keys = List[UInt64]()
        for p in range(len(graph.paths)):
            var letters = List[UInt8](capacity=len(graph.paths[p]))
            for node in graph.paths[p]:
                letters.append(graph.bases[Int(node)])
            keys.append((UInt64(aligner.distance(Span(letters), query)) << 32) | UInt64(p))
        sort(keys)
        for t in range(voters):
            members.append(Int(keys[t] & 0xFFFFFFFF))
    var n = len(query)
    var k = len(members)
    var node_of = List[Int](length=k * n, fill=-1)
    for v in range(k):
        ref member = graph.paths[members[v]]
        var letters = List[UInt8](capacity=len(member))
        for node in member:
            letters.append(graph.bases[Int(node)])
        var found = aligner.align(Span(letters), query, Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
        var runs = cigar_runs(found.cigar)
        var i = 0
        var j = 0
        for x in range(len(runs.ops)):
            for _ in range(runs.counts[x]):
                var op = runs.ops[x]
                if op == DELETION:
                    i += 1
                elif op == INSERTION:
                    j += 1
                else:
                    if letters[i] == query[j]:
                        node_of[v * n + j] = Int(member[i])
                    i += 1
                    j += 1
    var need = max(1, Int(share * Float64(k) + 0.999999))
    var lead = List[Int](length=n, fill=-1)
    for j in range(n):
        var top = -1
        var top_count = 0
        for v in range(k):
            var node = node_of[v * n + j]
            if node < 0 or node == top:
                continue
            var c = 0
            for u in range(k):
                if node_of[u * n + j] == node:
                    c += 1
            if c > top_count:
                top_count = c
                top = node
        if top_count >= need:
            lead[j] = top
    var out = List[Anchor]()
    var last_rank = -1
    var start = 0
    for j in range(n + 1):
        var extends = False
        if j < n and lead[j] >= 0 and j > start:
            var e = Int(graph.first_out[lead[j - 1]])
            while e >= 0:
                if Int(graph.edge_to[e]) == lead[j]:
                    extends = True
                    break
                e = Int(graph.next_out[e])
        if j < n and lead[j] >= 0 and (j == start or extends):
            continue
        # The run start ..< j ends here.
        if j - start >= least and lead[start] >= 0 and topology.rank[lead[start]] > last_rank:
            out.append(Anchor(len(path), start, j - start))
            for u in range(start, j):
                path.append(lead[u])
            last_rank = topology.rank[lead[j - 1]]
        start = j if j < n and lead[j] >= 0 else j + 1
    return out^


def _seeded(
    graph: Graph,
    topology: Topology,
    query: Span[UInt8, _],
    scoring: Scoring,
    band_b: Int,
    band_f: Float64,
    window: Int,
    trim: Int,
    lookback: Int,
    repeat_k: Int,
    repeat_distance: Int,
    split: Bool,
    agree: Bool,
    mut row_of: List[Int],
    mut scratch: WindowScratch,
    mut workspace: Workspace,
) raises -> GraphAlignment:
    """`query` aligned to the graph through anchors on the consensus path: each window between two anchors
    aligned to the subgraph between their nodes, the anchors' letters matched."""
    var path = graph.heaviest_path(topology.order)
    var agreed = List[Anchor]()
    var by_agreement = agree and len(graph.paths) > 0
    if by_agreement:
        path = List[Int]()
        agreed = _agreed(graph, topology, query, AGREE_SHARE, AGREE_RUN, AGREE_VOTERS, path)
    var letters = List[UInt8](capacity=len(path))
    for node in path:
        letters.append(graph.bases[node])
    # MUMs, placed by uniqueness, of 51 letters, then 31, then 21, each gap seeded again until it yields nothing;
    # none inside a tandem repeat of the path (its k-mers recurring nearby), where even a unique exact match may
    # pair different copies: there the DP aligns the whole array.
    var levels = List[Level]()
    for n in [51, 31, 21]:
        levels.append(Level(MUMS, n))
    # With `split`, `split_anchors` instead: the flanks, then inside the repeats the anchors it is sure of.
    var chained = agreed.copy() if by_agreement else (
        split_anchors(Span(letters), query) if split else backbone(
            Span(letters),
            query,
            levels,
            ChainCosts(substitution=0.1),
            window,
            True,
            lookback=lookback,
            repeat_k=repeat_k,
            repeat_distance=repeat_distance,
        )
    )
    # Each anchor loses `trim` letters at both ends: the DP, not the anchor, decides where a gap next to it goes.
    var anchors = List[Anchor]()
    for a in chained:
        if a.length > 2 * trim:
            anchors.append(
                Anchor(a.ref_pos + trim, a.query_pos + trim, a.length - 2 * trim, a.ref_copies, a.query_copies)
            )
    # The seed audit (PLEGMA_SEED_DUMP=<file>): each sequence's path, letters, anchors (as aligned, trimmed) and
    # where its letters landed on the path, appended for `anchors --seed-audit` to check against its optimal
    # alignments.
    var dump = getenv("PLEGMA_SEED_DUMP")
    if dump != "":
        var text = String(">path\n")
        text.write(StringSlice(unsafe_from_utf8=Span(letters)), "\n>query\n")
        var q = List[UInt8](capacity=len(query))
        for c in query:
            q.append(c)
        text.write(StringSlice(unsafe_from_utf8=Span(q)), "\n>anchors")
        for a in anchors:
            text.write(" ", a.ref_pos, ",", a.query_pos, ",", a.length)
        text.write("\n")
        with open(dump, "a") as f:
            f.write(text)
    var gaps = windows(anchors, len(letters), len(query))
    var aligned = List[Int](length=len(query), fill=-1)
    var score = 0
    var cells = 0
    var peak = 0
    for w in range(len(gaps)):
        var gap = gaps[w]
        var left = path[gap.ref_start - 1] if gap.ref_start > 0 else SOURCE
        var right = path[gap.ref_end] if gap.ref_end < len(path) else SINK
        if gap.query_end > gap.query_start:
            var nodes = window_nodes(graph, topology.rank, topology.order, left, right, scratch)
            var found = workspace.align(
                graph,
                topology,
                nodes,
                left,
                right,
                query[gap.query_start : gap.query_end],
                scoring,
                row_of,
                GLOBAL,
                band_b,
                band_f,
            )
            for j in range(len(found.nodes)):
                aligned[gap.query_start + j] = found.nodes[j]
            score += found.score
            cells += found.cells
            peak = max(peak, found.peak_bytes)
        if w < len(anchors):
            var a = anchors[w]
            for k in range(a.length):
                aligned[a.query_pos + k] = path[a.ref_pos + k]
            score += Int(scoring.match_score) * a.length
    if dump != "":
        # Each letter's place on the path (-1: aligned off it, or not aligned).
        var at = Dict[Int, Int]()
        for k in range(len(path)):
            at[path[k]] = k
        var text = String(">aligned")
        for j in range(len(query)):
            text.write(" ", at.get(aligned[j], -1) if aligned[j] >= 0 else -1)
        text.write("\n")
        with open(dump, "a") as f:
            f.write(text)
    return GraphAlignment(aligned^, score, cells, peak)


struct Options(Copyable, Movable):
    var mode: Int
    var scoring: Scoring
    var band_b: Int
    var band_f: Float64
    var longest_first: Bool
    var progressive: Bool
    var strand: Bool
    var seeding: Bool
    var split_seeds: Bool
    var voted: Bool
    var threads: Int
    var agree_seeds: Bool
    var window: Int
    var trim: Int
    var lookback: Int
    var repeat_k: Int
    var repeat_distance: Int
    var start_msa: String
    var result: Int
    var consensus_rule: Int
    var max_consensus: Int
    var rounds: Int
    var polish_rounds: Int
    var likely_rounds: Int
    var polish_backbone: Bool
    var sp_passes: Int
    var sp_band: Int
    var sp_gap_open: Int
    var sp_seeded: Bool
    var sp_max_band: Int
    var sp_set: Bool
    var group_band: Int
    var group_set: Bool
    var anchored: Bool
    var units: Bool
    var one_indel: Int
    var votes: Int
    var mismatch: Int
    var gap: Int
    var hard: Int
    var nearest: Bool
    var consistency: Bool
    var shuffle: Int
    var exact_cells: Int
    var interleave: Int
    var corridor: Int
    var subfamilies: Bool
    var tree_start: Bool
    var confidence_set: Bool
    var prices_set: Bool
    var min_freq: Float64
    var verbose: Bool

    def __init__(out self):
        self.mode = GLOBAL
        self.scoring = Scoring()
        self.band_b = 10
        self.band_f = 0.01
        self.longest_first = False
        self.progressive = False
        self.strand = False
        self.seeding = False
        self.split_seeds = False
        self.voted = False
        self.threads = 1
        self.agree_seeds = False
        self.window = 500
        self.trim = 10
        self.lookback = 200
        self.repeat_k = 12
        self.repeat_distance = 500
        self.start_msa = String()
        self.result = 0
        self.consensus_rule = 0
        self.max_consensus = 1
        self.rounds = 0
        self.polish_rounds = 50
        self.likely_rounds = 20
        self.polish_backbone = False
        self.sp_passes = 0
        self.sp_band = 64
        self.sp_gap_open = 1
        self.sp_seeded = False
        self.sp_max_band = 4096
        self.sp_set = False
        self.group_band = 0
        self.group_set = False
        self.anchored = False
        self.units = True
        self.one_indel = -1
        self.votes = 3
        self.mismatch = 1
        self.gap = 1
        self.hard = 0
        self.nearest = False
        self.consistency = False
        self.shuffle = 0
        self.exact_cells = 0
        self.interleave = 0
        self.corridor = -1
        self.subfamilies = False
        self.tree_start = False
        self.confidence_set = False
        self.prices_set = False
        self.min_freq = 0.25
        self.verbose = False


def _reverse_complement(seq: Span[UInt8, _]) -> List[UInt8]:
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


def _align(
    graph: Graph,
    seq: Span[UInt8, _],
    options: Options,
    mut row_of: List[Int],
    mut scratch: WindowScratch,
    mut workspace: Workspace,
) raises -> GraphAlignment:
    var topology = Topology(graph)
    row_of.resize(len(graph), -1)
    if options.seeding and options.mode == GLOBAL:
        return _seeded(
            graph,
            topology,
            seq,
            options.scoring,
            options.band_b,
            options.band_f,
            options.window,
            options.trim,
            options.lookback,
            options.repeat_k,
            options.repeat_distance,
            options.split_seeds,
            options.agree_seeds,
            row_of,
            scratch,
            workspace,
        )
    return workspace.align(
        graph,
        topology,
        topology.order,
        SOURCE,
        SINK,
        seq,
        options.scoring,
        row_of,
        options.mode,
        options.band_b,
        options.band_f,
    )


def _write_gfa(graph: Graph, shown: List[Int], consensus_path: List[Int], mut text: String) raises:
    """The graph as GFA 1.0 with every unbranched run of nodes one segment (a node with one out-edge into a node
    with one in-edge, no path starting or ending between them): the paths in `shown` order, then the consensus
    path when there is one."""
    var n = len(graph)
    var starts = List[Bool](length=n, fill=False)
    var ends = List[Bool](length=n, fill=False)
    for path in graph.paths:
        if len(path) > 0:
            starts[Int(path[0])] = True
            ends[Int(path[len(path) - 1])] = True
    if len(consensus_path) > 0:
        starts[consensus_path[0]] = True
        ends[consensus_path[len(consensus_path) - 1]] = True
    # A node continues its predecessor's segment when nothing branches or stops between them.
    var segment = List[Int32](length=n, fill=-1)
    var heads = List[Int]()
    var topology = graph.node_order()
    var order = List[Int](length=n, fill=0)
    for node in range(n):
        order[topology[node]] = node
    for node in order:
        if graph.in_degree[node] == 1 and not starts[node]:
            var p = Int(graph.edge_from[Int(graph.first_in[node])])
            if graph.out_degree[p] == 1 and not ends[p]:
                segment[node] = segment[p]
                continue
        segment[node] = Int32(len(heads))
        heads.append(node)
    text.write("H\tVN:Z:1.0\n")
    var last = List[Int](length=len(heads), fill=0)
    for k in range(len(heads)):
        var node = heads[k]
        text.write("S\t", k + 1, "\t")
        var letters = List[UInt8]()
        while True:
            letters.append(graph.bases[node])
            last[k] = node
            if graph.out_degree[node] != 1:
                break
            var s = Int(graph.edge_to[Int(graph.first_out[node])])
            if segment[s] != Int32(k) or s == node:
                break
            node = s
        text.write(StringSlice(unsafe_from_utf8=Span(letters)), "\n")
    for k in range(len(heads)):
        var e = Int(graph.first_out[last[k]])
        while e >= 0:
            text.write("L\t", k + 1, "\t+\t", Int(segment[Int(graph.edge_to[e])]) + 1, "\t+\t0M\n")
            e = Int(graph.next_out[e])
    for p in shown:
        text.write("P\t", graph.names[p], "\t")
        _write_steps(graph.paths[p], segment, heads, text)
        text.write("\t*\n")
    if len(consensus_path) > 0:
        var steps = List[Int32]()
        for node in consensus_path:
            steps.append(Int32(node))
        text.write("P\tConsensus_sequence\t")
        _write_steps(steps, segment, heads, text)
        text.write("\t*\n")


def _write_steps(path: List[Int32], segment: List[Int32], heads: List[Int], mut text: String):
    var first = True
    for node in path:
        var k = Int(segment[Int(node)])
        if heads[k] != Int(node):
            continue
        if not first:
            text.write(",")
        text.write(k + 1, "+")
        first = False


def _polished(consensus: String, seqs: List[List[UInt8]], members: List[Int], options: Options) raises -> String:
    if len(members) == 0:
        return consensus
    var letters = List[UInt8](consensus.as_bytes())
    return String(
        unsafe_from_utf8=polish(letters, seqs, members, options.polish_rounds, likely_rounds=options.likely_rounds)
    )


def _poa(input: String, options: Options) raises -> String:
    """The POA of the sequences of one file, as `options` asks for its output."""
    var t0 = perf_counter_ns()
    var records = read_fasta(input)
    var seqs = List[List[UInt8]]()
    for r in range(len(records)):
        # Upper case, as abPOA reads its input: a and A are one letter to the kernel and to the graph alike.
        var letters = List[UInt8](records.seq(r))
        for k in range(len(letters)):
            if letters[k] >= UInt8(ord("a")) and letters[k] <= UInt8(ord("z")):
                letters[k] -= 32
        seqs.append(letters^)
    var order = List[Int]()
    for r in range(len(seqs)):
        order.append(r)
    if options.progressive:
        order = guide_order(seqs)
    elif options.longest_first:
        # Longest first, file order among equals: sort keys (2^32 - length) << 32 | index.
        var keys = List[UInt64]()
        for r in range(len(seqs)):
            keys.append((UInt64((1 << 32) - len(seqs[r])) << 32) | UInt64(r))
        sort(keys)
        for k in range(len(keys)):
            order[k] = Int(keys[k] & 0xFFFFFFFF)
    var t1 = perf_counter_ns()
    var graph = Graph()
    # -A: the MSA by anchored progressive placement; its graph is built from it, and its rows are the output's.
    var anchored_rows = List[String]()
    if options.anchored and not options.voted:
        # Anchored on the first sequence, minimap2's asm5 scores as costs (match 1, mismatch 19, gap opens 39 and 81,
        # extensions 3 and 1, made costs: a mismatch 40, a gap of k letters min(78 + 7k, 162 + 3k)). The truth's
        # convention: on chr20 they genotype better than pgrealign's scores (two-piece 24, 16, 12, 96, 5; structural
        # F1 0.709 against 0.658, small 0.9715 against 0.9704).
        var reference_costs = Costs.two_piece(40, 78, 7, 162, 3)
        # PLEGMA_QCOSTS=x,o,e or x,o1,e1,o2,e2 (an experiment): other costs for the alignments to the reference.
        var q_costs = getenv("PLEGMA_QCOSTS")
        if q_costs != "":
            var parts = q_costs.split(",")
            if len(parts) == 3:
                reference_costs = Costs.affine(Int(String(parts[0])), Int(String(parts[1])), Int(String(parts[2])))
            else:
                reference_costs = Costs.two_piece(
                    Int(String(parts[0])),
                    Int(String(parts[1])),
                    Int(String(parts[2])),
                    Int(String(parts[3])),
                    Int(String(parts[4])),
                )
        # PLEGMA_TREE=<share> (an experiment): along a tree from the reference through the sequences at least that
        # share of the input spells.
        var tree = getenv("PLEGMA_TREE")
        var bytes = on_tree(seqs, reference_costs, Float64(tree)) if tree != "" else on_reference(
            seqs, reference_costs, getenv("PLEGMA_SNAP") == "1", options.threads, getenv("PLEGMA_CONSOLIDATE") == "1"
        )
        if options.sp_set and options.sp_passes > 0:
            # With -I, refined on the sum of pairs with the reference's row weighing 4 times the rows' number: the
            # other rows move only among the placements that keep their alignment with the reference as good, the
            # one that suits the others best (in a tandem repeat, which copy a deletion takes).
            var weights = List[Int](length=len(seqs), fill=1)
            weights[0] = 4 * len(seqs)
            bytes = refine(
                bytes,
                options.sp_passes,
                options.sp_band,
                options.sp_max_band,
                # PLEGMA_REFINE_UNANCHORED=1 (an experiment): each row's band centred where it is, no partners aligned.
                getenv("PLEGMA_REFINE_UNANCHORED") != "1",
                # The bands' centres come from pairwise alignments at -A's vote costs (vntr-eval's affine): under the
                # reference's two-piece costs karan may refuse a divergent pair at its default memory.
                Costs.affine(4, 6, 2),
                options.sp_gap_open,
                False,
                # No group moves by default: on -A's MSAs they left the objective as it was (or within 0.03%) at a tenth
                # to a half of the refinement's time.
                options.group_band if options.group_set else 0,
                0,
                weights,
                1000,
                options.mismatch if options.prices_set else 2,
                options.gap if options.prices_set else 2,
            )
        var names = List[String]()
        for r in range(len(seqs)):
            names.append(records.names[r])
            anchored_rows.append(String(unsafe_from_utf8=bytes[r].copy()))
        graph = Graph.from_msa(names, anchored_rows)
        order.clear()
    elif options.anchored:
        var bytes = progressive(
            seqs,
            options.sp_passes if options.sp_set else 20,
            options.sp_band,
            options.sp_max_band,
            options.sp_gap_open,
            options.group_band if options.group_set else 256,
            options.units,
            True,
            options.one_indel if options.confidence_set else -2,
            options.votes,
            options.mismatch if options.prices_set else 2,
            options.gap if options.prices_set else 2,
            options.hard,
            options.nearest,
            options.consistency,
            options.subfamilies,
            options.shuffle,
            options.tree_start,
            options.exact_cells,
            options.interleave,
            options.corridor,
        )
        var names = List[String]()
        for r in range(len(seqs)):
            names.append(records.names[r])
            anchored_rows.append(String(unsafe_from_utf8=bytes[r].copy()))
        graph = Graph.from_msa(names, anchored_rows)
        order.clear()
    elif options.start_msa.byte_length() > 0:
        var start = read_fasta(options.start_msa)
        var names = List[String]()
        var rows = List[String]()
        for r in range(len(start)):
            names.append(start.names[r])
            var letters = List[UInt8](start.seq(r))
            for k in range(len(letters)):
                if letters[k] >= UInt8(ord("a")) and letters[k] <= UInt8(ord("z")):
                    letters[k] -= 32
            rows.append(String(unsafe_from_utf8=letters))
        graph = Graph.from_msa(names, rows)
    # Outputs list the starting MSA's rows, then every record in input order, whatever order they joined in.
    var path_of = List[Int](length=len(seqs), fill=-1)
    var shown = List[Int]()
    for p in range(len(graph.paths)):
        shown.append(p)
    var cells = 0
    var peak = 0
    var reversed_count = 0
    var row_of = List[Int]()
    var scratch = WindowScratch()
    var workspace = Workspace()
    for k in range(len(order)):
        var r = order[k]
        var name = records.names[r]
        if len(graph) == 0:
            path_of[r] = graph.add(name, Span(seqs[r]))
            continue
        var found = _align(graph, Span(seqs[r]), options, row_of, scratch, workspace)
        cells += found.cells
        peak = max(peak, found.peak_bytes)
        if options.strand:
            # abPOA's rule: a score under a third of the best a match run of the shorter length could earn.
            var floor = Int(options.scoring.match_score) * min(len(seqs[r]), len(graph.paths[0])) // 3
            if found.score < floor:
                var flipped = _reverse_complement(Span(seqs[r]))
                var other = _align(graph, Span(flipped), options, row_of, scratch, workspace)
                cells += other.cells
                if other.score > found.score:
                    seqs[r] = flipped^
                    found = other^
                    name += "_reverse_complement"
                    reversed_count += 1
        path_of[r] = graph.fuse_nodes(name, Span(seqs[r]), found.nodes)
    # Refinement: a fresh graph from the consensus so far, every sequence aligned to it again in the same order.
    for _ in range(options.rounds if options.start_msa.byte_length() == 0 and not options.anchored else 0):
        var members = List[Int]()
        for r in range(len(seqs)):
            members.append(path_of[r])
        var topology = Topology(graph)
        var backbone_letters = List[UInt8]()
        for node in graph.heaviest_path_of(topology.order, members):
            backbone_letters.append(graph.bases[node])
        if options.polish_backbone and options.mode == GLOBAL:
            var everyone = List[Int]()
            for r in range(len(seqs)):
                everyone.append(r)
            backbone_letters = polish(
                backbone_letters, seqs, everyone, options.polish_rounds, likely_rounds=options.likely_rounds
            )
        var next = Graph()
        _ = next.add("plegma_backbone", Span(backbone_letters))
        for k in range(len(order)):
            var r = order[k]
            var found = _align(next, Span(seqs[r]), options, row_of, scratch, workspace)
            cells += found.cells
            peak = max(peak, found.peak_bytes)
            path_of[r] = next.fuse_nodes(graph.names[path_of[r]], Span(seqs[r]), found.nodes)
        graph = next^
    var t2 = perf_counter_ns()
    if not options.anchored:
        for r in range(len(seqs)):
            shown.append(path_of[r])

    var result = options.result
    var text = String()
    var refining = options.sp_passes > 0 and result != 0 and not options.anchored
    var need_rows = result == 1 or result == 2 or options.consensus_rule == 1 or options.max_consensus > 1 or refining
    var by_path = graph.msa() if need_rows and not options.anchored else List[String]()
    var rows = List[String]()
    if need_rows and options.anchored:
        rows = anchored_rows.copy()
    elif need_rows:
        for p in shown:
            rows.append(by_path[p])
    if refining:
        # The MSA refined on its sum of pairs; a graph output is the refined MSA's graph.
        var bytes = List[List[UInt8]]()
        for row in rows:
            bytes.append(List[UInt8](row.as_bytes()))
        var refined = refine(
            bytes,
            options.sp_passes,
            options.sp_band,
            anchor_costs=Costs.two_piece(4, 4, 2, 24, 1),
            gap_open=options.sp_gap_open,
            seeded=options.sp_seeded,
            max_band=options.sp_max_band,
            group_band=options.group_band,
            mismatch=options.mismatch,
            gap=options.gap,
        )
        var names = List[String]()
        for k in range(len(rows)):
            rows[k] = String(unsafe_from_utf8=refined[k].copy())
            names.append(graph.names[shown[k]])
        if result == 3 or result == 4:
            graph = Graph.from_msa(names, rows)
            shown.clear()
            for k in range(len(rows)):
                shown.append(k)
    var consensus: String
    var consensus_path = List[Int]()
    var clusters = List[List[Int]]()
    var polishing = (options.polish_rounds > 0 or options.likely_rounds > 0) and options.mode == GLOBAL
    if result != 1 and result != 3:
        if options.max_consensus > 1:
            clusters = cluster_rows(rows, options.max_consensus, options.min_freq)
        var topology = Topology(graph)
        if len(clusters) <= 1:
            if options.consensus_rule == 1:
                consensus = _consensus_most_frequent(rows)
            else:
                # The sequences' paths alone weigh the consensus (not a refinement round's backbone).
                consensus_path = graph.heaviest_path_of(topology.order, shown)
                var letters = List[UInt8]()
                for node in consensus_path:
                    letters.append(graph.bases[node])
                consensus = String(unsafe_from_utf8=letters)
            if result == 0 or result == 2:
                if polishing:
                    var members = List[Int]()
                    for r in range(len(seqs)):
                        members.append(r)
                    consensus = _polished(consensus, seqs, members, options)
                text.write(">Consensus_sequence\n", consensus, "\n")
        elif result == 0 or result == 2:
            for k in range(len(clusters)):
                var cons: String
                if options.consensus_rule == 1:
                    var mine = List[String]()
                    for p in clusters[k]:
                        mine.append(rows[p])
                    cons = _consensus_most_frequent(mine)
                else:
                    var letters = List[UInt8]()
                    var paths = List[Int]()
                    for member in clusters[k]:
                        paths.append(shown[member])
                    for node in graph.heaviest_path_of(topology.order, paths):
                        letters.append(graph.bases[node])
                    cons = String(unsafe_from_utf8=letters)
                if polishing:
                    # A cluster's members past the starting MSA's rows are input records.
                    var members = List[Int]()
                    for member in clusters[k]:
                        if member >= len(shown) - len(seqs):
                            members.append(member - (len(shown) - len(seqs)))
                    cons = _polished(cons, seqs, members, options)
                var ids = String()
                for i in range(len(clusters[k])):
                    if i > 0:
                        ids += ","
                    ids += String(clusters[k][i])
                text.write(">Consensus_sequence_", k + 1, " ", ids, "\n", cons, "\n")
    if result == 1 or result == 2:
        for k in range(len(shown)):
            text.write(">", graph.names[shown[k]], "\n", rows[k], "\n")
    if result == 3 or result == 4:
        if result == 4 and len(consensus_path) == 0:
            # The graph's consensus path is the heaviest bundle, whichever consensus the FASTA output used.
            consensus_path = graph.heaviest_path_of(Topology(graph).order, shown)
        _write_gfa(graph, shown, consensus_path.copy() if result == 4 else List[Int](), text)
    if options.verbose:
        var t3 = perf_counter_ns()
        print(
            t"plegma: {input}: {len(records)} sequences ({reversed_count} reverse-complemented), {len(graph)} nodes,",
            t" {cells} cells, DP peak {peak // 1024} KB; read {Float64(t1 - t0) / 1e9} s,",
            t" align {Float64(t2 - t1) / 1e9} s, write {Float64(t3 - t2) / 1e9} s",
            file=stderr,
        )
    return text^


def main() raises:
    var raw = argv()
    var args = List[String]()
    for i in range(1, len(raw)):
        args.append(String(raw[i]))
    var options = Options()
    var output = String()
    var listed = False
    var input = String()
    var i = 0
    while i < len(args):
        var a = args[i]
        if a == "-h" or a == "--help":
            print(USAGE)
            return
        elif a == "-m":
            options.mode = Int(_value(args, i))
            i += 1
        elif a == "-M":
            options.scoring.match_score = Int32(Int(_value(args, i)))
            i += 1
        elif a == "-X":
            options.scoring.mismatch = Int32(Int(_value(args, i)))
            i += 1
        elif a == "-O":
            var o = _pair(_value(args, i))
            options.scoring.open1 = Int32(o[0])
            if o[1] >= 0:
                options.scoring.open2 = Int32(o[1])
            i += 1
        elif a == "-E":
            var e = _pair(_value(args, i))
            options.scoring.extend1 = Int32(e[0])
            if e[1] >= 0:
                options.scoring.extend2 = Int32(e[1])
            i += 1
        elif a == "-b":
            options.band_b = Int(_value(args, i))
            i += 1
        elif a == "-f":
            options.band_f = Float64(_value(args, i))
            i += 1
        elif a == "-L":
            options.longest_first = True
        elif a == "-p":
            options.progressive = True
        elif a == "-s":
            options.strand = True
        elif a == "-S":
            options.seeding = True
        elif a == "-J":
            options.split_seeds = True
        elif a == "-Y":
            options.agree_seeds = True
        elif a == "-Q":
            # The default with -A; kept for older command lines.
            pass
        elif a == "-N":
            options.voted = True
        elif a == "-t":
            options.threads = Int(_value(args, i))
            i += 1
        elif a == "-n":
            options.window = Int(_value(args, i))
            i += 1
        elif a == "-T":
            options.trim = Int(_value(args, i))
            i += 1
        elif a == "-K":
            var both = _pair(_value(args, i).replace(":", ","))
            options.repeat_k = both[0]
            options.repeat_distance = max(both[1], 0)
            i += 1
        elif a == "-H":
            options.lookback = Int(_value(args, i))
            i += 1
        elif a == "-l":
            listed = True
        elif a == "-i":
            options.start_msa = _value(args, i)
            i += 1
        elif a == "-o":
            output = _value(args, i)
            i += 1
        elif a == "-r":
            options.result = Int(_value(args, i))
            i += 1
        elif a == "-a":
            options.consensus_rule = Int(_value(args, i))
            i += 1
        elif a == "-R":
            options.rounds = Int(_value(args, i))
            i += 1
        elif a == "-I":
            var parts = _value(args, i).split(",")
            options.sp_passes = Int(parts[0])
            if len(parts) > 1:
                options.sp_band = Int(parts[1])
            if len(parts) > 2:
                options.sp_gap_open = Int(parts[2])
            if len(parts) > 3:
                options.sp_seeded = Int(parts[3]) != 0
            if len(parts) > 4:
                options.sp_max_band = Int(parts[4])
            options.sp_set = True
            i += 1
        elif a == "-G":
            options.group_band = Int(_value(args, i))
            options.group_set = True
            i += 1
        elif a == "-A":
            options.anchored = True
        elif a == "-U":
            options.units = False
        elif a == "-W":
            var both = _pair(_value(args, i))
            options.mismatch = both[0]
            options.gap = both[1] if both[1] >= 0 else both[0]
            options.prices_set = True
            i += 1
        elif a == "-V":
            var parts = _value(args, i).split(",")
            options.one_indel = Int(parts[0])
            if len(parts) > 1:
                options.votes = Int(parts[1])
            if len(parts) > 2:
                options.hard = Int(parts[2])
            if len(parts) > 3:
                options.nearest = Int(parts[3]) == 1
            if len(parts) > 4:
                options.consistency = Int(parts[4]) == 1
            if len(parts) > 5:
                options.shuffle = Int(parts[5])
            if len(parts) > 6:
                options.exact_cells = Int(parts[6]) * 1_000_000
            if len(parts) > 7:
                options.interleave = Int(parts[7])
            if len(parts) > 8:
                options.corridor = Int(parts[8])
            if len(parts) > 9:
                options.subfamilies = Int(parts[9]) == 1
            if len(parts) > 10:
                options.tree_start = Int(parts[10]) == 1
            options.confidence_set = True
            i += 1
        elif a == "-B":
            options.polish_backbone = True
        elif a == "-P":
            var both = _pair(_value(args, i))
            options.polish_rounds = both[0]
            if both[1] >= 0:
                options.likely_rounds = both[1]
            i += 1
        elif a == "-d":
            options.max_consensus = Int(_value(args, i))
            i += 1
        elif a == "-q":
            options.min_freq = Float64(_value(args, i))
            i += 1
        elif a == "-v":
            options.verbose = True
        elif a.startswith("-") and a.byte_length() > 1:
            raise Error(t"unknown option {a}\n{USAGE}")
        else:
            input = a
        i += 1
    if input.byte_length() == 0:
        print(USAGE)
        exit(1)
    if options.anchored and options.start_msa.byte_length() > 0:
        raise Error("-A builds its own MSA: it takes no -i")
    if options.mode != GLOBAL and options.mode != LOCAL and options.mode != EXTENSION:
        raise Error(t"-m {options.mode}: the modes are 0, 1 and 2")
    var inputs = List[String]()
    if listed:
        with open(input, "r") as f:
            # The lines are views of the text: keep it alive while they are read. Split on newlines alone:
            # this nightly's splitlines also splits at tabs, which a file name may hold.
            var content = f.read()
            for line in content.split("\n"):
                var name = String(line.strip())
                if name.byte_length() > 0:
                    inputs.append(name)
    else:
        inputs.append(input)
    var text = String()
    for path in inputs:
        text += _poa(path, options)
    if output.byte_length() > 0:
        with open(output, "w") as f:
            f.write(text)
    else:
        print(text, end="")
