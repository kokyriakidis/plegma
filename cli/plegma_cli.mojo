"""The `plegma` command line: a partial order alignment of the sequences of one FASTA/FASTQ file, as abPOA runs
one (its flags, defaults and output modes), written as a consensus, an MSA or a GFA graph."""

from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from sporos import read_fasta

from plegma.anchors import Anchor
from plegma.backbone import Level, MUMS, backbone, windows
from plegma.chain import ChainCosts
from plegma.graph import Graph, GAP
from plegma.clusters import cluster_rows
from plegma.order import guide_order
from plegma.kernel import EXTENSION, GLOBAL, LOCAL, GraphAlignment, Scoring, Topology, Workspace
from plegma.poa import SINK, SOURCE, WindowScratch, window_nodes
from plegma.polish import polish


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
    -P INT      rounds of consensus polishing, at most [50] (global mode): each sequence aligned to the consensus, and
                every one-letter edit kept that lowers their summed edit distance to it; 0 off
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
    mut row_of: List[Int],
    mut scratch: WindowScratch,
    mut workspace: Workspace,
) raises -> GraphAlignment:
    """`query` aligned to the graph through anchors on the consensus path: each window between two anchors
    aligned to the subgraph between their nodes, the anchors' letters matched."""
    var path = graph.heaviest_path(topology.order)
    var letters = List[UInt8](capacity=len(path))
    for node in path:
        letters.append(graph.bases[node])
    # MUMs, placed by uniqueness, of 51 letters, then 31, then 21, each gap seeded again until it yields nothing;
    # none inside a tandem repeat of the path (its k-mers recurring nearby), where even a unique exact match may
    # pair different copies: there the DP aligns the whole array.
    var levels = List[Level]()
    for n in [51, 31, 21]:
        levels.append(Level(MUMS, n))
    var chained = backbone(
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
    # Each anchor loses `trim` letters at both ends: the DP, not the anchor, decides where a gap next to it goes.
    var anchors = List[Anchor]()
    for a in chained:
        if a.length > 2 * trim:
            anchors.append(
                Anchor(a.ref_pos + trim, a.query_pos + trim, a.length - 2 * trim, a.ref_copies, a.query_copies)
            )
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


def _polished(consensus: String, seqs: List[List[UInt8]], members: List[Int], rounds: Int) raises -> String:
    if len(members) == 0:
        return consensus
    var letters = List[UInt8](consensus.as_bytes())
    return String(unsafe_from_utf8=polish(letters, seqs, members, rounds))


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
    var graph = Graph()
    if options.start_msa.byte_length() > 0:
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
    var t1 = perf_counter_ns()
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
    for _ in range(options.rounds if options.start_msa.byte_length() == 0 else 0):
        var members = List[Int]()
        for r in range(len(seqs)):
            members.append(path_of[r])
        var topology = Topology(graph)
        var backbone_letters = List[UInt8]()
        for node in graph.heaviest_path_of(topology.order, members):
            backbone_letters.append(graph.bases[node])
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
    for r in range(len(seqs)):
        shown.append(path_of[r])

    var result = options.result
    var text = String()
    var need_rows = result == 1 or result == 2 or options.consensus_rule == 1 or options.max_consensus > 1
    var by_path = graph.msa() if need_rows else List[String]()
    var rows = List[String]()
    if need_rows:
        for p in shown:
            rows.append(by_path[p])
    var consensus: String
    var consensus_path = List[Int]()
    var clusters = List[List[Int]]()
    var polishing = options.polish_rounds > 0 and options.mode == GLOBAL
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
                    consensus = _polished(consensus, seqs, members, options.polish_rounds)
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
                    cons = _polished(cons, seqs, members, options.polish_rounds)
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
        elif a == "-P":
            options.polish_rounds = Int(_value(args, i))
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
