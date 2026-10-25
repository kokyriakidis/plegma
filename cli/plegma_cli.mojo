"""The `plegma` command line: two subcommands over one FASTA/FASTQ file.

- `plegma poa`: a partial order alignment as abPOA runs one (its flags and output modes): each sequence, in abPOA's
  guide order, pinned to the graph's consensus so far, aligned to the graph between the pins and fused into it.
- `plegma ref`: the MSA anchored on the file's first sequence, the reference (`plegma.on_reference`): every other
  aligned to it alone, at minimap2's asm5 costs, pinned as minimap2 pins it given the reference's flanks (`-F`).

Both write a consensus (polished), the MSA or the graph as GFA (`-r`)."""

from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from karan import Costs
from sporos import read_fasta

from plegma.clusters import cluster_rows
from plegma.graph import GAP, Graph
from plegma.kernel import EXTENSION, GLOBAL, LOCAL, GraphAlignment, Scoring, SINK, SOURCE, Topology, Workspace
from plegma.order import guide_order
from plegma.pinned import WindowScratch, align_pinned
from plegma.polish import polish
from plegma.reference import on_reference


comptime USAGE = """Usage: plegma <command> [options] <in.fa/fq> > cons.fa/msa.fa/plegma.gfa

Commands:
  poa           partial order alignment of the sequences, as abPOA runs one
  ref           the MSA anchored on the first sequence (the reference): every other aligned to it alone at
                minimap2's asm5 scores (as costs: mismatch 40, a gap of k letters min(78 + 7k, 162 + 3k)), gaps
                left-normalized, its letters in the reference's columns and the runs inserted at one place aligned to
                each other on the longest; reads of one sequence past asm5's divergence fall back to pgrealign's scores

poa options:
    -m INT      alignment mode [0]: 0 global, 1 local, 2 extension
    -M INT      match score [2]
    -X INT      mismatch penalty [4]
    -O INT[,INT] gap opening penalties O1,O2 [4,24]
    -E INT[,INT] gap extension penalties E1,E2 [2,1]
                a g-long gap costs min{O1+g*E1, O2+g*E2} (convex); O2 0: O1+g*E1 (affine); O1 0: g*E1 (linear)
    -b INT      adaptive band: extra band b [10]; < 0 disables the band
    -f FLOAT    adaptive band: extra band f [0.01]: b+f*L letters either side, L the sequence's length
    -p          progressive: abPOA's guide order, most similar minimizer multisets first [the default]
    -L          the sequences longest first instead
    -n          the sequences in input order instead
    -S          seeded (global mode): each sequence pinned to the consensus path as minimap2 pins an alignment
                (asm5's seeds and chain; pins only on 19-mers unique in both), the DP only between the pins
                [the default]
    -N          not seeded: each sequence aligned to the whole graph
    -s          ambiguous strand: a sequence scoring poorly is tried reverse-complemented, the better kept
    -i FILE     start from the graph of an existing MSA (FASTA, rows with '-')

ref options:
    -F FILE     the reference's flanks (FASTA: left, then right): each alignment to the reference pinned as minimap2
                pins it (-x asm5 as dipcall runs it; its seeds and chain over flank + sequence + flank), exact between
                the pins, so differences are written as dipcall's truth sets write them
    -t INT      threads for the alignments to the reference [1]

Output options (both):
    -r INT      output [0]: 0 consensus FASTA, 1 MSA FASTA, 2 both, 3 GFA, 4 GFA with the consensus path
    -a INT      consensus [0]: 0 heaviest path, 1 most frequent letter per column
    -d INT      at most INT consensus sequences [1]; at 2 or more the sequences are clustered on the MSA's
                heterozygous columns and each cluster gets its own
    -q FLOAT    smallest share of the sequences a cluster and a column's second allele must hold [0.25]
    -P INT[,INT] rounds of consensus polishing, at most [50,20] (global alignment): each sequence aligned to the
                consensus, and every one-letter or whole-repeat-unit edit kept that lowers their summed edit distance
                to it (the first), then that raises their likelihood under a pair HMM (the second); 0,0 off
    -l          the input is a list of sequence files, one per line, each aligned on its own
    -o FILE     output to FILE [stdout]
    -v          report timings on stderr
    -h          this help
"""


struct Options(Copyable, Movable):
    var reference: Bool
    """`ref`, else `poa`."""
    var mode: Int
    var scoring: Scoring
    var band_b: Int
    var band_f: Float64
    var longest_first: Bool
    var progressive: Bool
    var strand: Bool
    var seeded: Bool
    var start_msa: String
    var flanks: String
    var threads: Int
    var result: Int
    var consensus_rule: Int
    var max_consensus: Int
    var min_freq: Float64
    var polish_rounds: Int
    var likely_rounds: Int
    var verbose: Bool

    def __init__(out self, reference: Bool):
        self.reference = reference
        self.mode = GLOBAL
        self.scoring = Scoring()
        self.band_b = 10
        self.band_f = 0.01
        self.longest_first = False
        # The guide order and the pins by default: on the read sets the same consensus as the input order unpinned,
        # on chr20's allele sets better MSAs, on its largest VNTRs several times faster.
        self.progressive = True
        self.strand = False
        self.seeded = True
        self.start_msa = String()
        self.flanks = String()
        self.threads = 1
        self.result = 0
        self.consensus_rule = 0
        self.max_consensus = 1
        self.min_freq = 0.25
        self.polish_rounds = 50
        self.likely_rounds = 20
        self.verbose = False


struct Aligned(Movable):
    """What either command made of a file, for the outputs: the graph, its paths in output order (`shown`: a starting
    MSA's rows, then every record in input order), and the MSA's rows in that order when the command made them."""

    var graph: Graph
    var shown: List[Int]
    var rows: List[String]

    def __init__(out self, var graph: Graph, var shown: List[Int], var rows: List[String]):
        self.graph = graph^
        self.shown = shown^
        self.rows = rows^


def _value(args: List[String], i: Int) raises -> String:
    if i + 1 >= len(args):
        raise Error(t"{args[i]} takes a value")
    return args[i + 1]


def _pair(text: String) raises -> Tuple[Int, Int]:
    var parts = text.split(",")
    if len(parts) == 1:
        return (Int(parts[0]), -1)
    return (Int(parts[0]), Int(parts[1]))


def _upper(seq: Span[UInt8, _]) -> List[UInt8]:
    """Upper case, as abPOA reads its input: a and A are one letter to the kernel and to the graph alike."""
    var out = List[UInt8](seq)
    for k in range(len(out)):
        if out[k] >= UInt8(ord("a")) and out[k] <= UInt8(ord("z")):
            out[k] -= 32
    return out^


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


def _poa(mut names: List[String], mut seqs: List[List[UInt8]], options: Options) raises -> Aligned:
    """`plegma poa`: each sequence aligned to the whole graph so far (in input order, `-L` longest first, or `-p`'s
    guide order) and fused into it; with `-s`, one scoring poorly tried reverse-complemented (it and its name then
    changed in `seqs` and `names`)."""
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
        var rows = List[String]()
        for r in range(len(start)):
            rows.append(String(unsafe_from_utf8=_upper(start.seq(r))))
        graph = Graph.from_msa(start.names, rows)
    var shown = List[Int]()
    for p in range(len(graph.paths)):
        shown.append(p)
    var path_of = List[Int](length=len(seqs), fill=-1)
    var row_of = List[Int]()
    var scratch = WindowScratch()
    var workspace = Workspace()
    for r in order:
        if len(graph) == 0:
            path_of[r] = graph.add(names[r], Span(seqs[r]))
            continue
        var found = _align(graph, seqs[r], options, row_of, scratch, workspace)
        if options.strand:
            # abPOA's rule: a score under a third of the best a match run of the shorter length could earn.
            var floor = Int(options.scoring.match_score) * min(len(seqs[r]), len(graph.paths[0])) // 3
            if found.score < floor:
                var flipped = _reverse_complement(Span(seqs[r]))
                var other = _align(graph, flipped, options, row_of, scratch, workspace)
                if other.score > found.score:
                    seqs[r] = flipped^
                    found = other^
                    names[r] += "_reverse_complement"
        path_of[r] = graph.fuse_nodes(names[r], Span(seqs[r]), found.nodes)
    for r in range(len(seqs)):
        shown.append(path_of[r])
    return Aligned(graph^, shown^, List[String]())


def _align(
    graph: Graph,
    seq: List[UInt8],
    options: Options,
    mut row_of: List[Int],
    mut scratch: WindowScratch,
    mut workspace: Workspace,
) raises -> GraphAlignment:
    var topology = Topology(graph)
    row_of.resize(len(graph), -1)
    if options.seeded and options.mode == GLOBAL:
        return align_pinned(
            graph, topology, seq, options.scoring, options.band_b, options.band_f, row_of, scratch, workspace
        )
    return workspace.align(
        graph,
        topology,
        topology.order,
        SOURCE,
        SINK,
        Span(seq),
        options.scoring,
        row_of,
        options.mode,
        options.band_b,
        options.band_f,
    )


def _ref(names: List[String], seqs: List[List[UInt8]], options: Options) raises -> Aligned:
    """`plegma ref`: `on_reference` at minimap2's asm5 scores as costs (match 1, mismatch 19, gap opens 39 and 81,
    extensions 3 and 1, made costs), the truth's convention (on chr20 its graphs genotype better than at pgrealign's
    scores: structural F1 0.709 against 0.658); reads past asm5's divergence, which it would delete and insert again,
    at pgrealign's (abPOA's 8, 16, 16+8k, 96+k halved: a mismatch 24, a gap of k letters min(16 + 12k, 96 + 5k)),
    whose long gaps cost more. The graph is the MSA's."""
    var flanks = List[List[UInt8]]()
    if options.flanks != "":
        var records = read_fasta(options.flanks)
        if len(records) != 2:
            raise Error("-F takes a FASTA of two sequences, the reference's left and right flanks")
        for r in range(2):
            flanks.append(_upper(records.seq(r)))
    var bytes = on_reference(
        seqs,
        Costs.two_piece(40, 78, 7, 162, 3),
        options.threads,
        Costs.two_piece(24, 16, 12, 96, 5),
        flanks,
    )
    var rows = List[String]()
    var shown = List[Int]()
    for r in range(len(seqs)):
        rows.append(String(unsafe_from_utf8=bytes[r]))
        shown.append(r)
    return Aligned(Graph.from_msa(names, rows), shown^, rows^)


def _most_frequent(rows: List[String]) -> String:
    """The consensus of `-a 1`: each column's most frequent letter, the column skipped where a gap is."""
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


def _consensus(
    aligned: Aligned, rows: List[String], members: List[Int], seqs: List[List[UInt8]], options: Options
) raises -> String:
    """The consensus of the paths `members` (indices into `shown` and `rows`), polished against their sequences."""
    var consensus: String
    if options.consensus_rule == 1:
        var mine = List[String]()
        for m in members:
            mine.append(rows[m])
        consensus = _most_frequent(mine)
    else:
        var paths = List[Int]()
        for m in members:
            paths.append(aligned.shown[m])
        var letters = List[UInt8]()
        for node in aligned.graph.heaviest_path(Topology(aligned.graph).order, paths):
            letters.append(aligned.graph.bases[node])
        consensus = String(unsafe_from_utf8=letters)
    if (options.polish_rounds == 0 and options.likely_rounds == 0) or options.mode != GLOBAL:
        return consensus
    # A starting MSA's rows come first in `shown`; past them, member m is input record m - first.
    var first = len(aligned.shown) - len(seqs)
    var records = List[Int]()
    for m in members:
        if m >= first:
            records.append(m - first)
    if len(records) == 0:
        return consensus
    var polished = polish(
        List[UInt8](consensus.as_bytes()), seqs, records, options.polish_rounds, likely_rounds=options.likely_rounds
    )
    return String(unsafe_from_utf8=polished)


def _write(var aligned: Aligned, seqs: List[List[UInt8]], options: Options) raises -> String:
    """The outputs `-r` asks for: consensus FASTA (one per cluster with `-d`), MSA FASTA, GFA."""
    var result = options.result
    var text = String()
    var need_rows = result == 1 or result == 2 or options.consensus_rule == 1 or options.max_consensus > 1
    var rows = aligned.rows.copy()
    if need_rows and len(rows) == 0:
        var by_path = aligned.graph.msa()
        for p in aligned.shown:
            rows.append(by_path[p])
    var everyone = List[Int]()
    for m in range(len(aligned.shown)):
        everyone.append(m)
    if result == 0 or result == 2:
        var clusters = cluster_rows(rows, options.max_consensus, options.min_freq) if options.max_consensus > 1 else (
            List[List[Int]]()
        )
        if len(clusters) <= 1:
            text.write(">Consensus_sequence\n", _consensus(aligned, rows, everyone, seqs, options), "\n")
        else:
            for k in range(len(clusters)):
                var ids = String()
                for i in range(len(clusters[k])):
                    if i > 0:
                        ids += ","
                    ids += String(clusters[k][i])
                text.write(
                    ">Consensus_sequence_",
                    k + 1,
                    " ",
                    ids,
                    "\n",
                    _consensus(aligned, rows, clusters[k], seqs, options),
                    "\n",
                )
    if result == 1 or result == 2:
        for k in range(len(aligned.shown)):
            text.write(">", aligned.graph.names[aligned.shown[k]], "\n", rows[k], "\n")
    if result == 3 or result == 4:
        var consensus_path = List[Int]()
        if result == 4:
            # The graph's consensus path is the heaviest bundle, whichever consensus the FASTA output used.
            consensus_path = aligned.graph.heaviest_path(Topology(aligned.graph).order, aligned.shown)
        _write_gfa(aligned.graph, aligned.shown, consensus_path, text)
    return text^


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


def _run(input: String, options: Options) raises -> String:
    """One file through the command `options` names, as its outputs."""
    var t0 = perf_counter_ns()
    var records = read_fasta(input)
    var names = records.names.copy()
    var seqs = List[List[UInt8]]()
    for r in range(len(records)):
        seqs.append(_upper(records.seq(r)))
    var t1 = perf_counter_ns()
    var aligned = _ref(names, seqs, options) if options.reference else _poa(names, seqs, options)
    var nodes = len(aligned.graph)
    var t2 = perf_counter_ns()
    var text = _write(aligned^, seqs, options)
    if options.verbose:
        var t3 = perf_counter_ns()
        print(
            t"plegma: {input}: {len(records)} sequences, {nodes} nodes; read {Float64(t1 - t0) / 1e9} s,",
            t" align {Float64(t2 - t1) / 1e9} s, write {Float64(t3 - t2) / 1e9} s",
            file=stderr,
        )
    return text^


def main() raises:
    var raw = argv()
    var args = List[String]()
    for i in range(1, len(raw)):
        args.append(String(raw[i]))
    if len(args) == 0 or args[0] == "-h" or args[0] == "--help":
        print(USAGE)
        exit(0 if len(args) > 0 else 1)
    if args[0] != "poa" and args[0] != "ref":
        raise Error(t"unknown command {args[0]}: poa or ref\n{USAGE}")
    var reference = args[0] == "ref"
    var options = Options(reference)
    var output = String()
    var listed = False
    var input = String()
    var i = 1
    while i < len(args):
        var a = args[i]
        if a == "-h" or a == "--help":
            print(USAGE)
            return
        # poa's own.
        elif a == "-m" and not reference:
            options.mode = Int(_value(args, i))
            i += 1
        elif a == "-M" and not reference:
            options.scoring.match_score = Int32(Int(_value(args, i)))
            i += 1
        elif a == "-X" and not reference:
            options.scoring.mismatch = Int32(Int(_value(args, i)))
            i += 1
        elif a == "-O" and not reference:
            var o = _pair(_value(args, i))
            options.scoring.open1 = Int32(o[0])
            if o[1] >= 0:
                options.scoring.open2 = Int32(o[1])
            i += 1
        elif a == "-E" and not reference:
            var e = _pair(_value(args, i))
            options.scoring.extend1 = Int32(e[0])
            if e[1] >= 0:
                options.scoring.extend2 = Int32(e[1])
            i += 1
        elif a == "-b" and not reference:
            options.band_b = Int(_value(args, i))
            i += 1
        elif a == "-f" and not reference:
            options.band_f = Float64(_value(args, i))
            i += 1
        elif a == "-L" and not reference:
            options.longest_first = True
            options.progressive = False
        elif a == "-n" and not reference:
            options.longest_first = False
            options.progressive = False
        elif a == "-p" and not reference:
            options.progressive = True
            options.longest_first = False
        elif a == "-s" and not reference:
            options.strand = True
        elif a == "-S" and not reference:
            options.seeded = True
        elif a == "-N" and not reference:
            options.seeded = False
        elif a == "-i" and not reference:
            options.start_msa = _value(args, i)
            i += 1
        # ref's own.
        elif a == "-F" and reference:
            options.flanks = _value(args, i)
            i += 1
        elif a == "-t" and reference:
            options.threads = Int(_value(args, i))
            i += 1
        # Both.
        elif a == "-l":
            listed = True
        elif a == "-o":
            output = _value(args, i)
            i += 1
        elif a == "-r":
            options.result = Int(_value(args, i))
            i += 1
        elif a == "-a":
            options.consensus_rule = Int(_value(args, i))
            i += 1
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
            raise Error(t"{args[0]}: unknown option {a}\n{USAGE}")
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
        text += _run(path, options)
    if output.byte_length() > 0:
        with open(output, "w") as f:
            f.write(text)
    else:
        print(text, end="")
