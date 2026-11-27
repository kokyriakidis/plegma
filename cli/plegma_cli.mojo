"""The `plegma` command line: two subcommands over one FASTA/FASTQ file.

- `plegma poa`: a partial order alignment as abPOA runs one (its flags and output modes): each sequence, in abPOA's
  guide order, pinned to the graph's consensus so far, aligned to the graph between the pins and fused into it.
- `plegma ref`: the MSA anchored on the file's first sequence, the reference (`plegma.on_reference`): every other
  aligned to it alone, at minimap2's asm5 costs, pinned as minimap2 pins it given the reference's flanks (`-F`).

Both write a consensus (polished), the MSA or the graph as GFA (`-r`)."""

from std.atomic import Atomic
from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from karan import Costs
from karan.common import spread, thread_count
from sporos import read_fasta

from plegma.clusters import cluster_rows
from plegma.consensus import ConsensusScratch, polished
from plegma.graph import GAP, Graph
from plegma.homopolymer import STRAND_UNKNOWN
from plegma.kernel import BOTH
from plegma.poa import poa_with_scratch
from plegma.polish import polish
from plegma.reference import on_reference
from plegma.version import VERSION

from cli.gfa import write_gfa
from cli.options import USAGE, Options, parse
from cli.records import Records, upper


struct Aligned(Movable):
    """What either command made of a file, for the outputs: the graph, its paths in output order (`shown`: a starting
    MSA's rows, then every record in input order), the MSA's rows in that order when the command made them, the
    records a consensus may be polished against (`polish_with`) on their strands, the assembly's own consensus when
    reads are anchored at one end (`assembled`: one sequence, or the two sides' when they could not be joined), and
    the consensus `-r 0` writes when the library made it (`polished`)."""

    var graph: Graph
    var shown: List[Int]
    var rows: List[String]
    var polish_with: List[Int]
    var anchors: List[Int]
    var strands: List[Int]
    var assembled: List[List[UInt8]]
    var polished: List[List[UInt8]]

    def __init__(out self, var graph: Graph, var shown: List[Int], var rows: List[String]):
        self.graph = graph^
        self.shown = shown^
        self.rows = rows^
        self.polish_with = List[Int]()
        self.anchors = List[Int]()
        self.strands = List[Int]()
        self.assembled = List[List[UInt8]]()
        self.polished = List[List[UInt8]]()


def _poa(mut scratch: ConsensusScratch, mut records: Records, options: Options, run_threads: Int) raises -> Aligned:
    """`plegma poa`: `plegma.poa_with_scratch` with the command's options, each record anchored as its header says; the
    consensus `-r 0` writes as `plegma.consensus` makes it (`polished`), unless `-d` or `-a` asks for another."""
    if options.max_consensus > 1 or options.consensus_rule != 0:
        for anchor in records.anchors:
            if anchor != BOTH:
                raise Error(
                    "-d and -a take reads spanning both anchors: a read anchored at one end (anchor=left or right)"
                    " makes one assembly"
                )
    var start = Graph()
    if options.start_msa.byte_length() > 0:
        var msa = read_fasta(options.start_msa)
        var rows = List[String]()
        for r in range(len(msa)):
            rows.append(String(unsafe_from_utf8=upper(msa.seq(r))))
        start = Graph.from_msa(msa.names, rows)
    var cons_opts = options.consensus.copy()
    cons_opts.poa.threads = run_threads
    var made = poa_with_scratch(
        records.names, records.seqs, records.anchors, cons_opts.poa, start^, scratch.poa, scratch.aligner
    )
    var strands = records.strands.copy() if options.homopolymers.loaded() else List[Int]()
    var assembled = len(made.consensus) > 0
    var built_anchors = made.anchors.copy()
    var polished_parts = List[List[UInt8]]()
    if (options.result == 0 or options.result == 2) and (
        assembled or (options.max_consensus == 1 and options.consensus_rule == 0)
    ):
        polished_parts = polished(
            scratch, made, records.seqs, records.anchors, strands, cons_opts, options.homopolymers
        )
    # A read `-s` reverse-complemented is on the other strand to the consensus now (for -d's consensus sequences).
    for r in made.reversed:
        if r < len(strands) and strands[r] != STRAND_UNKNOWN:
            strands[r] = 1 - strands[r]
    var foreign = made.foreign.copy()
    var out = Aligned(Graph(), List[Int](), List[String]())
    var spanning = List[Int]()
    made^.parts(out.graph, out.shown, spanning, out.assembled, foreign)
    out.strands = strands^
    out.polished = polished_parts^
    if assembled:
        # Reads anchored at one end: every read but another copy's, each as it is anchored.
        var f = 0
        for r in range(len(records)):
            if f < len(foreign) and foreign[f] == r:
                f += 1
                continue
            out.polish_with.append(r)
        out.anchors = built_anchors^ if len(built_anchors) == len(records) else records.anchors.copy()
    else:
        out.polish_with = spanning^
    return out^


def _ref(records: Records, options: Options, run_threads: Int) raises -> Aligned:
    """`plegma ref`: `on_reference` at minimap2's asm5 scores as costs (match 1, mismatch 19, gap opens 39 and 81,
    extensions 3 and 1, made costs), the truth's convention (on chr20 its graphs genotype better than at pgrealign's
    scores: structural F1 0.709 against 0.658); reads past asm5's divergence, which it would delete and insert again,
    at pgrealign's (abPOA's 8, 16, 16+8k, 96+k halved: a mismatch 24, a gap of k letters min(16 + 12k, 96 + 5k)),
    whose long gaps cost more. The graph is the MSA's."""
    var flanks = List[List[UInt8]]()
    if options.flanks != "":
        var given = read_fasta(options.flanks)
        if len(given) != 2:
            raise Error("-F takes a FASTA of two sequences, the reference's left and right flanks")
        for r in range(2):
            flanks.append(upper(given.seq(r)))
    var bytes = on_reference(
        records.seqs,
        Costs.two_piece(40, 78, 7, 162, 3),
        run_threads,
        Costs.two_piece(24, 16, 12, 96, 5),
        flanks,
    )
    var rows = List[String]()
    var shown = List[Int]()
    for r in range(len(records)):
        rows.append(String(unsafe_from_utf8=bytes[r]))
        shown.append(r)
    var need_graph = (options.result == 0 and options.consensus_rule == 0) or options.result >= 2
    var graph: Graph
    if need_graph:
        graph = Graph.from_msa(records.names, rows)
    else:
        graph = Graph()
        graph.names = records.names.copy()
    var out = Aligned(graph^, shown^, rows^)
    for r in range(len(records)):
        out.polish_with.append(r)
    return out^


def _most_frequent(rows: List[String]) -> String:
    """The consensus of `-a 1`: each column's most frequent letter, the column skipped where a gap is."""
    var out = List[UInt8]()
    if len(rows) == 0:
        return String()
    var width = rows[0].byte_length()
    var counts = List[Int](length=256, fill=0)
    var touched = List[Int](capacity=8)
    for c in range(width):
        touched.clear()
        for r in rows:
            var b = Int(r.as_bytes()[c])
            if counts[b] == 0:
                touched.append(b)
            counts[b] += 1
        for a in range(1, len(touched)):
            var key = touched[a]
            var j = a
            while j > 0 and touched[j - 1] > key:
                touched[j] = touched[j - 1]
                j -= 1
            touched[j] = key
        var best = Int(GAP)
        for b in touched:
            if counts[b] > counts[best]:
                best = b
        for b in touched:
            counts[b] = 0
        if best != Int(GAP):
            out.append(UInt8(best))
    return String(unsafe_from_utf8=out)


def _consensus(
    aligned: Aligned,
    rows: List[String],
    members: List[Int],
    seqs: List[List[UInt8]],
    options: Options,
    run_threads: Int,
    order: List[Int],
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
        for node in aligned.graph.heaviest_path(order, paths):
            letters.append(aligned.graph.bases[node])
        consensus = String(unsafe_from_utf8=letters)
    if not options.polishing():
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
        List[UInt8](consensus.as_bytes()),
        seqs,
        records,
        options.consensus.polish_rounds,
        likely_rounds=options.consensus.likely_rounds,
        homopolymers=options.homopolymers,
        strands=aligned.strands,
        threads=run_threads,
    )
    return String(unsafe_from_utf8=polished)


def _write(var aligned: Aligned, seqs: List[List[UInt8]], options: Options, run_threads: Int) raises -> String:
    """The outputs `-r` asks for: consensus FASTA (one per cluster with `-d`), MSA FASTA, GFA."""
    var result = options.result
    var text = String()
    var need_rows = result == 1 or result == 2 or options.consensus_rule == 1 or options.max_consensus > 1
    var rows = List[String]()
    swap(rows, aligned.rows)
    if need_rows and len(rows) == 0:
        var by_path = aligned.graph.msa()
        for p in aligned.shown:
            rows.append(by_path[p])
    var everyone = List[Int]()
    for m in range(len(aligned.shown)):
        everyone.append(m)
    # A starting MSA's rows come first in `shown`; past them, record r is at first + r.
    var first = len(aligned.shown) - len(seqs)
    var graph_order = List[Int]()
    var need_order = (
        (result == 0 or result == 2) and len(aligned.polished) == 0 and options.consensus_rule == 0
    ) or result == 4
    if need_order and len(aligned.graph) > 0:
        var rank = aligned.graph.node_order()
        graph_order.resize(len(rank), 0)
        for node in range(len(rank)):
            graph_order[rank[node]] = node
    if (result == 0 or result == 2) and len(aligned.polished) > 0:
        # The library's consensus: one sequence, or the two sides of an assembly that could not be joined.
        if len(aligned.polished) == 1:
            text.write(">Consensus_sequence\n", StringSlice(unsafe_from_utf8=Span(aligned.polished[0])), "\n")
        else:
            text.write(">Consensus_sequence_left\n", StringSlice(unsafe_from_utf8=Span(aligned.polished[0])), "\n")
            text.write(">Consensus_sequence_right\n", StringSlice(unsafe_from_utf8=Span(aligned.polished[1])), "\n")
    elif result == 0 or result == 2:
        var clusters = cluster_rows(
            rows, options.max_consensus, options.min_freq, run_threads
        ) if options.max_consensus > 1 else (List[List[Int]]())
        if len(clusters) <= 1:
            text.write(
                ">Consensus_sequence\n",
                _consensus(aligned, rows, everyone, seqs, options, run_threads, graph_order),
                "\n",
            )
        else:
            var cons_list = List[String](length=len(clusters), fill=String())
            var workers = thread_count(run_threads, len(clusters))
            if workers > 1 and len(clusters) >= 2:
                var sub_threads = max(1, run_threads // workers)
                var failed = List[Bool](length=len(clusters), fill=False)
                var cons_p = cons_list.unsafe_ptr()
                var failed_p = failed.unsafe_ptr()
                var taken = Atomic[Int64](0)
                var n_clusters = len(clusters)

                def polish_cluster(
                    slot: Int,
                ) {
                    mut taken,
                    imm n_clusters,
                    imm aligned,
                    imm rows,
                    imm clusters,
                    imm seqs,
                    imm options,
                    imm sub_threads,
                    imm graph_order,
                    imm cons_p,
                    imm failed_p,
                }:
                    while True:
                        var k = Int(taken.fetch_add(1))
                        if k >= n_clusters:
                            return
                        try:
                            cons_p[unsafe_offset=k] = _consensus(
                                aligned, rows, clusters[k], seqs, options, sub_threads, graph_order
                            )
                        except:
                            failed_p[unsafe_offset=k] = True

                spread(polish_cluster, workers, workers)
                for k in range(len(clusters)):
                    if failed[k]:
                        cons_list[k] = _consensus(aligned, rows, clusters[k], seqs, options, run_threads, graph_order)
            else:
                for k in range(len(clusters)):
                    cons_list[k] = _consensus(aligned, rows, clusters[k], seqs, options, run_threads, graph_order)
            for k in range(len(clusters)):
                var ids = String()
                for i in range(len(clusters[k])):
                    if i > 0:
                        ids += ","
                    ids += String(clusters[k][i])
                text.write(">Consensus_sequence_", k + 1, " ", ids, "\n", cons_list[k], "\n")
    if result == 1 or result == 2:
        for k in range(len(aligned.shown)):
            text.write(">", aligned.graph.names[aligned.shown[k]], "\n", rows[k], "\n")
    if result == 3 or result == 4:
        var consensus_path = List[Int]()
        if result == 4 and len(aligned.assembled) == 0:
            # The graph's consensus path is the heaviest bundle, whichever consensus the FASTA output used.
            consensus_path = aligned.graph.heaviest_path(graph_order, aligned.shown)
        elif result == 4:
            # Reads anchored at one end: the bundle between the spanning reads' ends; with none, the assembled
            # sequence the reads were aligned to (its path the graph's first, not shown).
            var spanning = List[Int]()
            for r in aligned.polish_with:
                if len(aligned.anchors) == 0 or aligned.anchors[r] == BOTH:
                    spanning.append(aligned.shown[first + r])
            if len(spanning) > 0:
                consensus_path = aligned.graph.heaviest_path_spanning(graph_order, spanning)
            elif len(aligned.graph.paths) > 0:
                for node in aligned.graph.paths[0]:
                    consensus_path.append(Int(node))
        write_gfa(aligned.graph, aligned.shown, consensus_path, text)
    return text^


def _run(mut scratch: ConsensusScratch, input: String, options: Options, run_threads: Int) raises -> String:
    """One file through the command `options` names, as its outputs."""
    var t0 = perf_counter_ns()
    var records = Records(input, parse_headers=not options.reference)
    var t1 = perf_counter_ns()
    if len(records) == 0:
        return String(">Consensus_sequence\n\n") if options.result == 0 or options.result == 2 else String()
    var aligned = _ref(records, options, run_threads) if options.reference else _poa(
        scratch, records, options, run_threads
    )
    var nodes = len(aligned.graph)
    var t2 = perf_counter_ns()
    var text = _write(aligned^, records.seqs, options, run_threads)
    if options.verbose:
        var t3 = perf_counter_ns()
        var msg = String(
            "plegma: ",
            input,
            ": ",
            len(records),
            " sequences, ",
            nodes,
            " nodes; read ",
            Float64(t1 - t0) / 1e9,
            " s, align ",
            Float64(t2 - t1) / 1e9,
            " s, write ",
            Float64(t3 - t2) / 1e9,
            " s",
        )
        print(msg, file=stderr)
    return text^


def _inputs(options: Options) raises -> List[String]:
    """The files to align: the input, or with `-l` the files it lists, one per line."""
    var out = List[String]()
    if not options.listed:
        out.append(options.input)
        return out^
    with open(options.input, "r") as f:
        # The lines are views of the text: keep it alive while they are read. Split on newlines alone: this
        # nightly's splitlines also splits at tabs, which a file name may hold.
        var content = f.read()
        for line in content.split("\n"):
            var name = String(line.strip())
            if name.byte_length() > 0:
                out.append(name)
    return out^


def _main() raises:
    var raw = argv()
    var args = List[String]()
    for i in range(1, len(raw)):
        args.append(String(raw[i]))
    if len(args) == 0:
        print(USAGE, file=stderr)
        exit(1)
    for a in args:
        if a == "-h" or a == "--help":
            print(USAGE)
            return
    if args[0] == "--version":
        print("plegma", VERSION)
        return
    var options = parse(args)
    if options.input.byte_length() == 0:
        print(USAGE, file=stderr)
        exit(1)
    var paths = _inputs(options)
    var outputs = List[String](length=len(paths), fill=String())
    var workers = thread_count(options.threads, len(paths))
    if workers > 1 and len(paths) >= 2:
        var sub_threads = max(1, options.threads // workers)
        var failed = List[Bool](length=len(paths), fill=False)
        var outputs_p = outputs.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)
        var n_paths = len(paths)

        def run_file(
            slot: Int,
        ) {mut taken, imm n_paths, imm paths, imm options, imm sub_threads, imm outputs_p, imm failed_p}:
            var local_scratch = ConsensusScratch()
            while True:
                var idx = Int(taken.fetch_add(1))
                if idx >= n_paths:
                    return
                try:
                    outputs_p[unsafe_offset=idx] = _run(local_scratch, paths[idx], options, sub_threads)
                except:
                    failed_p[unsafe_offset=idx] = True

        spread(run_file, workers, workers)
        var scratch = ConsensusScratch()
        for idx in range(len(paths)):
            if failed[idx]:
                outputs[idx] = _run(scratch, paths[idx], options, options.threads)
    else:
        var scratch = ConsensusScratch()
        for idx in range(len(paths)):
            outputs[idx] = _run(scratch, paths[idx], options, options.threads)
    if options.output.byte_length() > 0:
        with open(options.output, "w") as f:
            for idx in range(len(outputs)):
                f.write(outputs[idx])
    else:
        for idx in range(len(outputs)):
            print(outputs[idx], end="")


def main():
    """Errors reach the user as one line on stderr and exit status 1."""
    try:
        _main()
    except e:
        print("plegma:", e, file=stderr)
        exit(1)
