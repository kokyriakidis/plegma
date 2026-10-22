"""The `plegma` command line: a partial order alignment of the sequences of one FASTA/FASTQ file, as abPOA runs
one (its flags, defaults and output modes), written as a consensus, an MSA or a GFA graph."""

from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from sporos import read_fasta

from plegma.graph import Graph, GAP
from plegma.kernel import EXTENSION, GLOBAL, LOCAL, Scoring, Topology, align_to_graph


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
  Order:
    -L          align the sequences longest first
  Input/Output:
    -i FILE     start from the graph of an existing MSA (FASTA, rows with '-')
    -o FILE     output to FILE [stdout]
    -r INT      output [0]: 0 consensus FASTA, 1 MSA FASTA, 2 both, 3 GFA, 4 GFA with the consensus path
    -a INT      consensus [0]: 0 heaviest path, 1 most frequent letter per column
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


def main() raises:
    var raw = argv()
    var args = List[String]()
    for i in range(1, len(raw)):
        args.append(String(raw[i]))
    var mode = GLOBAL
    var scoring = Scoring()
    var band_b = 10
    var band_f = 0.01
    var longest_first = False
    var start_msa = String()
    var output = String()
    var result = 0
    var consensus_rule = 0
    var verbose = False
    var input = String()
    var i = 0
    while i < len(args):
        var a = args[i]
        if a == "-h" or a == "--help":
            print(USAGE)
            return
        elif a == "-m":
            mode = Int(_value(args, i))
            i += 1
        elif a == "-M":
            scoring.match_score = Int32(Int(_value(args, i)))
            i += 1
        elif a == "-X":
            scoring.mismatch = Int32(Int(_value(args, i)))
            i += 1
        elif a == "-O":
            var o = _pair(_value(args, i))
            scoring.open1 = Int32(o[0])
            if o[1] >= 0:
                scoring.open2 = Int32(o[1])
            i += 1
        elif a == "-E":
            var e = _pair(_value(args, i))
            scoring.extend1 = Int32(e[0])
            if e[1] >= 0:
                scoring.extend2 = Int32(e[1])
            i += 1
        elif a == "-b":
            band_b = Int(_value(args, i))
            i += 1
        elif a == "-f":
            band_f = Float64(_value(args, i))
            i += 1
        elif a == "-L":
            longest_first = True
        elif a == "-i":
            start_msa = _value(args, i)
            i += 1
        elif a == "-o":
            output = _value(args, i)
            i += 1
        elif a == "-r":
            result = Int(_value(args, i))
            i += 1
        elif a == "-a":
            consensus_rule = Int(_value(args, i))
            i += 1
        elif a == "-v":
            verbose = True
        elif a.startswith("-") and a.byte_length() > 1:
            raise Error(t"unknown option {a}\n{USAGE}")
        else:
            input = a
        i += 1
    if input.byte_length() == 0:
        print(USAGE)
        exit(1)
    if mode != GLOBAL and mode != LOCAL and mode != EXTENSION:
        raise Error(t"-m {mode}: the modes are 0, 1 and 2")

    var t0 = perf_counter_ns()
    var records = read_fasta(input)
    var order = List[Int]()
    for r in range(len(records)):
        order.append(r)
    if longest_first:
        # Longest first, file order among equals: sort keys (2^32 - length) << 32 | index.
        var keys = List[UInt64]()
        for r in range(len(records)):
            keys.append((UInt64((1 << 32) - len(records.seq(r))) << 32) | UInt64(r))
        sort(keys)
        for k in range(len(keys)):
            order[k] = Int(keys[k] & 0xFFFFFFFF)
    var graph = Graph()
    if start_msa.byte_length() > 0:
        var start = read_fasta(start_msa)
        var names = List[String]()
        var rows = List[String]()
        for r in range(len(start)):
            names.append(start.names[r])
            rows.append(String(unsafe_from_utf8=List[UInt8](start.seq(r))))
        graph = Graph.from_msa(names, rows)
    var t1 = perf_counter_ns()
    var cells = 0
    var peak = 0
    for k in range(len(order)):
        var r = order[k]
        var name = records.names[r]
        if len(graph) == 0:
            _ = graph.add(name, records.seq(r))
            continue
        var topology = Topology(graph)
        var found = align_to_graph(graph, topology, records.seq(r), scoring, mode, band_b, band_f)
        cells += found.cells
        peak = max(peak, found.peak_bytes)
        _ = graph.fuse_nodes(name, records.seq(r), found.nodes)
    var t2 = perf_counter_ns()

    var text = String()
    var need_rows = result == 1 or result == 2 or consensus_rule == 1
    var rows = graph.msa() if need_rows else List[String]()
    var consensus = String()
    var consensus_path = List[Int]()
    if result != 1 and result != 3:
        if consensus_rule == 1:
            consensus = _consensus_most_frequent(rows)
        else:
            var topology = Topology(graph)
            consensus_path = graph.heaviest_path(topology.order)
            var letters = List[UInt8]()
            for node in consensus_path:
                letters.append(graph.bases[node])
            consensus = String(unsafe_from_utf8=letters)
    if result == 0 or result == 2:
        text += ">Consensus_sequence\n" + consensus + "\n"
    if result == 1 or result == 2:
        for p in range(len(graph.paths)):
            text += ">" + graph.names[p] + "\n" + rows[p] + "\n"
    if result == 3 or result == 4:
        text += "H\tVN:Z:1.0\n"
        for node in range(len(graph)):
            text.write("S\t", node + 1, "\t", chr(Int(graph.bases[node])), "\n")
        for node in range(len(graph)):
            for s in graph.successors[node]:
                text.write("L\t", node + 1, "\t+\t", s + 1, "\t+\t0M\n")
        for p in range(len(graph.paths)):
            text += "P\t" + graph.names[p] + "\t"
            for k in range(len(graph.paths[p])):
                if k > 0:
                    text += ","
                text += String(graph.paths[p][k] + 1) + "+"
            text += "\t*\n"
        if result == 4:
            text += "P\tConsensus_sequence\t"
            for k in range(len(consensus_path)):
                if k > 0:
                    text += ","
                text += String(consensus_path[k] + 1) + "+"
            text += "\t*\n"
    if output.byte_length() > 0:
        with open(output, "w") as f:
            f.write(text)
    else:
        print(text, end="")
    if verbose:
        var t3 = perf_counter_ns()
        print(
            t"plegma: {len(records)} sequences, {len(graph)} nodes, {cells} cells, DP peak {peak // 1024} KB;",
            t" read {Float64(t1 - t0) / 1e9} s,",
            t" align {Float64(t2 - t1) / 1e9} s, write {Float64(t3 - t2) / 1e9} s",
            file=stderr,
        )
