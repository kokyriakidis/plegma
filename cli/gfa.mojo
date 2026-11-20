"""`plegma poa -r 3` and `-r 4`: the graph as GFA 1.0."""

from plegma.graph import Graph


def write_gfa(graph: Graph, shown: List[Int], consensus_path: List[Int], mut text: String) raises:
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
