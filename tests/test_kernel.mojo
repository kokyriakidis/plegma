from std.testing import assert_equal, assert_true, TestSuite

from karan import Costs, distance

from plegma.graph import Graph
from plegma.kernel import (
    EXTENSION,
    GLOBAL,
    LOCAL,
    PREFIX,
    SINK,
    SOURCE,
    SUFFIX,
    Scoring,
    Topology,
    Workspace,
    align_to_graph,
)


def _random(n: Int, seed: Int) -> String:
    var out = String()
    var state = UInt64(seed) * 0x9E3779B97F4A7C15 + 1
    for _ in range(n):
        state ^= state >> 12
        state ^= state << 25
        state ^= state >> 27
        out += "ACGT"[byte=Int((state * 0x2545F4914F6CDD1D) >> 62)]
    return out^


def _mutated(src: String, every: Int, seed: Int) -> String:
    var out = String()
    var state = UInt64(seed) * 0x9E3779B97F4A7C15 + 3
    for c in src.codepoint_slices():
        state ^= state >> 12
        state ^= state << 25
        state ^= state >> 27
        var r = state * 0x2545F4914F6CDD1D
        if Int(r % UInt64(every)) == 0:
            var kind = (r >> 40) % 3
            if kind == 0:
                out += "ACGT"[byte=Int(r >> 62)]
            elif kind == 2:
                out += c
                out += "ACGT"[byte=Int(r >> 62)]
            continue
        out += c
    return out^


def _score_from_cost(a: String, b: String, s: Scoring) raises -> Int:
    """abPOA's global score of a pair, from karan's least cost: with a match earning M, the score is
    M (n + m) / 2 less the cost of mismatches at X + M and gap letters at E + M / 2."""
    var half = Int(s.match_score) // 2
    var costs = Costs.two_piece(
        Int(s.mismatch) + Int(s.match_score),
        Int(s.open1),
        Int(s.extend1) + half,
        Int(s.open2),
        Int(s.extend2) + half,
    )
    return half * (a.byte_length() + b.byte_length()) - distance(a, b, costs)


def test_one_path_matches_karan() raises:
    var s = Scoring()
    for case_id in range(12):
        var a = _random(40 + 37 * case_id, case_id)
        var b = _mutated(a, 6 + case_id % 5, case_id + 99)
        var graph = Graph()
        _ = graph.add("a", a.as_bytes())
        var topology = Topology(graph)
        var want = _score_from_cost(a, b, s)
        var full = align_to_graph(graph, topology, b.as_bytes(), s, GLOBAL, -1)
        assert_equal(full.score, want)
        var banded = align_to_graph(graph, topology, b.as_bytes(), s, GLOBAL, 10, 0.01)
        assert_true(banded.score <= want)
        if case_id < 8:
            assert_equal(banded.score, want)
        # The traceback lands every letter on a distinct node, in order.
        var last = -1
        for node in full.nodes:
            if node >= 0:
                assert_true(topology.rank[node] > last)
                last = topology.rank[node]


def test_more_paths_never_score_less() raises:
    var s = Scoring()
    var a = _random(300, 5)
    var b = _mutated(a, 8, 6)
    var c = _mutated(a, 8, 7)
    var graph = Graph()
    _ = graph.add("a", a.as_bytes())
    var topology = Topology(graph)
    var found = align_to_graph(graph, topology, b.as_bytes(), s, GLOBAL, -1)
    _ = graph.fuse_nodes("b", b.as_bytes(), found.nodes)
    var both = Topology(graph)
    var to_both = align_to_graph(graph, both, c.as_bytes(), s, GLOBAL, -1)
    assert_true(to_both.score >= _score_from_cost(a, c, s))
    assert_true(to_both.score >= _score_from_cost(b, c, s))
    # A sequence already in the graph aligns at full match score.
    assert_equal(align_to_graph(graph, both, b.as_bytes(), s, GLOBAL, -1).score, 2 * b.byte_length())


def test_local_and_extension() raises:
    var s = Scoring()
    var a = _random(400, 11)
    var graph = Graph()
    _ = graph.add("a", a.as_bytes())
    var topology = Topology(graph)
    var inner = String(a[byte=100:220])
    var flanked = _random(30, 12) + inner + _random(25, 13)
    var local = align_to_graph(graph, topology, flanked.as_bytes(), s, LOCAL)
    assert_true(local.score >= 2 * 120)
    for j in range(30, 150):
        assert_equal(local.nodes[j], Int(graph.paths[0][100 + j - 30]))
    var prefix = String(a[byte=0:150]) + _random(60, 14)
    var extended = align_to_graph(graph, topology, prefix.as_bytes(), s, EXTENSION)
    assert_true(extended.score >= 2 * 150)
    assert_equal(extended.nodes[0], Int(graph.paths[0][0]))


def _paths(graph: Graph, node: Int, mut prefix: List[UInt8], mut out: List[String]):
    prefix.append(graph.bases[node])
    if len(graph.successors(node)) == 0:
        out.append(String(unsafe_from_utf8=prefix.copy()))
    for s in graph.successors(node):
        _paths(graph, s, prefix, out)
    _ = prefix.pop()


def test_graph_score_is_the_best_path_score() raises:
    var s = Scoring()
    for case_id in range(10):
        var base = _random(60, 100 + case_id)
        var graph = Graph()
        _ = graph.add("s0", base.as_bytes())
        for k in range(1, 4):
            var next = _mutated(base, 25, 200 + 10 * case_id + k)
            var topology = Topology(graph)
            var found = align_to_graph(graph, topology, next.as_bytes(), s, GLOBAL, -1)
            _ = graph.fuse_nodes(String(t"s{k}"), next.as_bytes(), found.nodes)
        var query = _mutated(base, 7, 300 + case_id)
        var spelled = List[String]()
        for node in range(len(graph)):
            if len(graph.predecessors(node)) == 0:
                var prefix = List[UInt8]()
                _paths(graph, node, prefix, spelled)
        var want = Int.MIN
        for path in spelled:
            want = max(want, _score_from_cost(path, query, s))
        var topology = Topology(graph)
        assert_equal(align_to_graph(graph, topology, query.as_bytes(), s, GLOBAL, -1).score, want)
        assert_equal(align_to_graph(graph, topology, query.as_bytes(), s, GLOBAL, 10, 0.01).score, want)


def test_16_and_32_bit_lanes_agree() raises:
    for mode in [GLOBAL, LOCAL, EXTENSION]:
        for case_id in range(6):
            var base = _random(150 + 90 * case_id, 400 + case_id)
            var graph = Graph()
            _ = graph.add("s0", base.as_bytes())
            for k in range(1, 5):
                var next = _mutated(base, 7, 500 + 10 * case_id + k)
                var topology = Topology(graph)
                var found = align_to_graph(graph, topology, next.as_bytes(), Scoring(), GLOBAL)
                _ = graph.fuse_nodes(String(t"s{k}"), next.as_bytes(), found.nodes)
            var query = _mutated(base, 6, 600 + case_id)
            var topology = Topology(graph)
            var row_of = List[Int](length=len(graph), fill=-1)
            var narrow = Workspace()
            var wide = Workspace()
            var a = narrow._align[DType.int16](
                graph, topology, topology.order, SOURCE, SINK, query.as_bytes(), Scoring(), row_of, mode
            )
            var b = wide._align[DType.int32](
                graph, topology, topology.order, SOURCE, SINK, query.as_bytes(), Scoring(), row_of, mode
            )
            assert_equal(a.score, b.score)
            for j in range(len(a.nodes)):
                assert_equal(a.nodes[j], b.nodes[j])


def test_a_long_deletion_reaches_the_end() raises:
    """A read holding only a long path's two ends: one deletion of most of the path, a score far below what the
    matches earn, which 16-bit lanes must still tell from an unreached cell."""
    var path = _random(12000, 41)
    var letters = List[UInt8](path.as_bytes()[0:200])
    letters.extend(List[UInt8](path.as_bytes()[11800:12000]))
    var query = String(unsafe_from_utf8=letters)
    var graph = Graph()
    _ = graph.add("path", path.as_bytes())
    var topology = Topology(graph)
    var found = align_to_graph(graph, topology, query.as_bytes(), Scoring(), GLOBAL, -1)
    assert_equal(len(found.nodes), 400)
    assert_equal(found.score, _score_from_cost(path, query, Scoring()))


def test_a_capped_band_keeps_the_alignment() raises:
    """The band capped (`band_max`) far below abPOA's b + f L on a 3 kb read: the same score as the uncapped band in
    fewer cells, and a deletion of most of a long path still found, the band following the path's diagonal."""
    var base = _random(3000, 71)
    var graph = Graph()
    _ = graph.add("s0", base.as_bytes())
    for k in range(1, 5):
        var next = _mutated(base, 9, 700 + k)
        var topology = Topology(graph)
        var found = align_to_graph(graph, topology, next.as_bytes(), Scoring(), GLOBAL)
        _ = graph.fuse_nodes(String(t"s{k}"), next.as_bytes(), found.nodes)
    var letters = List[UInt8](_mutated(base, 8, 800).as_bytes())
    var query = List[UInt8](letters[0:1500])
    query.extend(List[UInt8](letters[1540 : len(letters)]))
    var topology = Topology(graph)
    var row_of = List[Int](length=len(graph), fill=-1)
    var workspace = Workspace()
    var capped = workspace.align(
        graph, topology, topology.order, SOURCE, SINK, Span(query), Scoring(), row_of, GLOBAL, 10, 0.01, 16
    )
    var uncapped = workspace.align(
        graph, topology, topology.order, SOURCE, SINK, Span(query), Scoring(), row_of, GLOBAL, 10, 0.01, -1
    )
    assert_equal(capped.score, uncapped.score)
    assert_true(capped.cells < uncapped.cells)

    var path = _random(12000, 43)
    var ends = List[UInt8](path.as_bytes()[0:200])
    ends.extend(List[UInt8](path.as_bytes()[11800:12000]))
    var line = Graph()
    _ = line.add("path", path.as_bytes())
    var order = Topology(line)
    var rows = List[Int](length=len(line), fill=-1)
    var found = workspace.align(
        line, order, order.order, SOURCE, SINK, Span(ends), Scoring(), rows, GLOBAL, 10, 0.01, 16
    )
    assert_equal(len(found.nodes), 400)
    assert_equal(found.score, _score_from_cost(path, String(unsafe_from_utf8=ends), Scoring()))


def test_prefix_and_suffix_modes_place_partial_reads() raises:
    """A read fixed at the graph's start ends anywhere, its letters past a sink inserted (an overhang); one fixed at
    the end starts anywhere, its letters before a source inserted."""
    var path = _random(3000, 47)
    var graph = Graph()
    _ = graph.add("path", path.as_bytes())
    var topology = Topology(graph)
    var bytes = path.as_bytes()
    # An exact prefix: letter j on node j.
    var head = List[UInt8](bytes[0:1200])
    var found = align_to_graph(graph, topology, Span(head), Scoring(), PREFIX)
    for j in range(len(head)):
        assert_equal(found.nodes[j], j)
    # Past the end: the path's letters on its nodes, the overhang inserted.
    var long = List[UInt8](bytes[0:3000])
    long.extend(List[UInt8](_random(200, 53).as_bytes()))
    found = align_to_graph(graph, topology, Span(long), Scoring(), PREFIX)
    for j in range(3000):
        assert_equal(found.nodes[j], j)
    for j in range(3000, 3200):
        assert_equal(found.nodes[j], -1)
    # An exact suffix: letter j on node 1800 + j; before the start, an overhang inserted.
    var tail = List[UInt8](bytes[1800:3000])
    found = align_to_graph(graph, topology, Span(tail), Scoring(), SUFFIX)
    for j in range(len(tail)):
        assert_equal(found.nodes[j], 1800 + j)
    var early = List[UInt8](_random(150, 59).as_bytes())
    early.extend(List[UInt8](bytes[0:3000]))
    found = align_to_graph(graph, topology, Span(early), Scoring(), SUFFIX)
    for j in range(150):
        assert_equal(found.nodes[j], -1)
    for j in range(3000):
        assert_equal(found.nodes[150 + j], j)


def test_16_bit_rows_far_apart_in_base_keep_the_alignment() raises:
    """A 20 kb path and a 30-letter source chain joining it near its end, the path aligned to it in SUFFIX mode: the
    chain's rows have bases some 38,000 below the path's at the join, past what 16 bits hold. Narrowed, the
    difference wrapped and the chain's scores came in far too high (67,510, past the 40,000 of matching every
    letter); 16-bit lanes must give 32-bit lanes' alignment."""
    var a = _random(20000, 5)
    var graph = Graph()
    _ = graph.add("a", a.as_bytes())
    var b = List[UInt8](_random(30, 6).as_bytes())
    var aligned = List[Int](length=30, fill=-1)
    for i in range(19000, 20000):
        b.append(a.as_bytes()[i])
        aligned.append(Int(graph.paths[0][i]))
    _ = graph.fuse_nodes("b", Span(b), aligned)
    var topology = Topology(graph)
    var row_of = List[Int](length=len(graph), fill=-1)
    var narrow = Workspace()
    var wide = Workspace()
    var got = narrow.align(graph, topology, topology.order, SOURCE, SINK, a.as_bytes(), Scoring(), row_of, SUFFIX)
    var want = wide._align[DType.int32](
        graph, topology, topology.order, SOURCE, SINK, a.as_bytes(), Scoring(), row_of, SUFFIX, 10, 0.01
    )
    assert_equal(got.score, want.score)
    assert_equal(got.score, 2 * 20000)


def test_a_large_mismatch_penalty_takes_32_bit_lanes() raises:
    """The profile holds -mismatch: past some 16,000 it wrapped 16-bit lanes in every mode."""
    var base = _random(40, 11)
    var query = _mutated(base, 7, 12)
    var graph = Graph()
    _ = graph.add("a", base.as_bytes())
    var topology = Topology(graph)
    var row_of = List[Int](length=len(graph), fill=-1)
    for mismatch in [20000, 40000]:
        var s = Scoring()
        s.mismatch = Int32(mismatch)
        for mode in [GLOBAL, LOCAL, EXTENSION, PREFIX, SUFFIX]:
            var picked = Workspace()
            var wide = Workspace()
            var got = picked.align(graph, topology, topology.order, SOURCE, SINK, query.as_bytes(), s, row_of, mode)
            var want = wide._align[DType.int32](
                graph, topology, topology.order, SOURCE, SINK, query.as_bytes(), s, row_of, mode, 10, 0.01
            )
            assert_equal(got.score, want.score)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
