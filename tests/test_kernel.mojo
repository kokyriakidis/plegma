from std.testing import assert_equal, assert_true, TestSuite

from karan import Costs, distance

from plegma.graph import Graph
from plegma.kernel import EXTENSION, GLOBAL, LOCAL, Scoring, Topology, align_to_graph


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


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
