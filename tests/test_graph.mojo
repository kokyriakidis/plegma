from std.testing import assert_equal, assert_true, TestSuite

from karan import Aligner, Costs, distance

from plegma.align import anchored_alignment
from plegma.backbone import Level, MUMS
from plegma.chain import ChainCosts
from plegma.graph import Graph, cigar_runs
from plegma.poa import GapAffine, SINK, SOURCE, align_window


def _bytes(text: StringSlice) -> List[UInt8]:
    return List[UInt8](text.as_bytes())


def _ungapped(row: String) -> String:
    var out = List[UInt8]()
    for c in row.as_bytes():
        if c != UInt8(ord("-")):
            out.append(c)
    return String(unsafe_from_utf8=out)


def test_cigar_runs() raises:
    var runs = cigar_runs("4=1D2=12I6X")
    assert_equal(len(runs.ops), 5)
    assert_equal(runs.counts[3], 12)
    assert_equal(runs.ops[4], UInt8(ord("X")))


def test_star_msa_spells_every_sequence() raises:
    var unit = "TGATGGTGGTGGTTGGTAG"
    var flank_left = "CCATTAGGACTTGACCGATAGGCAT"
    var flank_right = "GGATCCTTAGAGCTCAATCGGA"
    var sequences = List[List[UInt8]]()
    for copies in [9, 6, 7, 9, 4]:
        var text = String(flank_left)
        for c in range(copies):
            text += unit if c != 2 else "TGATGGTGCTGGTTGGTAG"
        text += flank_right
        sequences.append(_bytes(text))
    sequences[3][40] = UInt8(ord("A"))
    var levels: List[Level] = [Level(MUMS, 12)]
    for config in range(2):
        var graph = Graph()
        var aligner = Aligner()
        _ = graph.add("s0", Span(sequences[0]))
        for i in range(1, len(sequences)):
            var found = anchored_alignment(
                aligner,
                Span(sequences[0]),
                Span(sequences[i]),
                levels.copy() if config == 1 else List[Level](),
                ChainCosts.chainx(),
                stay=True,
            )
            _ = graph.fuse(String(t"s{i}"), Span(sequences[i]), 0, found.ops, found.counts)
        var rows = graph.msa()
        assert_equal(len(rows), len(sequences))
        for i in range(len(rows)):
            assert_equal(rows[i].byte_length(), rows[0].byte_length())
            assert_true(_ungapped(rows[i]) == String(unsafe_from_utf8=sequences[i].copy()))
        # The shared flanks stack: the first columns hold one letter each across all rows.
        for c in range(10):
            for i in range(1, len(rows)):
                assert_equal(rows[i].as_bytes()[c], rows[0].as_bytes()[c])


def _order(graph: Graph) raises -> Tuple[List[Int], List[Int]]:
    var rank = graph.node_order()
    var order = List[Int](length=len(rank), fill=0)
    for node in range(len(rank)):
        order[rank[node]] = node
    return (rank^, order^)


def test_window_on_one_path_matches_karan() raises:
    var texts: List[String] = [
        "ACGTTGCAAGGCTTAGGACCATGACGTA",
        "ACGTGCAAGGCTTTTAGGACCATGACGTA",
        "TTTTACGTTGCAAGG",
        "ACGTTGCAAGGCTTAGGACCATGACGTAACGTTGCAAGGCTTAGGACCATGACGTA",
        "AC",
    ]
    for a in texts:
        for b in texts:
            var graph = Graph()
            _ = graph.add("a", a.as_bytes())
            var ro = _order(graph)
            var found = align_window(graph, ro[0], ro[1], SOURCE, SINK, b.as_bytes(), GapAffine(4, 6, 2), 1 << 30)
            assert_equal(found.value().cost, distance(a, b, Costs.affine(4, 6, 2)))
            # Through the window's own nodes: the same alignment between the first and last letters.
            if a.byte_length() > 2 and b.byte_length() > 2:
                var path = graph.paths[0].copy()
                var inner = align_window(
                    graph,
                    ro[0],
                    ro[1],
                    path[0],
                    path[len(path) - 1],
                    b.as_bytes()[1 : b.byte_length() - 1],
                    GapAffine(4, 6, 2),
                    1 << 30,
                )
                var a_inner = String(a[byte = 1 : a.byte_length() - 1])
                var b_inner = String(b[byte = 1 : b.byte_length() - 1])
                assert_equal(inner.value().cost, distance(a_inner, b_inner, Costs.affine(4, 6, 2)))


def test_window_takes_the_better_path() raises:
    var a = "GGCATTAGCATCAGGGATTACCTTAGCAGTCCA"
    var b = "GGCATTAGCATCATTTTTTTTTTTTCAGTCCA"
    var graph = Graph()
    var aligner = Aligner()
    _ = graph.add("a", a.as_bytes())
    var found = anchored_alignment(
        aligner, a.as_bytes(), b.as_bytes(), List[Level](), ChainCosts.chainx(), Costs.affine(4, 6, 2)
    )
    _ = graph.fuse("b", b.as_bytes(), 0, found.ops, found.counts)
    var ro = _order(graph)
    for query in [a, b, String("GGCATTAGCATCATTTTTTTTTTTTCAGTCCA")]:
        var got = align_window(graph, ro[0], ro[1], SOURCE, SINK, query.as_bytes(), GapAffine(4, 6, 2), 1 << 30)
        var best = min(distance(a, query, Costs.affine(4, 6, 2)), distance(b, query, Costs.affine(4, 6, 2)))
        assert_true(got.value().cost <= best)
        assert_equal(got.value().cost, 0)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
