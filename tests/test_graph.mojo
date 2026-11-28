from std.testing import assert_equal, assert_true, TestSuite

from karan import Aligner, Costs, Mode, Ties

from plegma.graph import DELETION, INSERTION, Graph, cigar_runs


def _bytes(text: StringSpan) -> List[UInt8]:
    return List[UInt8](text.as_bytes())


def _ungapped(row: String) -> String:
    var out = List[UInt8]()
    for c in row.as_bytes():
        if c != UInt8(ord("-")):
            out.append(c)
    return String(unsafe_from_utf8=out)


def _sequences() -> List[List[UInt8]]:
    """A tandem array between two flanks, in five copy numbers, one copy changed in each and one letter in one."""
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
    return sequences^


def test_cigar_runs() raises:
    var runs = cigar_runs("4=1D2=12I6X")
    assert_equal(len(runs.ops), 5)
    assert_equal(runs.counts[3], 12)
    assert_equal(runs.ops[4], UInt8(ord("X")))


def test_star_msa_spells_every_sequence() raises:
    """Each sequence fused along its alignment to the first one's path: every row spells its sequence and the
    shared flanks stack."""
    var sequences = _sequences()
    var graph = Graph()
    var aligner = Aligner()
    _ = graph.add("s0", Span(sequences[0]))
    for k in range(1, len(sequences)):
        var found = aligner.align(Span(sequences[0]), Span(sequences[k]), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
        var runs = cigar_runs(found.cigar)
        var aligned = List[Int]()
        var i = 0
        for x in range(len(runs.ops)):
            for _ in range(runs.counts[x]):
                if runs.ops[x] == DELETION:
                    i += 1
                elif runs.ops[x] == INSERTION:
                    aligned.append(-1)
                else:
                    aligned.append(Int(graph.paths[0][i]))
                    i += 1
        _ = graph.fuse_nodes(String(t"s{k}"), Span(sequences[k]), aligned)
    var rows = graph.msa()
    assert_equal(len(rows), len(sequences))
    for k in range(len(rows)):
        assert_equal(rows[k].byte_length(), rows[0].byte_length())
        assert_true(_ungapped(rows[k]) == String(unsafe_from_utf8=sequences[k].copy()))
    # The shared flanks stack: the first columns hold one letter each across all rows.
    for c in range(10):
        for k in range(1, len(rows)):
            assert_equal(rows[k].as_bytes()[c], rows[0].as_bytes()[c])


def test_from_msa_keeps_the_rows() raises:
    """An MSA's graph spells each row as its path, and its own MSA is the same rows."""
    var names: List[String] = ["a", "b", "c"]
    var rows: List[String] = ["ACG-TTA", "A-GCTTA", "ACGCT-A"]
    var graph = Graph.from_msa(names, rows)
    assert_equal(len(graph.paths), 3)
    var again = graph.msa()
    for k in range(3):
        assert_true(again[k] == rows[k])
    # A column holds one node per letter: ACG-TTA and ACGCT-A share the C of column 1.
    assert_equal(graph.paths[0][1], graph.paths[2][1])


def test_heaviest_path_follows_the_majority() raises:
    var names: List[String] = ["a", "b", "c"]
    var rows: List[String] = ["ACGTA", "ACTTA", "ACGTA"]
    var graph = Graph.from_msa(names, rows)
    var rank = graph.node_order()
    var order = List[Int](length=len(rank), fill=0)
    for node in range(len(rank)):
        order[rank[node]] = node
    var everyone: List[Int] = [0, 1, 2]
    var letters = List[UInt8]()
    for node in graph.heaviest_path(order, everyone):
        letters.append(graph.bases[node])
    assert_true(String(unsafe_from_utf8=letters) == "ACGTA")
    # Only b's path weighs: its own letters.
    var b_only: List[Int] = [1]
    letters.clear()
    for node in graph.heaviest_path(order, b_only):
        letters.append(graph.bases[node])
    assert_true(String(unsafe_from_utf8=letters) == "ACTTA")


def test_from_msa_refuses_rows_of_different_lengths() raises:
    """A shorter row would be read past its end, a longer one cut short: both raise."""
    var names: List[String] = ["a", "b"]
    for rows in [List[String](["ACGTAC", "ACG"]), List[String](["ACG", "ACGTAC"])]:
        var raised = False
        try:
            _ = Graph.from_msa(names, rows)
        except:
            raised = True
        assert_true(raised)


def test_heaviest_path_spanning_counts_past_32_bits() raises:
    """Edge weights scaled by the paths' count pass 32 bits at some 46,000 paths: the heaviest path stays the
    spanning reads' all the same."""
    var graph = Graph()
    var spanning = List[Int]()
    for _ in range(47000):
        spanning.append(graph.add("r", "ACGTACGTAC".as_bytes()))
    var left: List[Int] = [0, 1, 2, 3, 4]
    _ = graph.fuse_nodes("l", "ACGTA".as_bytes(), left)
    var rank = graph.node_order()
    var order = List[Int](length=len(rank), fill=0)
    for node in range(len(rank)):
        order[rank[node]] = node
    assert_equal(len(graph.heaviest_path_spanning(order, spanning)), 10)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
