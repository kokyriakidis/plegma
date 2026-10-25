from std.random import random_ui64, seed
from std.testing import assert_equal, assert_true, TestSuite

from karan import distance

from plegma.graph import Graph
from plegma.kernel import BOTH, LEFT, RIGHT
from plegma.poa import PoaOptions, poa
from plegma.polish import polish


def _random(length: Int) -> List[UInt8]:
    var out = List[UInt8](capacity=length)
    var letters = "ACGT".as_bytes()
    for _ in range(length):
        out.append(letters[Int(random_ui64(0, 3))])
    return out^


def _mutated(seq: List[UInt8], rate: Float64) -> List[UInt8]:
    var out = List[UInt8]()
    var letters = "ACGT".as_bytes()
    for c in seq:
        var r = Float64(random_ui64(0, 999_999)) / 1e6
        if r < rate / 3:
            out.append(letters[Int(random_ui64(0, 3))])
        elif r < 2 * rate / 3:
            continue
        elif r < rate:
            out.append(c)
            out.append(letters[Int(random_ui64(0, 3))])
        else:
            out.append(c)
    return out^


def _slice(seq: List[UInt8], start: Int, end: Int) -> List[UInt8]:
    var out = List[UInt8](capacity=end - start)
    for i in range(start, end):
        out.append(seq[i])
    return out^


def _assembled(truth: List[UInt8], n_both: Int, n_left: Int, n_right: Int, rate: Float64) raises -> List[UInt8]:
    """The polished consensus of reads of `truth`: `n_both` spanning it, `n_left` fixed at its start, `n_right` at
    its end."""
    var seqs = List[List[UInt8]]()
    var names = List[String]()
    var anchors = List[Int]()
    var length = len(truth)
    for k in range(n_both + n_left + n_right):
        var anchor = BOTH if k < n_both else (LEFT if k < n_both + n_left else RIGHT)
        var start = 0
        var end = length
        if anchor == LEFT:
            end = Int(random_ui64(UInt64(length * 3 // 10), UInt64(length * 95 // 100)))
        elif anchor == RIGHT:
            start = Int(random_ui64(UInt64(length // 20), UInt64(length * 7 // 10)))
        seqs.append(_mutated(_slice(truth, start, end), rate))
        names.append(String(t"r{k}"))
        anchors.append(anchor)
    var made = poa(names, seqs, anchors, PoaOptions(), Graph())
    # Every read is in the graph and its row spells it.
    var rows = made.graph.msa()
    for k in range(len(seqs)):
        var spelled = List[UInt8]()
        for c in rows[made.shown[k]].as_bytes():
            if c != UInt8(ord("-")):
                spelled.append(c)
        assert_true(spelled == seqs[k])
    assert_equal(len(made.consensus), 1)
    var everyone = List[Int]()
    for k in range(len(seqs)):
        everyone.append(k)
    return polish(made.consensus[0], seqs, everyone, 50, anchors=anchors)


def test_few_spanning_reads_assemble_the_truth() raises:
    seed(89)
    var truth = _random(3000)
    var consensus = _assembled(truth, 1, 8, 8, 0.05)
    assert_true(distance(String(unsafe_from_utf8=consensus), String(unsafe_from_utf8=truth)) <= 3)


def test_no_spanning_read_joins_the_two_sides() raises:
    seed(97)
    var truth = _random(3000)
    var consensus = _assembled(truth, 0, 10, 10, 0.05)
    assert_true(distance(String(unsafe_from_utf8=consensus), String(unsafe_from_utf8=truth)) <= 3)


def test_one_side_assembles_a_prefix_where_reads_agree() raises:
    seed(101)
    var truth = _random(3000)
    var consensus = _assembled(truth, 0, 12, 0, 0.05)
    # The consensus stops where fewer than three reads go on: a prefix of the truth, most of it.
    assert_true(len(consensus) >= 1000)
    var prefix = _slice(truth, 0, len(consensus))
    assert_true(distance(String(unsafe_from_utf8=consensus), String(unsafe_from_utf8=prefix)) <= 10)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
