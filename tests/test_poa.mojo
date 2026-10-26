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


def _assembled(
    truth: List[UInt8], n_both: Int, n_left: Int, n_right: Int, rate: Float64, overhang: Int = 0, pinned: Bool = True
) raises -> List[UInt8]:
    """The polished consensus of reads of `truth`: `n_both` spanning it, `n_left` fixed at its start, `n_right` at
    its end; with `overhang`, the one-sided reads run that many letters past the far end instead (into sequence
    beyond it), as a local assembler cuts them."""
    var seqs = List[List[UInt8]]()
    var names = List[String]()
    var anchors = List[Int]()
    var length = len(truth)
    var before = _random(overhang)
    var after = _random(overhang)
    for k in range(n_both + n_left + n_right):
        var anchor = BOTH if k < n_both else (LEFT if k < n_both + n_left else RIGHT)
        var start = 0
        var end = length
        if anchor == LEFT:
            end = Int(random_ui64(UInt64(length * 3 // 10), UInt64(length * 95 // 100)))
        elif anchor == RIGHT:
            start = Int(random_ui64(UInt64(length // 20), UInt64(length * 7 // 10)))
        var piece = _slice(truth, start, end)
        if overhang > 0 and anchor == LEFT:
            piece = truth + after.copy()
        elif overhang > 0 and anchor == RIGHT:
            piece = before + truth.copy()
        seqs.append(_mutated(piece, rate))
        names.append(String(t"r{k}"))
        anchors.append(anchor)
    var options = PoaOptions()
    options.seeded = pinned
    var made = poa(names, seqs, anchors, options, Graph())
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


def test_reads_past_the_far_anchor_do_not_extend_the_consensus() raises:
    # More one-sided reads than spanning ones, each running past the far anchor: the consensus still ends where the
    # spanning reads do.
    seed(103)
    var truth = _random(3000)
    var consensus = _assembled(truth, 3, 8, 8, 0.05, overhang=330)
    assert_true(distance(String(unsafe_from_utf8=consensus), String(unsafe_from_utf8=truth)) <= 3)


def test_reads_past_the_far_anchor_unpinned() raises:
    # The same, each read aligned to the whole graph: the left reads' letters past the end join the graph after it,
    # and the right reads still end at the right anchor, not at the end of those letters.
    seed(103)
    var truth = _random(3000)
    var consensus = _assembled(truth, 3, 8, 8, 0.05, overhang=330, pinned=False)
    assert_true(distance(String(unsafe_from_utf8=consensus), String(unsafe_from_utf8=truth)) <= 3)


def test_another_copys_one_sided_reads_get_no_vote() raises:
    # Five spanning reads of the truth, two left reads of it, and eight left reads of another copy (2% of its letters
    # substituted), outvoting them: the consensus is the truth's all the same.
    seed(107)
    var truth = _random(2000)
    var other = truth.copy()
    var letters = "ACGT".as_bytes()
    for i in range(0, len(other), 50):
        other[i] = letters[(Int(random_ui64(1, 3)) + Int(_code(other[i]))) % 4]
    var seqs = List[List[UInt8]]()
    var names = List[String]()
    var anchors = List[Int]()
    for k in range(15):
        if k >= 5 and k < 13:
            seqs.append(_mutated(other, 0.01))
        else:
            seqs.append(_mutated(truth, 0.01))
        names.append(String(t"r{k}"))
        anchors.append(BOTH if k < 5 else LEFT)
    var made = poa(names, seqs, anchors, PoaOptions(), Graph())
    assert_equal(len(made.foreign), 8)
    for r in made.foreign:
        assert_true(r >= 5 and r < 13)
    assert_true(distance(String(unsafe_from_utf8=made.consensus[0]), String(unsafe_from_utf8=truth)) <= 3)


def test_a_spanning_read_wins_a_tie() raises:
    # One read on both anchors, one fixed at the right anchor alone that differs from it at three places (a paralog's,
    # missing the left anchor): one vote each, and the read on both anchors is the surer.
    seed(109)
    var truth = _random(400)
    var other = truth.copy()
    var letters = "ACGT".as_bytes()
    for i in [250, 300, 350]:
        other[i] = letters[(_code(other[i]) + 1) % 4]
    var seqs: List[List[UInt8]] = [truth.copy(), _slice(other, 150, 400)]
    var names: List[String] = ["both", "right"]
    var anchors: List[Int] = [BOTH, RIGHT]
    var made = poa(names, seqs, anchors, PoaOptions(), Graph())
    assert_true(made.consensus[0] == truth)


def _code(letter: UInt8) -> Int:
    if letter == UInt8(ord("A")):
        return 0
    if letter == UInt8(ord("C")):
        return 1
    return 2 if letter == UInt8(ord("G")) else 3


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
