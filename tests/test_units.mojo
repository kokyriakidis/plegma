from std.random import random_ui64, seed
from std.testing import assert_equal, assert_true, TestSuite

from karan import Costs

from plegma.progressive import on_reference, on_tree, progressive
from plegma.refine import gap_opens, sum_of_pairs
from plegma.units import tandem_layout, tandem_repeat


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


def _ungapped(row: List[UInt8]) -> List[UInt8]:
    var out = List[UInt8]()
    for x in row:
        if x != UInt8(ord("-")):
            out.append(x)
    return out^


def _array(left: List[UInt8], unit: List[UInt8], copies: Int, right: List[UInt8], rate: Float64) -> List[UInt8]:
    var out = left.copy()
    for _ in range(copies):
        out.extend(_mutated(unit, rate))
    out.extend(right.copy())
    return out^


def test_the_period_of_an_array_between_flanks() raises:
    seed(31)
    var left = _random(150)
    var right = _random(150)
    for period in [7, 27, 60]:
        var unit = _random(period)
        var seq = _array(left, unit, 1500 // period, right, 0.02)
        var r = tandem_repeat(seq)
        assert_equal(r.period, period)
        assert_true(r.start >= 140 and r.start <= 150 + 3 * period, String(r.start))
        assert_true(r.identity > 0.8)


def test_the_smallest_period_not_a_multiple() raises:
    # A unit of 10 repeated exactly matches itself 20 and 30 letters on as well: the period is 10.
    seed(37)
    var unit = _random(10)
    var seq = _array(_random(100), unit, 80, _random(100), 0.0)
    assert_equal(tandem_repeat(seq).period, 10)


def test_no_period_in_random_sequence() raises:
    seed(41)
    var seq = _random(2000)
    var r = tandem_repeat(seq)
    assert_true(r.period == 0 or r.end - r.start < 100)


def test_layout_spells_each_array_and_cuts_copies_in_phase() raises:
    seed(43)
    var left = _random(120)
    var right = _random(120)
    var unit = _random(23)
    var seqs = List[List[UInt8]]()
    for copies in [30, 34, 31, 38, 30, 36]:
        seqs.append(_array(left, unit, copies, right, 0.04))
    var layout = tandem_layout(seqs)
    assert_equal(layout.period, 23)
    assert_equal(len(layout.rows), len(seqs))
    for k in range(len(seqs)):
        assert_equal(len(layout.rows[k]), len(layout.rows[0]))
        var spelled = _ungapped(layout.rows[k])
        assert_equal(len(spelled), layout.ends[k] - layout.starts[k])
        for t in range(len(spelled)):
            assert_equal(spelled[t], seqs[k][layout.starts[k] + t])


def test_no_layout_without_an_array() raises:
    seed(47)
    var truth = _random(800)
    var seqs = List[List[UInt8]]()
    for _ in range(5):
        seqs.append(_mutated(truth, 0.05))
    assert_equal(len(tandem_layout(seqs).rows), 0)


def test_progressive_keeps_the_sequences_and_the_copy_start_helps() raises:
    seed(53)
    var left = _random(100)
    var right = _random(100)
    var unit = _random(19)
    var seqs = List[List[UInt8]]()
    for copies in [12, 14, 13, 16, 12]:
        seqs.append(_array(left, unit, copies, right, 0.04))
    seqs.append(seqs[1].copy())
    var both = progressive(seqs, 3, group_band=32)
    var placed = progressive(seqs, 3, group_band=32, units=False)
    for k in range(len(seqs)):
        assert_true(_ungapped(both[k]) == seqs[k])
        assert_equal(len(both[k]), len(both[0]))
        assert_true(_ungapped(placed[k]) == seqs[k])
    assert_true(both[1] == both[len(seqs) - 1])
    # The copy-level start is kept only when cheaper.
    assert_true(sum_of_pairs(both) + gap_opens(both) <= sum_of_pairs(placed) + gap_opens(placed))


def _gap_runs(a: List[UInt8], b: List[UInt8]) -> Int:
    """The runs of letters of `a` facing gaps in `b`, over the columns where either holds a letter."""
    var runs = 0
    var inside = False
    for c in range(len(a)):
        if a[c] == UInt8(ord("-")) and b[c] == UInt8(ord("-")):
            continue
        var gap = a[c] != UInt8(ord("-")) and b[c] == UInt8(ord("-"))
        if gap and not inside:
            runs += 1
        inside = gap
    return runs


def test_on_reference_keeps_each_deletion_whole_against_the_reference() raises:
    seed(59)
    var left = _random(80)
    var right = _random(80)
    var unit = _random(23)
    var reference = _array(left, unit, 20, right, 0.0)
    var seqs = List[List[UInt8]]()
    seqs.append(reference.copy())
    # Sequences with copies deleted, and some with copies added: against the reference one deletion each.
    for copies in [14, 9, 20, 26, 24, 11, 14]:
        seqs.append(_array(left, unit, copies, right, 0.0))
    var rows = on_reference(seqs, Costs.two_piece(24, 16, 12, 96, 5))
    for k in range(len(seqs)):
        assert_true(_ungapped(rows[k]) == seqs[k])
        assert_equal(len(rows[k]), len(rows[0]))
    for k in range(1, len(seqs)):
        assert_true(_gap_runs(rows[0], rows[k]) <= 1)
    # The same sequence twice, the same row.
    assert_true(rows[1] == rows[7])
    # The copies inserted at one place take no more columns than the longest insertion there.
    var extra = len(rows[0]) - len(reference)
    assert_true(extra <= 6 * 23)


def test_on_tree_keeps_each_allele_as_simple_against_the_reference() raises:
    seed(61)
    var left = _random(80)
    var right = _random(80)
    var unit = _random(23)
    var extra = _random(30)
    var reference = _array(left, unit, 20, right, 0.0)
    var family = _array(left, unit, 14, right, 0.0)
    # A member of the family with a stretch of its own inside the array: the family's deletion and an insertion.
    var member = _array(left, unit, 7, List[UInt8](), 0.0)
    member.extend(extra.copy())
    member.extend(_array(List[UInt8](), unit, 7, right, 0.0))
    var seqs = List[List[UInt8]]()
    seqs.append(reference.copy())
    seqs.append(member.copy())
    for _ in range(3):
        seqs.append(family.copy())
    for share in [0.0, 0.5]:
        var rows = on_tree(seqs, Costs.two_piece(24, 16, 12, 96, 5), share)
        for k in range(len(seqs)):
            assert_true(_ungapped(rows[k]) == seqs[k])
            assert_equal(len(rows[k]), len(rows[0]))
        assert_true(rows[2] == rows[4])
        # Against the reference the family is one deletion, the member that deletion and its insertion.
        assert_equal(_gap_runs(rows[0], rows[2]) + _gap_runs(rows[2], rows[0]), 1)
        assert_equal(_gap_runs(rows[0], rows[1]) + _gap_runs(rows[1], rows[0]), 2)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
