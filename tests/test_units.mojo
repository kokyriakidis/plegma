from std.random import random_ui64, seed
from std.testing import assert_equal, assert_true, TestSuite

from plegma.progressive import progressive
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


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
