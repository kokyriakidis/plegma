from std.random import random_ui64, seed
from std.testing import assert_equal, assert_true, TestSuite

from plegma.refine import gap_opens, objective, refine, sum_of_pairs


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


def _flush_left(seqs: List[List[UInt8]]) -> List[List[UInt8]]:
    """A poor MSA: every sequence from the first column, padded with gaps on the right."""
    var width = 0
    for s in seqs:
        width = max(width, len(s))
    var rows = List[List[UInt8]]()
    for s in seqs:
        var row = s.copy()
        while len(row) < width:
            row.append(UInt8(ord("-")))
        rows.append(row^)
    return rows^


def test_refinement_never_raises_the_sum_and_keeps_the_sequences() raises:
    seed(17)
    for trial in range(30):
        var truth = _random(Int(random_ui64(20, 300)))
        var seqs = List[List[UInt8]]()
        for _ in range(Int(random_ui64(3, 12))):
            seqs.append(_mutated(truth, 0.1))
        seqs.append(seqs[0].copy())
        var rows = _flush_left(seqs)
        var refined = refine(rows, 10, 32)
        assert_true(sum_of_pairs(refined) <= sum_of_pairs(rows), String(trial))
        for k in range(len(seqs)):
            assert_true(_ungapped(refined[k]) == seqs[k], String(trial))
            assert_equal(len(refined[k]), len(refined[0]))
        # Identical sequences share their row.
        assert_true(refined[0] == refined[len(seqs) - 1])


def test_refinement_mends_a_shifted_row() raises:
    # Three copies of one sequence, one shifted right by a gap run in front: refined, all three stack.
    var seq = List[UInt8]("ACGTTGCAACGT".as_bytes())
    var shifted = List[UInt8]("---".as_bytes())
    shifted.extend(seq.copy())
    var padded = seq.copy()
    padded.extend(List[UInt8]("---".as_bytes()))
    var rows = List[List[UInt8]]()
    rows.append(padded.copy())
    rows.append(shifted^)
    var other = List[UInt8]("ACGTTGCTACGT---".as_bytes())
    rows.append(other^)
    var refined = refine(rows, 5, 8)
    assert_equal(sum_of_pairs(refined), 2)


def test_refinement_past_the_partner_cap() raises:
    # Thirty distinct sequences: each one's anchors come from its nearest 24 alone.
    seed(19)
    var truth = _random(200)
    var seqs = List[List[UInt8]]()
    for _ in range(30):
        seqs.append(_mutated(truth, 0.08))
    var rows = _flush_left(seqs)
    var refined = refine(rows, 5, 32, gap_open=1)
    assert_true(sum_of_pairs(refined) < sum_of_pairs(rows))
    for k in range(len(seqs)):
        assert_true(_ungapped(refined[k]) == seqs[k])


def test_group_moves_shift_rows_together() raises:
    # Four sequences stacked one way and two another, the four ten columns right of the two: the refinement with
    # group moves stacks all six.
    seed(23)
    var core = _random(60)
    var seqs = List[List[UInt8]]()
    for k in range(6):
        var s = core.copy()
        var at = k * 7 if k < 4 else 40 + k * 3
        s[at] = UInt8(ord("A")) if s[at] != UInt8(ord("A")) else UInt8(ord("C"))
        seqs.append(s^)
    var rows = List[List[UInt8]]()
    for k in range(6):
        var pad = List[UInt8](length=10, fill=UInt8(ord("-")))
        var row = List[UInt8]()
        if k < 4:
            row.extend(pad.copy())
            row.extend(seqs[k].copy())
        else:
            row.extend(seqs[k].copy())
            row.extend(pad.copy())
        rows.append(row^)
    var together = refine(rows, 5, 32, gap_open=1, group_band=32)
    assert_equal(sum_of_pairs(together), sum_of_pairs(seqs))
    for k in range(6):
        assert_true(_ungapped(together[k]) == seqs[k])


def test_group_moves_never_raise_the_objective() raises:
    seed(29)
    for trial in range(20):
        var truth = _random(Int(random_ui64(20, 200)))
        var seqs = List[List[UInt8]]()
        for _ in range(Int(random_ui64(3, 10))):
            seqs.append(_mutated(truth, 0.1))
        var rows = _flush_left(seqs)
        var refined = refine(rows, 5, 16, gap_open=1, group_band=16)
        assert_true(sum_of_pairs(refined) + gap_opens(refined) <= sum_of_pairs(rows) + gap_opens(rows), String(trial))
        for k in range(len(seqs)):
            assert_true(_ungapped(refined[k]) == seqs[k], String(trial))
            assert_equal(len(refined[k]), len(refined[0]))


def test_copy_aware_moves_keep_the_sequences_and_never_raise_the_objective() raises:
    # Arrays of one unit differing by whole copies, flush left: with the period, the refinement's bands reach a copy
    # and its group moves split the rows by copy number.
    seed(59)
    var unit = _random(11)
    var left = _random(30)
    var right = _random(30)
    for trial in range(6):
        var seqs = List[List[UInt8]]()
        for _ in range(Int(random_ui64(3, 7))):
            var s = left.copy()
            for _ in range(Int(random_ui64(5, 12))):
                s.extend(_mutated(unit, 0.05))
            s.extend(right.copy())
            seqs.append(s^)
        var rows = _flush_left(seqs)
        var refined = refine(rows, 4, 8, gap_open=1, period=11)
        assert_true(sum_of_pairs(refined) + gap_opens(refined) <= sum_of_pairs(rows) + gap_opens(rows), String(trial))
        for k in range(len(seqs)):
            assert_true(_ungapped(refined[k]) == seqs[k], String(trial))
            assert_equal(len(refined[k]), len(refined[0]))


def test_priced_refinement_never_raises_its_objective() raises:
    # Mismatches, gap letters and gap openings priced apart: the refinement lowers that objective, and at 1, 1 the
    # objective is the unit sum of pairs and the gaps opened.
    seed(61)
    for trial in range(12):
        var truth = _random(Int(random_ui64(30, 200)))
        var seqs = List[List[UInt8]]()
        for _ in range(Int(random_ui64(3, 8))):
            seqs.append(_mutated(truth, 0.1))
        var rows = _flush_left(seqs)
        assert_equal(objective(rows, 1), sum_of_pairs(rows) + gap_opens(rows))
        for prices in [(2, 2, 1), (4, 2, 6), (10, 7, 9)]:
            var refined = refine(rows, 3, 16, gap_open=prices[2], mismatch=prices[0], gap=prices[1], group_band=16)
            assert_true(
                objective(refined, prices[2], List[Int](), prices[0], prices[1])
                <= objective(rows, prices[2], List[Int](), prices[0], prices[1]),
                String(trial),
            )
            for k in range(len(seqs)):
                assert_true(_ungapped(refined[k]) == seqs[k], String(trial))


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
