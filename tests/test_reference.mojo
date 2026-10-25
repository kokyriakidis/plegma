from std.collections import Set
from std.random import random_ui64, seed
from std.testing import assert_equal, assert_true, TestSuite

from karan import Aligner, Costs

from plegma.pins import pinned_pairing, pins_in_context
from plegma.reference import on_reference


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


def _long_change_places(rows: List[List[UInt8]], least: Int) -> Int:
    """The distinct columns where some row's insertion or deletion of `least` letters or more against `rows[0]`
    starts (an insertion: letters facing the reference's gaps)."""
    var gap = UInt8(ord("-"))
    var starts = Set[Int]()
    for k in range(1, len(rows)):
        var c = 0
        while c < len(rows[0]):
            var deleting = rows[0][c] != gap and rows[k][c] == gap
            var inserting = rows[0][c] == gap and rows[k][c] != gap
            if not deleting and not inserting:
                c += 1
                continue
            var end = c
            var letters = 0
            while end < len(rows[0]):
                if deleting and rows[0][end] != gap and rows[k][end] == gap:
                    letters += 1
                elif inserting and rows[0][end] == gap and rows[k][end] != gap:
                    letters += 1
                elif not (rows[0][end] == gap and rows[k][end] == gap):
                    break
                end += 1
            if letters >= least:
                starts.add(c * 2 + Int(inserting))
            c = end
    return len(starts)


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
    var rows = on_reference(seqs, Costs.two_piece(40, 78, 7, 162, 3))
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


def test_on_reference_aligns_distant_sequences_at_the_fallback_costs() raises:
    seed(71)
    # Reads of one sequence at 12% error, about a quarter apart pairwise: past asm5's divergence, whose costs delete
    # and insert each again whole; the fallback costs align them letter by letter.
    var truth = _random(600)
    var seqs = List[List[UInt8]]()
    for _ in range(6):
        seqs.append(_mutated(truth, 0.12))
    var asm5 = Costs.two_piece(40, 78, 7, 162, 3)
    var plain = on_reference(seqs, asm5)
    var rows = on_reference(seqs, asm5, fallback=Costs.two_piece(24, 16, 12, 96, 5))
    for k in range(len(seqs)):
        assert_true(_ungapped(rows[k]) == seqs[k])
        assert_equal(len(rows[k]), len(rows[0]))
    assert_true(_long_change_places(plain, 50) > 0)
    assert_equal(_long_change_places(rows, 50), 0)
    assert_true(len(rows[0]) < len(plain[0]))
    # Alleles: one of six with a stretch of its own in place of the reference's keeps asm5's alignment.
    var left = _random(300)
    var right = _random(300)
    var alleles = List[List[UInt8]]()
    alleles.append(_array(left, _random(120), 1, right, 0.0))
    for _ in range(5):
        alleles.append(_mutated(alleles[0], 0.005))
    alleles.append(_array(left, _random(120), 1, right, 0.0))
    var kept = on_reference(alleles, asm5)
    var with_fallback = on_reference(alleles, asm5, fallback=Costs.two_piece(24, 16, 12, 96, 5))
    assert_true(_long_change_places(kept, 50) > 0)
    for k in range(len(alleles)):
        assert_true(with_fallback[k] == kept[k])


def test_pinned_alignment_spells_and_pins_identical_sequences_on_the_diagonal() raises:
    seed(73)
    var left = _random(1500)
    var right = _random(1500)
    var unit = _random(37)
    var reference = _array(_random(300), unit, 30, _random(300), 0.02)
    var found = pins_in_context(reference, reference, left, right)
    assert_true(len(found.r) > 0)
    for m in range(len(found.r)):
        assert_equal(found.r[m], found.q[m])
    # Copies added and removed: each pairing runs forward through both, every sequence spelled by its row.
    var seqs = List[List[UInt8]]()
    seqs.append(reference.copy())
    for copies in [22, 41, 30, 35]:
        seqs.append(_array(_random(300), unit, copies, _random(300), 0.02))
    var flanks = List[List[UInt8]]()
    flanks.append(left.copy())
    flanks.append(right.copy())
    var aligner = Aligner()
    for k in range(1, len(seqs)):
        var pairing = pinned_pairing(aligner, seqs[0], seqs[k], left, right, Costs.two_piece(40, 78, 7, 162, 3))
        assert_equal(len(pairing), len(seqs[k]))
        var last = -1
        for j in range(len(pairing)):
            if pairing[j] >= 0:
                assert_true(pairing[j] > last)
                last = pairing[j]
    var rows = on_reference(seqs, Costs.two_piece(40, 78, 7, 162, 3), flanks=flanks)
    for k in range(len(seqs)):
        assert_true(_ungapped(rows[k]) == seqs[k])
        assert_equal(len(rows[k]), len(rows[0]))


def test_on_reference_rows_do_not_depend_on_the_threads() raises:
    seed(79)
    var left = _random(1500)
    var right = _random(1500)
    var unit = _random(29)
    var seqs = List[List[UInt8]]()
    seqs.append(_array(_random(200), unit, 25, _random(200), 0.02))
    for copies in [18, 25, 31, 22, 27, 12]:
        seqs.append(_array(_random(200), unit, copies, _random(200), 0.02))
    var flanks = List[List[UInt8]]()
    flanks.append(left.copy())
    flanks.append(right.copy())
    var asm5 = Costs.two_piece(40, 78, 7, 162, 3)
    var exact = on_reference(seqs, asm5)
    var pinned = on_reference(seqs, asm5, flanks=flanks)
    for threads in [2, 4]:
        assert_true(on_reference(seqs, asm5, threads=threads) == exact)
        assert_true(on_reference(seqs, asm5, threads=threads, flanks=flanks) == pinned)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
