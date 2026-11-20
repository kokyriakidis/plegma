from std.collections import Set
from std.random import random_ui64, seed
from std.testing import assert_equal, assert_true, TestSuite

from karan import Aligner, Costs

from sporos import read_fasta

from plegma.pins import _left_aligned, pinned_pairing, pins_in_context
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


def _b(text: StringSlice) -> List[UInt8]:
    return List[UInt8](text.as_bytes())


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


def test_pinned_past_the_pins_reach_raises() raises:
    """Past 2^25 letters with the flanks the seeds' positions do not fit: pinning raises rather than leaving one
    exact alignment of the whole (some 2 GiB of DP)."""
    var long = List[UInt8](length=1 << 25, fill=UInt8(ord("A")))
    var flank = List[UInt8](length=100, fill=UInt8(ord("C")))
    var aligner = Aligner(warm=False)
    var costs = Costs.two_piece(40, 78, 7, 162, 3)
    var raised = False
    try:
        _ = pinned_pairing(aligner, long, long, flank, flank, costs)
    except e:
        raised = "2^25" in String(e)
    assert_true(raised)


def _cigar(reference_row: String, row: String) -> String:
    """The pair's alignment as M/I/D runs, from two MSA rows (columns gapped in both dropped)."""
    var ops = List[UInt8]()
    var a = reference_row.as_bytes()
    var b = row.as_bytes()
    for c in range(len(a)):
        var gap_a = a[c] == UInt8(ord("-"))
        var gap_b = b[c] == UInt8(ord("-"))
        if gap_a and gap_b:
            continue
        ops.append(UInt8(ord("I")) if gap_a else (UInt8(ord("D")) if gap_b else UInt8(ord("M"))))
    var out = String()
    var k = 0
    while k < len(ops):
        var e = k
        while e < len(ops) and ops[e] == ops[k]:
            e += 1
        out.write(e - k, chr(Int(ops[k])))
        k = e
    return out^


def test_pinned_gaps_shift_left_as_minimap2_writes_them() raises:
    """chr20:60,425,968's two alleles pinned between their flanks: minimap2 2.31 (-x asm5, dipcall's settings)
    writes the 78-letter deletion at the reference's letter 517 (2017M78D1965M over flank + sequence + flank, flanks
    1,500), shifted left through the matched run before it past the pin at 586 (mm_fix_cigar); aligned between the
    pins alone it sat at 586."""
    var given = read_fasta("tests/data/chr20_left_shift/alleles.fa")
    var flanking = read_fasta("tests/data/chr20_left_shift/flanks.fa")
    var seqs = List[List[UInt8]]()
    for r in range(len(given)):
        seqs.append(List[UInt8](given.seq(r)))
    var flanks = List[List[UInt8]]()
    for r in range(len(flanking)):
        flanks.append(List[UInt8](flanking.seq(r)))
    var rows = on_reference(seqs, Costs.two_piece(40, 78, 7, 162, 3), 1, None, flanks)
    assert_equal(_cigar(String(unsafe_from_utf8=rows[0]), String(unsafe_from_utf8=rows[1])), "517M78D465M")


def test_left_aligned_shifts_each_gap_through_the_run_before_it() raises:
    """mm_fix_cigar's first loop on a pairing: a deletion or an insertion between two matched runs moves left while
    the letters allow; one at the sequence's end only with a right flank after it; none before the first match."""
    var reference = _b("GATTTTC")
    # GATTT_C against GATTTTC: the deletion of the last T moves to the first.
    var deleted: List[Int] = [0, 1, 2, 3, 4, 6]
    assert_true(_left_aligned(deleted, reference, _b("GATTTC"), True) == [0, 1, 3, 4, 5, 6])
    # An inserted T after the run moves to the run's start.
    var query = _b("GATTTTTC")
    var inserted: List[Int] = [0, 1, 2, 3, 4, 5, -1, 6]
    assert_true(_left_aligned(inserted, reference, query, True) == [0, 1, -1, 2, 3, 4, 5, 6])
    # A deletion at the end: shifted only when the right flank follows.
    var tail: List[Int] = [0, 1, 2, 3, 4, 5]
    assert_true(_left_aligned(tail, reference, _b("GATTTT"), False) == [0, 1, 2, 3, 4, 5])
    var shifted: List[Int] = [0, 1, 2, 3, 4]
    assert_true(_left_aligned(shifted, _b("GATTTT"), _b("GATTT"), True) == [0, 1, 3, 4, 5])


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
