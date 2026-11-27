"""Unit tests for `plegma/bands.mojo`, `plegma/placement.mojo`, and `plegma/runs.mojo`."""

from std.math import log10
from std.testing import TestSuite, assert_equal, assert_false, assert_true

from karan import Aligner, Costs, Mode, Ties

from plegma.bands import REPEAT_DISTANCE, REPEAT_K, RepeatMask, _Paths, _Realigned, _periods, _segments
from plegma.homopolymer import HomopolymerModel
from plegma.kernel import BOTH, LEFT, RIGHT, Scoring
from plegma.placement import _place_one, alignment_costs
from plegma.runs import _align_group_reads, _length_in, _set_runs, _two_lengths, _two_lengths_from


def _rng_seq(seed: UInt64, length: Int) -> List[UInt8]:
    comptime ALPH = "ACGT"
    var s = seed
    var out = List[UInt8](capacity=length)
    for _ in range(length):
        s = s * 6364136223846793005 + 1442695040888963407
        out.append(ALPH.as_bytes()[Int((s >> 33) & 3)])
    return out^


def _toy_model() raises -> HomopolymerModel:
    var csv = String("base,strand,left,right,n,m,count,probability,logProbabilityDb\n")
    var bases = "ACGT".as_bytes()
    for b in range(4):
        var ch = String(unsafe_from_utf8=Span(bases)[b : b + 1])
        for strand in range(2):
            for n in range(2, 10):
                for m in range(0, 12):
                    var p = 0.85 if m == n else (0.06 if abs(m - n) == 1 else 0.003)
                    csv += String(
                        ch,
                        ",",
                        strand,
                        ",*,*,",
                        n,
                        ",",
                        m,
                        ",",
                        100 if m == n else 1,
                        ",",
                        p,
                        ",",
                        10.0 * log10(p),
                        "\n",
                    )
    return HomopolymerModel.parse(csv)


def test_repeat_mask_segments_and_periods() raises:
    var left = _rng_seq(11, 160)
    var unit = List[UInt8]("ACGTAG".as_bytes())
    var right = _rng_seq(12, 160)
    var seq = List[UInt8]()
    seq.extend(left.copy())
    for _ in range(15):
        seq.extend(unit.copy())
    seq.extend(right.copy())

    var mask = RepeatMask(Span(seq), REPEAT_K, REPEAT_DISTANCE)
    assert_false(mask.covers(10, 50))
    assert_true(mask.covers(160, 160 + 90))
    assert_false(mask.covers(len(seq) - 50, len(seq) - 10))

    # `_segments` must never place a boundary inside the tandem repeat array.
    var discord = List[Int32](length=len(seq) + 1, fill=0)
    var bounds = _segments(seq, discord, mask)
    assert_equal(bounds[0], 0)
    assert_equal(bounds[len(bounds) - 1], len(seq))
    assert_true(len(bounds) >= 3)
    for k in range(1, len(bounds) - 1):
        assert_false(mask.covers(max(0, bounds[k] - REPEAT_K), min(len(seq), bounds[k] + REPEAT_K)))

    # `_periods` detects the 6-bp unit (and its 12-bp multiple) inside the tandem array, capped at 2.
    var periods = List[Int]()
    _periods(seq, 160 + 30, 160, periods)
    assert_equal(len(periods), 2)
    assert_equal(periods[0], 6)
    assert_equal(periods[1], 12)


def test_paths_fill_weights_anchors_and_parallel() raises:
    var consensus = _rng_seq(21, 300)
    var seqs = List[List[UInt8]]()
    # Distinct reads with BOTH, LEFT, RIGHT anchors, plus exact duplicate weighted reads.
    var r0 = consensus.copy()
    r0[50] = UInt8(ord("A")) if r0[50] != UInt8(ord("A")) else UInt8(ord("C"))
    var r1 = List[UInt8](Span(consensus)[0:180])
    var r2 = List[UInt8](Span(consensus)[120:])
    var r3 = consensus.copy()
    r3.insert(100, consensus[100])
    seqs.append(r0.copy())
    seqs.append(r1.copy())
    seqs.append(r2.copy())
    seqs.append(r3.copy())
    var anchors: List[Int] = [BOTH, LEFT, RIGHT, BOTH]
    var members: List[Int] = [0, 1, 2, 3]
    var weights: List[Int] = [3, 2, 2, 1]
    var same_as: List[Int] = [0, 1, 2, 3]

    var a1 = Aligner(warm=False)
    var a4 = Aligner(warm=False)
    var real1 = _Realigned(same_as.copy())
    var real4 = _Realigned(same_as.copy())
    var p1 = _Paths()
    var p4 = _Paths()
    p1.fill(a1, real1, consensus, seqs, members, anchors, weights, threads=1)
    p4.fill(a4, real4, consensus, seqs, members, anchors, weights, threads=4)

    assert_equal(p4.total, p1.total)
    assert_equal(p4.spanned, p1.spanned)
    assert_equal(p4.matched, p1.matched)
    assert_equal(p4.substituted, p1.substituted)
    assert_equal(p4.inserted, p1.inserted)
    assert_equal(p4.deleted, p1.deleted)
    assert_equal(p4.copying, p1.copying)
    assert_equal(p4.enter, p1.enter)
    assert_equal(p4.leave, p1.leave)
    assert_equal(p4.discord, p1.discord)
    assert_equal(p4.cover, p1.cover)
    assert_true(p1.substituted == 3)
    assert_true(p1.inserted == 1 and p1.copying == 1)


def test_placement_costs_and_place_one_margin() raises:
    var affine_scoring = Scoring(2, 4, 4, 2, 0, 0)
    var two_piece_scoring = Scoring(2, 4, 4, 2, 24, 1)
    var flat_second_scoring = Scoring(2, 4, 4, 2, 24, 0)
    var c_aff = alignment_costs(affine_scoring)
    var c_two = alignment_costs(two_piece_scoring)
    var c_flat = alignment_costs(flat_second_scoring)
    assert_equal(c_aff.mismatch, 12)
    assert_equal(c_aff.opening, 8)
    assert_equal(c_aff.extension, 6)
    assert_equal(c_two.opening2, 48)
    assert_equal(c_two.extension2, 4)
    assert_equal(c_flat.opening2, 48)
    assert_equal(c_flat.extension2, 2)

    var hap0 = _rng_seq(31, 120)
    var hap1 = hap0.copy()
    hap1[40] = UInt8(ord("A")) if hap1[40] != UInt8(ord("A")) else UInt8(ord("T"))
    hap1[80] = UInt8(ord("A")) if hap1[80] != UInt8(ord("A")) else UInt8(ord("T"))
    var made = List[List[UInt8]]()
    made.append(hap0.copy())
    made.append(hap1.copy())

    var aligner = Aligner(warm=False)
    # Exact match to hap0 -> cluster 0 with positive margin == 2 mismatches (24).
    var p0 = _place_one(aligner, made, hap0, -1, BOTH, c_two)
    assert_equal(p0.cluster, 0)
    assert_equal(p0.cost, 0)
    assert_equal(p0.margin, 24)
    # Forced cluster -> margin is 0 even when another cluster differs.
    var forced = _place_one(aligner, made, hap0, 1, BOTH, c_two)
    assert_equal(forced.cluster, 1)
    assert_equal(forced.margin, 0)
    # Equidistant read (1 mismatch to each) -> margin is 0 and tie picks cluster 0.
    var mid = hap0.copy()
    mid[40] = hap1[40]
    var tied = _place_one(aligner, made, mid, -1, BOTH, c_two)
    assert_equal(tied.cluster, 0)
    assert_equal(tied.margin, 0)


def test_runs_align_group_length_in_set_runs_and_two_lengths() raises:
    var prefix = List[UInt8]("ACGTACGTAC".as_bytes())
    var suffix = List[UInt8]("GCTAGCTAGC".as_bytes())
    var hap0 = List[UInt8]()
    hap0.extend(prefix.copy())
    hap0.extend(List[UInt8]("TTTT".as_bytes()))
    hap0.extend(suffix.copy())
    var hap1 = List[UInt8]()
    hap1.extend(prefix.copy())
    hap1.extend(List[UInt8]("TTTTTTT".as_bytes()))
    hap1.extend(suffix.copy())

    # `_set_runs` updates run lengths cleanly.
    var updated = _set_runs(hap0, [10], [7])
    assert_equal(updated, hap1)

    # `_length_in` reports exact length, 0 for fully deleted run, and -1 when corrupted or at boundary.
    var aligner = Aligner(warm=False)
    var a_same = aligner.align(Span(hap0), Span(hap0), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
    assert_equal(_length_in(a_same.cigar, a_same.reference_start, hap0, 10, 14, UInt8(ord("T")), hap0), 4)
    var del_read = List[UInt8]()
    del_read.extend(prefix.copy())
    del_read.extend(suffix.copy())
    var a_del = aligner.align(Span(hap0), Span(del_read), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
    assert_equal(_length_in(a_del.cigar, a_del.reference_start, del_read, 10, 14, UInt8(ord("T")), hap0), 0)
    var bad_read = hap0.copy()
    bad_read[12] = UInt8(ord("G"))
    var a_bad = aligner.align(Span(hap0), Span(bad_read), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
    assert_equal(_length_in(a_bad.cigar, a_bad.reference_start, bad_read, 10, 14, UInt8(ord("T")), hap0), -1)
    assert_equal(_length_in(a_same.cigar, a_same.reference_start, hap0, 0, 1, hap0[0], hap0), -1)

    # `_align_group_reads` 1T == 4T and `_two_lengths_from` matches `_two_lengths` with strong separation.
    var seqs = List[List[UInt8]]()
    var anchors = List[Int]()
    var strands = List[Int]()
    var g0 = List[Int]()
    var g1 = List[Int]()
    for k in range(6):
        seqs.append(hap0.copy())
        anchors.append(BOTH)
        strands.append(k % 2)
        g0.append(k)
    for k in range(6, 12):
        seqs.append(hap1.copy())
        anchors.append(BOTH)
        strands.append(k % 2)
        g1.append(k)

    var al0_1 = _align_group_reads(aligner, hap0, seqs, g0, anchors, threads=1)
    var al0_4 = _align_group_reads(aligner, hap0, seqs, g0, anchors, threads=4)
    assert_equal(len(al0_1), len(al0_4))
    for k in range(len(al0_1)):
        assert_equal(al0_4[k].cost, al0_1[k].cost)
        assert_equal(al0_4[k].cigar, al0_1[k].cigar)
    var al1_4 = _align_group_reads(aligner, hap1, seqs, g1, anchors, threads=4)

    var model = _toy_model()
    var made = List[List[UInt8]]()
    made.append(hap0.copy())
    made.append(hap1.copy())
    var groups = List[List[Int]]()
    groups.append(g0.copy())
    groups.append(g1.copy())
    var starts: List[Int] = [10, 10]
    var s1 = 0
    var l1 = 0
    var one1 = -1
    var said1 = List[List[Int]]()
    var gain_cached = _two_lengths_from(made, starts, groups, al0_4, al1_4, seqs, strands, model, s1, l1, one1, said1)
    var s2 = 0
    var l2 = 0
    var one2 = -1
    var said2 = List[List[Int]]()
    var gain_direct = _two_lengths(
        aligner, made, starts, groups, seqs, anchors, strands, model, s2, l2, one2, said2, threads=1
    )
    assert_equal(s1, 4)
    assert_equal(l1, 7)
    assert_equal(s1, s2)
    assert_equal(l1, l2)
    assert_equal(one1, one2)
    assert_equal(gain_cached, gain_direct)
    assert_true(gain_cached > 50.0)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
