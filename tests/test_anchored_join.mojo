"""Unit tests for `plegma/anchored.mojo` and `plegma/join.mojo` (serial and parallel)."""

from std.testing import TestSuite, assert_equal, assert_true

from karan import Aligner

from plegma.anchored import _STRAY, ends_checked, find_foreign, trim_consensus
from plegma.join import join_sides, repeat_period
from plegma.kernel import BOTH, LEFT, RIGHT


def _rng_seq(seed: UInt64, length: Int) -> List[UInt8]:
    comptime ALPH = "ACGT"
    var s = seed
    var out = List[UInt8](capacity=length)
    for _ in range(length):
        s = s * 6364136223846793005 + 1442695040888963407
        out.append(ALPH.as_bytes()[Int((s >> 33) & 3)])
    return out^


def _slice(seq: List[UInt8], start: Int, end: Int) -> List[UInt8]:
    var out = List[UInt8](capacity=end - start)
    for i in range(start, end):
        out.append(seq[i])
    return out^


def test_repeat_period_detects_tandem_and_rejects_random() raises:
    var unit = _rng_seq(9, 27)
    var array = List[UInt8]()
    for _ in range(15):
        for b in unit:
            array.append(b)
    var p = repeat_period(array, 0, len(array))
    assert_true(p >= 26 and p <= 28)
    var rand = _rng_seq(19, 400)
    assert_equal(repeat_period(rand, 0, len(rand)), 0)


def test_join_sides_serial_matches_parallel_and_recovers_truth() raises:
    var truth = _rng_seq(31, 2200)
    var left = _slice(truth, 0, 1400)
    var right = _slice(truth, 800, 2200)
    var left_reach: List[Int] = [1400, 1300, 1200, 1100, 1000]
    var right_reach: List[Int] = [1400, 1300, 1200, 1100, 1000]

    var aligner = Aligner(warm=False)
    var j1 = join_sides(aligner, left, right, left_reach, right_reach, threads=1)
    var j4 = join_sides(aligner, left, right, left_reach, right_reach, threads=4)
    assert_equal(j1, truth)
    assert_equal(j4, j1)


def test_join_sides_tandem_copy_shifts_parallel_matches_serial() raises:
    var flank_l = _rng_seq(501, 400)
    var flank_r = _rng_seq(502, 400)
    var unit = _rng_seq(503, 30)
    # Build a tandem array with distinguishable mutated copies so the true copy number wins uniquely.
    var truth = List[UInt8]()
    for b in flank_l:
        truth.append(b)
    for c in range(25):
        for t in range(len(unit)):
            var ch = unit[t]
            if t == c % len(unit):
                ch = UInt8(ord("A")) if ch != UInt8(ord("A")) else UInt8(ord("C"))
            truth.append(ch)
    for b in flank_r:
        truth.append(b)
    var left = _slice(truth, 0, 1050)
    var right = _slice(truth, 450, len(truth))
    var reach: List[Int] = [1000, 950, 900, 850, 800]

    var aligner = Aligner(warm=False)
    var j1 = join_sides(aligner, left, right, reach, reach, threads=1)
    var j4 = join_sides(aligner, left, right, reach, reach, threads=4)
    assert_equal(j1, truth)
    assert_equal(j4, j1)


def test_ends_checked_and_find_foreign_parallel_matches_serial() raises:
    var truth = _rng_seq(701, 600)
    var foreign_seq = _rng_seq(702, 600)
    var seqs = List[List[UInt8]]()
    var anchors = List[Int]()
    # 0: spanning ending short (should become LEFT); 1: stray (should become _STRAY); 2..5: left; 6..9: right.
    seqs.append(_slice(truth, 0, 350))
    anchors.append(BOTH)
    seqs.append(foreign_seq.copy())
    anchors.append(BOTH)
    for _ in range(4):
        seqs.append(_slice(truth, 0, 300))
        anchors.append(LEFT)
    for _ in range(4):
        seqs.append(_slice(truth, 300, 600))
        anchors.append(RIGHT)

    var aligner = Aligner(warm=False)
    var c1 = ends_checked(aligner, seqs, anchors, threads=1)
    var c4 = ends_checked(aligner, seqs, anchors, threads=4)
    assert_equal(c1[0], LEFT)
    assert_equal(c1[1], _STRAY)
    assert_equal(c4, c1)

    # Foreign detection: 3 spanning matching truth, 2 true left, 2 foreign right.
    var fseqs = List[List[UInt8]]()
    for _ in range(3):
        fseqs.append(truth.copy())
    fseqs.append(_slice(truth, 0, 300))
    fseqs.append(_slice(truth, 0, 280))
    fseqs.append(_slice(foreign_seq, 0, 300))
    fseqs.append(_slice(foreign_seq, 200, 600))
    var sp: List[Int] = [0, 1, 2]
    var lf: List[Int] = [3, 4]
    var rt: List[Int] = [5, 6]
    var f1 = find_foreign(aligner, truth, fseqs, sp, lf, rt, threads=1)
    var f4 = find_foreign(aligner, truth, fseqs, sp, lf, rt, threads=4)
    assert_equal(f1, [5, 6])
    assert_equal(f4, f1)


def test_trim_consensus_left_and_right() raises:
    var left_seq = _rng_seq(81, 100)
    var reach: List[Int] = [95, 80, 70, 60, 40]
    var expected_left = _slice(left_seq, 0, 70)
    trim_consensus(left_seq, reach, 3, True)
    assert_equal(left_seq, expected_left)

    var right_seq = _rng_seq(82, 100)
    var expected_right = _slice(right_seq, 30, 100)
    trim_consensus(right_seq, reach, 3, False)
    assert_equal(right_seq, expected_right)
    assert_true(len(right_seq) == 70)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
