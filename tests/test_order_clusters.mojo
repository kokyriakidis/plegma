"""Unit tests for `plegma/order.mojo`, `plegma/clusters.mojo`, and `plegma/bands.mojo`."""

from std.testing import TestSuite, assert_equal, assert_true

from plegma.bands import REPEAT_DISTANCE, REPEAT_K, RepeatMask, _periods, _segments
from plegma.clusters import cluster_rows
from plegma.order import guide_order


def _rng_seq(seed: UInt64, length: Int) -> List[UInt8]:
    comptime ALPH = "ACGT"
    var s = seed
    var out = List[UInt8](capacity=length)
    for _ in range(length):
        s = s * 6364136223846793005 + 1442695040888963407
        out.append(ALPH.as_bytes()[Int((s >> 33) & 3)])
    return out^


def _mutate(seq: List[UInt8], seed: UInt64) -> List[UInt8]:
    comptime ALPH = "ACGT"
    var s = seed
    var out = List[UInt8](capacity=len(seq))
    for i in range(len(seq)):
        s = s * 6364136223846793005 + 1442695040888963407
        var r = Int((s >> 33) % 40)
        if r == 0:
            out.append(ALPH.as_bytes()[Int((s >> 17) & 3)])
        else:
            out.append(seq[i])
    return out^


def test_guide_order_parallel_matches_serial_and_is_permutation() raises:
    var base = _rng_seq(17, 600)
    var seqs = List[List[UInt8]]()
    for i in range(28):
        seqs.append(_mutate(base, UInt64(1000 + i)))
    var o1 = guide_order(seqs, threads=1)
    var o4 = guide_order(seqs, threads=4)
    assert_equal(len(o1), len(seqs))
    assert_equal(o4, o1)
    var seen = List[Bool](length=len(seqs), fill=False)
    for x in o1:
        assert_true(not seen[x])
        seen[x] = True

    # Duplicate-sequence tie-breaking must not depend on input order.
    var s_dup = _rng_seq(99, 200)
    var s_a = _rng_seq(100, 200)
    var s_b = _rng_seq(101, 200)
    var fwd: List[List[UInt8]] = [s_dup.copy(), s_b.copy(), s_dup.copy(), s_a.copy()]
    var rev: List[List[UInt8]] = [s_a.copy(), s_dup.copy(), s_b.copy(), s_dup.copy()]
    var ofwd = guide_order(fwd, threads=1)
    var orev = guide_order(rev, threads=1)
    for i in range(4):
        assert_equal(fwd[ofwd[i]], rev[orev[i]])


def test_cluster_rows_parallel_matches_serial() raises:
    var rows = List[String]()
    for i in range(48):
        if i % 2 == 0:
            rows.append("ACGTACGTACGTACGTACGT")
        else:
            rows.append("ACGTTCGATCGAACGTACGT")
    var c1 = cluster_rows(rows, 2, 0.25, threads=1)
    var c4 = cluster_rows(rows, 2, 0.25, threads=4)
    assert_equal(len(c1), 2)
    assert_equal(len(c4), 2)
    assert_equal(c4[0], c1[0])
    assert_equal(c4[1], c1[1])


def test_repeat_mask_segments_and_periods() raises:
    var unit = _rng_seq(55, 15)
    var seq = List[UInt8]()
    for b in _rng_seq(56, 120):
        seq.append(b)
    var rep_start = len(seq)
    for _ in range(8):
        for b in unit:
            seq.append(b)
    var rep_end = len(seq)
    for b in _rng_seq(57, 120):
        seq.append(b)
    var mask = RepeatMask(Span(seq), REPEAT_K, REPEAT_DISTANCE)
    assert_true(mask.covers(rep_start, rep_end))
    assert_true(not mask.covers(0, 80))

    var discord = List[Int32](length=len(seq) + 1, fill=0)
    var bounds = _segments(seq, discord, mask)
    assert_equal(bounds[0], 0)
    assert_equal(bounds[len(bounds) - 1], len(seq))
    # No segment boundary placed deep inside the tandem repeat.
    for k in range(1, len(bounds) - 1):
        assert_true(bounds[k] <= rep_start or bounds[k] >= rep_end)

    var periods = List[Int]()
    _periods(seq, rep_start + 30, 0, periods)
    assert_true(len(periods) >= 1)
    assert_equal(periods[0], 15)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
