"""Unit tests for `plegma/seeds.mojo`, `plegma/chain.mojo`, and `plegma/pins.mojo`."""

from std.testing import TestSuite, assert_equal, assert_true

from karan import Aligner, Costs

from plegma.chain import primary_chain, rmq_chains
from plegma.pins import pinned_pairing, pinned_pairing_with_index, pins, pins_with_index
from plegma.seeds import ReferenceIndex, collect_anchors, radix_sort, sketch


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
        var r = Int((s >> 33) % 100)
        if r == 0:
            out.append(ALPH.as_bytes()[Int((s >> 17) & 3)])
        elif r == 1:
            continue
        elif r == 2:
            out.append(seq[i])
            out.append(ALPH.as_bytes()[Int((s >> 17) & 3)])
        else:
            out.append(seq[i])
    return out^


def test_radix_sort_small_and_large() raises:
    # Small path (<= 32: insertion sort).
    var small = List[UInt64]()
    for i in range(19):
        small.append(UInt64((19 - i) * 7 + 3))
    radix_sort(small)
    for i in range(1, len(small)):
        assert_true(small[i - 1] <= small[i])

    # Large path (> 32: multi-digit LSD radix sort with sparse high bits).
    var large = List[UInt64]()
    var s = UInt64(1234567)
    for _ in range(257):
        s = s * 6364136223846793005 + 1442695040888963407
        large.append(s)
    var expected = large.copy()
    sort(expected)
    radix_sort(large)
    assert_equal(len(large), len(expected))
    for i in range(len(large)):
        assert_equal(large[i], expected[i])


def test_reference_index_matches_unindexed_anchors_and_pins() raises:
    var ref_seq = _rng_seq(101, 3000)
    var query = _mutate(ref_seq, 202)
    var index = ReferenceIndex(ref_seq)

    var a_direct = collect_anchors(ref_seq, query)
    var a_indexed = collect_anchors(index, query)
    assert_equal(len(a_direct), len(a_indexed))
    for i in range(len(a_direct)):
        assert_equal(a_direct[i].x, a_indexed[i].x)
        assert_equal(a_direct[i].q, a_indexed[i].q)
        assert_equal(a_direct[i].span, a_indexed[i].span)
        assert_equal(a_direct[i].flags, a_indexed[i].flags)

    var c_direct = primary_chain(ref_seq, query)
    var c_indexed = primary_chain(index, query)
    assert_true(len(c_direct) > 20)
    assert_equal(len(c_direct), len(c_indexed))

    var p_direct = pins(ref_seq, query)
    var p_indexed = pins_with_index(index, query)
    assert_true(len(p_direct.r) > 5)
    assert_equal(len(p_direct.r), len(p_indexed.r))
    for i in range(len(p_direct.r)):
        assert_equal(p_direct.r[i], p_indexed.r[i])
        assert_equal(p_direct.q[i], p_indexed.q[i])


def test_sketch_and_rmq_chains_basic() raises:
    var seq = _rng_seq(77, 1500)
    var sk = sketch(seq)
    assert_true(len(sk.x) > 50)
    assert_equal(len(sk.x), len(sk.y))
    var anchors = collect_anchors(seq, seq)
    var chains = rmq_chains(anchors, 2000)
    assert_true(len(chains) >= 1)
    assert_true(len(chains[0]) > 20)
    for i in range(1, len(chains[0])):
        assert_true(chains[0][i].x > chains[0][i - 1].x)
        assert_true(chains[0][i].q > chains[0][i - 1].q)


def test_pinned_pairing_parallel_matches_serial_and_indexed() raises:
    var left = _rng_seq(11, 400)
    var right = _rng_seq(22, 400)
    var ref_seq = _rng_seq(33, 2500)
    var query = _mutate(ref_seq, 44)
    var costs = Costs.two_piece(40, 78, 7, 162, 3)

    var whole_ref = List[UInt8]()
    for x in left:
        whole_ref.append(x)
    for x in ref_seq:
        whole_ref.append(x)
    for x in right:
        whole_ref.append(x)
    var idx = ReferenceIndex(whole_ref)

    var aligner = Aligner(warm=False)
    var p1 = pinned_pairing(aligner, ref_seq, query, left, right, costs, threads=1)
    var p4 = pinned_pairing(aligner, ref_seq, query, left, right, costs, threads=4)
    var pi4 = pinned_pairing_with_index(aligner, idx, ref_seq, query, left, right, costs, threads=4)
    assert_equal(len(p1), len(query))
    assert_equal(len(p4), len(p1))
    assert_equal(len(pi4), len(p1))
    for i in range(len(p1)):
        assert_equal(p4[i], p1[i])
        assert_equal(pi4[i], p1[i])


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
