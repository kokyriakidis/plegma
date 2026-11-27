"""Tests for `ConsensusScratch` reuse and bit-exact parallel determinism across every major plegma API."""

from std.testing import TestSuite, assert_equal, assert_raises, assert_true

from karan import Costs

from plegma.clustered import clustered_consensus, clustered_consensus_batch
from plegma.consensus import ConsensusOptions, ConsensusScratch, consensus, consensus_batch
from plegma.graph import Graph
from plegma.homopolymer import FLANK_UNKNOWN, STRAND_UNKNOWN, HomopolymerModel
from plegma.kernel import BOTH, INSIDE, LEFT, RIGHT, Scoring, Topology, Workspace
from plegma.phased import phased_consensus, phased_consensus_batch
from plegma.pinned import WindowScratch, align_pinned
from plegma.poa import LONGEST_FIRST, PoaOptions, _longest_first, poa
from plegma.polish import PolishScratch, polish, polish_with_scratch
from plegma.reference import on_reference
from plegma.runs import _best_length


def _rng_seq(seed: UInt64, length: Int) -> List[UInt8]:
    comptime ALPH = "ACGT"
    var s = seed
    var out = List[UInt8](capacity=length)
    for _ in range(length):
        s = s * 6364136223846793005 + 1442695040888963407
        out.append(ALPH.as_bytes()[Int((s >> 33) & 3)])
    return out^


def _mutate(seq: List[UInt8], seed: UInt64, rate_pct: Int = 3) -> List[UInt8]:
    comptime ALPH = "ACGT"
    var s = seed
    var out = List[UInt8](capacity=len(seq))
    for i in range(len(seq)):
        s = s * 6364136223846793005 + 1442695040888963407
        var r = Int((s >> 33) % 100)
        if r < rate_pct:
            var kind = Int((s >> 21) % 3)
            if kind == 0:
                out.append(ALPH.as_bytes()[Int((s >> 17) & 3)])
            elif kind == 1:
                continue
            else:
                out.append(seq[i])
                out.append(ALPH.as_bytes()[Int((s >> 17) & 3)])
        else:
            out.append(seq[i])
    return out^


def test_consensus_scratch_reuse_matches_fresh_calls() raises:
    var scratch = ConsensusScratch()
    var opts = ConsensusOptions()
    opts.poa.seeded = False

    for c in range(4):
        var truth = _rng_seq(UInt64(100 + c), 180 + c * 90)
        var reads = List[List[UInt8]]()
        for r in range(12):
            reads.append(_mutate(truth, UInt64(1000 * c + r), 3))
        var fresh = consensus(reads.copy(), options=opts)
        var reused = consensus(scratch, reads.copy(), options=opts)
        assert_equal(len(fresh), len(reused))
        for p in range(len(fresh)):
            assert_equal(reused[p], fresh[p])


def test_align_pinned_parallel_matches_serial() raises:
    var truth = _rng_seq(202, 2400)
    var graph = Graph()
    _ = graph.add("r0", Span(truth))
    var q = _mutate(truth, 303, 2)
    var topo = Topology(graph)
    var scoring = Scoring()
    var row1 = List[Int](length=len(graph), fill=-1)
    var row4 = List[Int](length=len(graph), fill=-1)
    var ws1 = Workspace()
    var ws4 = Workspace()
    var sc1 = WindowScratch()
    var sc4 = WindowScratch()
    var a1 = align_pinned(graph, topo, q, scoring, 10, 0.01, row1, sc1, ws1, threads=1)
    var a4 = align_pinned(graph, topo, q, scoring, 10, 0.01, row4, sc4, ws4, threads=4)
    assert_equal(a4.score, a1.score)
    assert_equal(a4.cells, a1.cells)
    assert_equal(a4.nodes, a1.nodes)


def test_polish_parallel_single_and_multi_segment_bit_exact() raises:
    # Single-segment case (short sequence < 2*SEGMENT) tests read-level parallelism in both edit & likelihood stages.
    var short_truth = _rng_seq(401, 90)
    var short_reads = List[List[UInt8]]()
    var members = List[Int]()
    for r in range(20):
        short_reads.append(_mutate(short_truth, UInt64(4000 + r), 4))
        members.append(r)
    var rough_short = _mutate(short_truth, 4999, 6)
    var p1_s = polish(rough_short, short_reads, members, 20, likely_rounds=10, threads=1)
    var p4_s = polish(rough_short, short_reads, members, 20, likely_rounds=10, threads=4)
    assert_equal(p4_s, p1_s)

    # Multi-segment case (~1200 bp) tests segment-level parallelism in both edit & likelihood stages.
    var long_truth = _rng_seq(501, 1200)
    var long_reads = List[List[UInt8]]()
    for r in range(16):
        long_reads.append(_mutate(long_truth, UInt64(5000 + r), 3))
    var lmembers = List[Int]()
    for r in range(16):
        lmembers.append(r)
    var rough_long = _mutate(long_truth, 5999, 4)
    var p1_l = polish(rough_long, long_reads, lmembers, 20, likely_rounds=10, threads=1)
    var p4_l = polish(rough_long, long_reads, lmembers, 20, likely_rounds=10, threads=4)
    assert_equal(p4_l, p1_l)


def test_poa_and_consensus_parallel_bit_exact() raises:
    # Two-sided assembly with no spanning read: exercises concurrent left/right _one_side, parallel join_sides, polish.
    var truth = _rng_seq(601, 1500)
    var reads = List[List[UInt8]]()
    var anchors = List[Int]()
    for r in range(10):
        var sl = List[UInt8]()
        for i in range(0, 950):
            sl.append(truth[i])
        reads.append(_mutate(sl, UInt64(6000 + r), 2))
        anchors.append(LEFT)
    for r in range(10):
        var sr = List[UInt8]()
        for i in range(550, 1500):
            sr.append(truth[i])
        reads.append(_mutate(sr, UInt64(6100 + r), 2))
        anchors.append(RIGHT)

    var o1 = ConsensusOptions()
    o1.poa.threads = 1
    var o4 = ConsensusOptions()
    o4.poa.threads = 4
    var c1 = consensus(reads.copy(), anchors, options=o1)
    var c4 = consensus(reads.copy(), anchors, options=o4)
    assert_equal(len(c1), 1)
    assert_equal(len(c4), 1)
    assert_equal(c4[0], c1[0])


def test_clustered_and_phased_parallel_and_scratch_bit_exact() raises:
    var hap0 = _rng_seq(701, 500)
    var hap1 = hap0.copy()
    for pos in [80, 160, 240, 320, 400]:
        hap1[pos] = UInt8(ord("A")) if hap1[pos] != UInt8(ord("A")) else UInt8(ord("G"))
    var reads = List[List[UInt8]]()
    var labels = List[Int]()
    for r in range(24):
        var h = r % 2
        reads.append(_mutate(hap0 if h == 0 else hap1, UInt64(7000 + r), 2))
        labels.append(h if r < 16 else -1)

    var o1 = ConsensusOptions()
    o1.poa.seeded = False
    o1.poa.threads = 1
    var o4 = ConsensusOptions()
    o4.poa.seeded = False
    o4.poa.threads = 4

    var scratch = ConsensusScratch()
    var cl1 = clustered_consensus(scratch, reads.copy(), options=o1)
    var cl4 = clustered_consensus(reads.copy(), options=o4)
    assert_equal(len(cl1.consensus), 2)
    assert_equal(len(cl4.consensus), len(cl1.consensus))
    for k in range(len(cl1.consensus)):
        assert_equal(cl4.consensus[k], cl1.consensus[k])
        assert_equal(cl4.members[k], cl1.members[k])
    for r in range(len(reads)):
        assert_equal(cl4.placements[r].cluster, cl1.placements[r].cluster)
        assert_equal(cl4.placements[r].cost, cl1.placements[r].cost)
        assert_equal(cl4.placements[r].cigar, cl1.placements[r].cigar)
        assert_equal(cl4.placements[r].margin, cl1.placements[r].margin)

    var ph1 = phased_consensus(scratch, reads.copy(), labels, options=o1)
    var ph4 = phased_consensus(reads.copy(), labels, options=o4)
    assert_equal(len(ph1.consensus), 2)
    for k in range(2):
        assert_equal(ph4.consensus[k], ph1.consensus[k])
        assert_equal(ph4.members[k], ph1.members[k])
    for r in range(len(reads)):
        assert_equal(ph4.placements[r].cluster, ph1.placements[r].cluster)
        assert_equal(ph4.placements[r].cost, ph1.placements[r].cost)
        assert_equal(ph4.placements[r].cigar, ph1.placements[r].cigar)


def test_on_reference_parallel_fallback_and_insertions_bit_exact() raises:
    var ref_seq = _rng_seq(801, 400)
    var seqs = List[List[UInt8]]()
    seqs.append(ref_seq.copy())
    # Create multiple distinct sequences with multi-base insertions at several positions so _on_reference_columns
    # runs parallel insertion-slot sub-alignments.
    for k in range(12):
        var s = List[UInt8]()
        for i in range(len(ref_seq)):
            s.append(ref_seq[i])
            if i == 100:
                var ins = _rng_seq(UInt64(8100 + k), 12 + (k % 5))
                for b in ins:
                    s.append(b)
            if i == 250:
                var ins = _rng_seq(UInt64(8200 + k), 10 + (k % 4))
                for b in ins:
                    s.append(b)
        seqs.append(s^)
    var costs = Costs.two_piece(40, 78, 7, 162, 3)
    var fb = Costs.two_piece(24, 16, 12, 96, 5)
    var r1 = on_reference(seqs, costs, threads=1, fallback=fb)
    var r4 = on_reference(seqs, costs, threads=4, fallback=fb)
    assert_equal(len(r1), len(seqs))
    assert_equal(len(r4), len(r1))
    for i in range(len(r1)):
        assert_equal(r4[i], r1[i])


def test_speculative_seeding_and_snapshot_batch_poa_bit_exact() raises:
    # 30 reads of 1800 bp: exercises speculative parallel seeding (batch_size=1) and snapshot mini-batch POA
    # (batch_size=8) for both pinned and unpinned modes, verifying 1T == 4T bit-exact and accurate recovery.
    var truth = _rng_seq(901, 1800)
    var reads = List[List[UInt8]]()
    for r in range(30):
        reads.append(_mutate(truth, UInt64(9000 + r), 3))

    # Pinned strict sequential (batch_size=1): speculative seeding 4T must match 1T bit-for-bit.
    var seq1 = ConsensusOptions()
    seq1.poa.seeded = True
    seq1.poa.threads = 1
    var seq4 = ConsensusOptions()
    seq4.poa.seeded = True
    seq4.poa.threads = 4
    var c_seq1 = consensus(reads.copy(), options=seq1)
    var c_seq4 = consensus(reads.copy(), options=seq4)
    assert_equal(c_seq4[0], c_seq1[0])
    assert_equal(c_seq4[0], truth)

    # Snapshot mini-batch POA (batch_size=8): 1T == 4T bit-exact for BOTH pinned and unpinned modes.
    for seeded in [False, True]:
        var b1 = ConsensusOptions()
        b1.poa.seeded = seeded
        b1.poa.batch_size = 8
        b1.poa.threads = 1
        var b4 = ConsensusOptions()
        b4.poa.seeded = seeded
        b4.poa.batch_size = 8
        b4.poa.threads = 4
        var cb1 = consensus(reads.copy(), options=b1)
        var cb4 = consensus(reads.copy(), options=b4)
        assert_equal(cb4[0], cb1[0])
        assert_equal(cb4[0], truth)


def test_multi_locus_batch_apis_bit_exact() raises:
    var seqs_batch = List[List[List[UInt8]]]()
    var haps_batch = List[List[Int]]()
    for locus in range(6):
        var hap0 = _rng_seq(UInt64(1100 + locus), 300 + locus * 50)
        var hap1 = hap0.copy()
        for pos in [40, 100, 160, 220]:
            hap1[pos] = UInt8(ord("A")) if hap1[pos] != UInt8(ord("A")) else UInt8(ord("C"))
        var reads = List[List[UInt8]]()
        var haps = List[Int]()
        for r in range(16):
            var h = r % 2
            reads.append(_mutate(hap0 if h == 0 else hap1, UInt64(12000 + locus * 100 + r), 2))
            haps.append(h if r < 12 else -1)
        seqs_batch.append(reads^)
        haps_batch.append(haps^)

    var o1 = ConsensusOptions()
    o1.poa.seeded = False
    o1.poa.threads = 1
    var o4 = ConsensusOptions()
    o4.poa.seeded = False
    o4.poa.threads = 4

    var cb1 = consensus_batch(seqs_batch, options=o1)
    var cb4 = consensus_batch(seqs_batch, options=o4)
    assert_equal(len(cb4), len(cb1))
    for i in range(len(cb1)):
        assert_equal(cb4[i][0], cb1[i][0])

    var clb1 = clustered_consensus_batch(seqs_batch, options=o1)
    var clb4 = clustered_consensus_batch(seqs_batch, options=o4)
    assert_equal(len(clb4), len(clb1))
    for i in range(len(clb1)):
        assert_equal(len(clb4[i].consensus), len(clb1[i].consensus))
        for k in range(len(clb1[i].consensus)):
            assert_equal(clb4[i].consensus[k], clb1[i].consensus[k])
            assert_equal(clb4[i].members[k], clb1[i].members[k])

    var phb1 = phased_consensus_batch(seqs_batch, haps_batch, options=o1)
    var phb4 = phased_consensus_batch(seqs_batch, haps_batch, options=o4)
    assert_equal(len(phb4), len(phb1))
    for i in range(len(phb1)):
        for k in range(2):
            assert_equal(phb4[i].consensus[k], phb1[i].consensus[k])
            assert_equal(phb4[i].members[k], phb1[i].members[k])


def test_audit_fixes_and_flattened_polish_bit_exact() raises:
    # 1. `_longest_first` must never sort empty reads before non-empty reads.
    var mixed: List[List[UInt8]] = [List[UInt8](), List[UInt8]("ACGTACGT".as_bytes()), List[UInt8]("AC".as_bytes())]
    var lf = _longest_first(mixed, [0, 1, 2])
    assert_equal(lf[0], 1)
    assert_equal(lf[1], 2)
    assert_equal(lf[2], 0)

    # 2. Flattened `(active_seg, read_k)` edit + likelihood polishing with duplicate reads & 2..4 segments: 1T == 8T.
    var truth = _rng_seq(404, 420)
    var bad = truth.copy()
    bad[80] = UInt8(ord("A")) if bad[80] != UInt8(ord("A")) else UInt8(ord("T"))
    bad[240] = UInt8(ord("G")) if bad[240] != UInt8(ord("G")) else UInt8(ord("C"))
    bad[350] = UInt8(ord("C")) if bad[350] != UInt8(ord("C")) else UInt8(ord("A"))
    var reads = List[List[UInt8]]()
    var members = List[Int]()
    for r in range(24):
        var seq = truth.copy() if r % 3 == 0 else _mutate(truth, UInt64(5000 + r // 2), 2)
        reads.append(seq^)
        members.append(r)
    var p1 = polish(bad, reads, members, 20, likely_rounds=10, threads=1)
    var p8 = polish(bad, reads, members, 20, likely_rounds=10, threads=8)
    assert_equal(p8, p1)
    assert_equal(p1, truth)

    # 3. `HomopolymerModel.most_likely` with `STRAND_UNKNOWN` uses precomputed SIMD slice & matches bounds safety.
    var csv = String("base,strand,left,right,n,m,count,probability,logProbabilityDb\n")
    for base_ch in ["A", "C", "G", "T"]:
        for strand in range(2):
            for n in range(2, 5):
                for m in range(5):
                    var lp = -0.5 if m == n else -15.0
                    csv += String(base_ch, ",", strand, ",*,*,", n, ",", m, ",100,0.5,", lp, "\n")
    var model = HomopolymerModel.parse(csv)
    assert_equal(
        model.most_likely(0, FLANK_UNKNOWN, FLANK_UNKNOWN, [3, 3, 3, 4], [STRAND_UNKNOWN, 0, 1, STRAND_UNKNOWN]), 3
    )
    assert_equal(model.most_likely(4, FLANK_UNKNOWN, FLANK_UNKNOWN, [3], [0]), -1)
    # `_best_length` with 1 valid read on uncovered `n` must return -1 (not accept -Float64.MAX tie).
    var sc = 0.0
    assert_equal(_best_length(model, 0, FLANK_UNKNOWN, FLANK_UNKNOWN, [3], [0], [0], 10, 12, sc), -1)

    # 4. `clustered_consensus` & `phased_consensus` validate options even on empty input, and odd threads match 1T.
    var bad_opts = ConsensusOptions()
    bad_opts.poa.threads = 0
    with assert_raises():
        _ = clustered_consensus(List[List[UInt8]](), options=bad_opts)
    with assert_raises():
        _ = phased_consensus(List[List[UInt8]](), List[Int](), options=bad_opts)
    var o3 = ConsensusOptions()
    o3.poa.seeded = False
    o3.poa.threads = 3
    var o1 = ConsensusOptions()
    o1.poa.seeded = False
    o1.poa.threads = 1
    var haps = List[Int](capacity=len(reads))
    for r in range(len(reads)):
        haps.append(0 if r < 16 else (1 if r < 20 else -1))
    var ph1 = phased_consensus(reads.copy(), haps, options=o1)
    var ph3 = phased_consensus(reads.copy(), haps, options=o3)
    assert_equal(ph3.consensus[0], ph1.consensus[0])
    assert_equal(ph3.consensus[1], ph1.consensus[1])
    assert_equal(ph3.members[0], ph1.members[0])
    assert_equal(ph3.members[1], ph1.members[1])


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
