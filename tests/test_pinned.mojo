from std.random import random_ui64, seed
from std.testing import TestSuite, assert_equal, assert_false, assert_true

from karan import distance

from plegma.graph import Graph
from plegma.kernel import Scoring, Topology, Workspace
from plegma.pinned import UNIQUE_K, WindowScratch, align_pinned, unique_pins
from plegma.pins import Pins


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


def _ungapped(row: String) -> String:
    var out = List[UInt8]()
    for c in row.as_bytes():
        if c != UInt8(ord("-")):
            out.append(c)
    return String(unsafe_from_utf8=out)


def test_pinned_alignment_spells_every_sequence_and_keeps_the_consensus() raises:
    """Noisy copies of one sequence, a copy-number change among them, each pinned to the consensus so far: every
    row spells its sequence, the graph stays acyclic, and the consensus is the sequence's."""
    seed(83)
    var unit = _random(37)
    var truth = _random(1500)
    for _ in range(12):
        truth.extend(unit.copy())
    truth.extend(_random(1500))
    var seqs = List[List[UInt8]]()
    for k in range(6):
        seqs.append(_mutated(truth, 0.02))
    # One copy-number change: four copies fewer.
    var shorter = truth.copy()
    shorter.resize(1500 + 8 * 37, 0)
    var tail = List[UInt8]()
    for i in range(1500 + 12 * 37, len(truth)):
        tail.append(truth[i])
    shorter.extend(tail^)
    seqs.append(_mutated(shorter, 0.02))
    var graph = Graph()
    var row_of = List[Int]()
    var scratch = WindowScratch()
    var workspace = Workspace()
    _ = graph.add("s0", Span(seqs[0]))
    for k in range(1, len(seqs)):
        var topology = Topology(graph)
        row_of.resize(len(graph), -1)
        var found = align_pinned(graph, topology, seqs[k], Scoring(), 10, 0.01, row_of, scratch, workspace)
        _ = graph.fuse_nodes(String(t"s{k}"), Span(seqs[k]), found.nodes)
    _ = graph.node_order()
    var rows = graph.msa()
    for k in range(len(seqs)):
        assert_equal(rows[k].byte_length(), rows[0].byte_length())
        assert_true(_ungapped(rows[k]) == String(unsafe_from_utf8=seqs[k].copy()))
    var letters = List[UInt8]()
    for node in graph.heaviest_path(Topology(graph).order):
        letters.append(graph.bases[node])
    assert_true(distance(String(unsafe_from_utf8=letters), String(unsafe_from_utf8=truth)) <= 20)


def test_unique_pins_bloom_filter_tandem_non_acgt_and_boundaries() raises:
    """`unique_pins` keeps pins whose centered 19-mer is unique in both sequences and rejects tandem repeats,
    ambiguous bases, mismatched k-mers, and boundary pins."""
    seed(19)
    var unique_left = _random(120)
    var repeat_unit = _random(30)
    var unique_mid = _random(120)
    var unique_right = _random(120)
    var ref_seq = List[UInt8]()
    ref_seq.extend(unique_left.copy())
    ref_seq.extend(repeat_unit.copy())
    ref_seq.extend(repeat_unit.copy())
    ref_seq.extend(unique_mid.copy())
    ref_seq.extend(unique_right.copy())
    var query = ref_seq.copy()
    # Place an N in the middle region and a mismatch in the right region.
    var n_pos = 120 + 60 + 60
    ref_seq[n_pos] = UInt8(ord("N"))
    query[n_pos] = UInt8(ord("N"))
    var mm_pos = 120 + 60 + 120 + 60
    query[mm_pos + 2] = UInt8(ord("A")) if ref_seq[mm_pos + 2] != UInt8(ord("A")) else UInt8(ord("C"))

    var found = Pins()
    # 0: boundary (too close to start, < UNIQUE_K // 2)
    found.r.append(3)
    found.q.append(3)
    # 1: valid unique pin in unique_left
    found.r.append(60)
    found.q.append(60)
    # 2: pin inside tandem repeat (occurs twice in both)
    found.r.append(120 + 15)
    found.q.append(120 + 15)
    # 3: pin whose 19-mer contains N
    found.r.append(n_pos)
    found.q.append(n_pos)
    # 4: pin whose 19-mer differs between ref and query
    found.r.append(mm_pos)
    found.q.append(mm_pos)
    # 5: valid unique pin past mismatch
    found.r.append(mm_pos + 30)
    found.q.append(mm_pos + 30)
    # 6: boundary (too close to end)
    found.r.append(len(ref_seq) - 2)
    found.q.append(len(query) - 2)

    var keep = unique_pins(ref_seq, query, found)
    assert_equal(len(keep), 7)
    assert_false(keep[0])
    assert_true(keep[1])
    assert_false(keep[2])
    assert_false(keep[3])
    assert_false(keep[4])
    assert_true(keep[5])
    assert_false(keep[6])


def test_window_scratch_cache_hit_invalidation_and_shrink() raises:
    """Reusing `WindowScratch` across unchanged consensus paths (cache hit), mutated consensus paths (cache
    invalidation), and smaller graphs matches fresh `WindowScratch` calls bit-for-bit."""
    seed(41)
    var s0 = _random(1400)
    var q1 = _mutated(s0, 0.02)
    var q2 = _mutated(s0, 0.02)
    var g1 = Graph()
    _ = g1.add("s0", Span(s0))
    var t1 = Topology(g1)
    var scoring = Scoring()

    var reused = WindowScratch()
    var ws_reused = Workspace()
    var row_reused = List[Int](length=len(g1), fill=-1)

    # First alignment populates the cache.
    var a1 = align_pinned(g1, t1, q1, scoring, 10, 0.01, row_reused, reused, ws_reused, threads=4)
    assert_true(reused.has_cached_index)
    assert_equal(reused.cached_letters, s0)

    # Second alignment on same graph hits the cache and matches fresh scratch.
    var fresh_sc = WindowScratch()
    var fresh_ws = Workspace()
    var fresh_row = List[Int](length=len(g1), fill=-1)
    var a2 = align_pinned(g1, t1, q2, scoring, 10, 0.01, row_reused, reused, ws_reused, threads=4)
    var a2_fresh = align_pinned(g1, t1, q2, scoring, 10, 0.01, fresh_row, fresh_sc, fresh_ws, threads=1)
    assert_equal(a2.score, a2_fresh.score)
    assert_equal(a2.nodes, a2_fresh.nodes)

    # Third alignment on a smaller, distinct graph invalidates the cache and shrinks worker row_of safely.
    var small = _random(700)
    var q_small = _mutated(small, 0.02)
    var g2 = Graph()
    _ = g2.add("small", Span(small))
    var t2 = Topology(g2)
    row_reused.resize(len(g2), -1)
    var a3 = align_pinned(g2, t2, q_small, scoring, 10, 0.01, row_reused, reused, ws_reused, threads=4)
    assert_equal(reused.cached_letters, small)
    var fresh_sc2 = WindowScratch()
    var fresh_ws2 = Workspace()
    var fresh_row2 = List[Int](length=len(g2), fill=-1)
    var a3_fresh = align_pinned(g2, t2, q_small, scoring, 10, 0.01, fresh_row2, fresh_sc2, fresh_ws2, threads=1)
    assert_equal(a3.score, a3_fresh.score)
    assert_equal(a3.nodes, a3_fresh.nodes)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
