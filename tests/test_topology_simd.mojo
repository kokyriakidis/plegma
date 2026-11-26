"""Unit tests for `plegma/topology.mojo`, `plegma/kernel_simd.mojo`, and `Graph.add_twin_path`/`node_order_into`."""

from std.sys import simd_width_of
from std.testing import TestSuite, assert_equal, assert_true

from plegma.graph import Graph
from plegma.kernel_simd import _first_equal, _last, _last_equal, _moved, _neg, _prefix_max, _rebased, _up
from plegma.topology import Topology


def _bytes(s: String) -> List[UInt8]:
    return List[UInt8](s.as_bytes())


def test_add_twin_path_matches_fuse_nodes() raises:
    var g1 = Graph()
    var g2 = Graph()
    var s0 = _bytes("ACGTACGT")
    var s1 = _bytes("ACGTTCGT")
    var p0_1 = g1.add("r0", Span(s0))
    var p0_2 = g2.add("r0", Span(s0))
    var aligned: List[Int] = [0, 1, 2, 3, 4, 5, 6, 7]
    _ = g1.fuse_nodes("r1", Span(s1), aligned)
    _ = g2.fuse_nodes("r1", Span(s1), aligned)

    # Add identical copy of r0 via fuse_nodes vs add_twin_path.
    var nodes = List[Int]()
    for node in g1.paths[p0_1]:
        nodes.append(Int(node))
    _ = g1.fuse_nodes("r2", Span(s0), nodes)
    _ = g2.add_twin_path("r2", p0_2)

    assert_equal(len(g1), len(g2))
    assert_equal(len(g1.edge_weight), len(g2.edge_weight))
    for e in range(len(g1.edge_weight)):
        assert_equal(g1.edge_weight[e], g2.edge_weight[e])
    assert_equal(g1.msa(), g2.msa())


def test_reusable_topology_matches_fresh_across_graphs() raises:
    var topo = Topology()
    var incoming = List[Int32]()
    var ready = List[Int]()
    var full = List[Int]()

    # First graph: larger with branch.
    var g1 = Graph()
    _ = g1.add("a", Span(_bytes("ACGTACGTACGT")))
    var aligned1: List[Int] = [0, 1, 2, -1, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    _ = g1.fuse_nodes("b", Span(_bytes("ACGTTACGTACGT")), aligned1)
    topo.compute(g1, incoming, ready)
    var fresh1 = Topology(g1)
    assert_equal(topo.order, fresh1.order)
    assert_equal(topo.rank, fresh1.rank)
    assert_equal(topo.remaining, fresh1.remaining)

    # Reuse on a smaller graph with an allowed mask.
    var g2 = Graph()
    _ = g2.add("c", Span(_bytes("GATTACA")))
    var allowed: List[Bool] = [True, True, False, True, True, True, False]
    topo.compute_allowed(g2, allowed, incoming, ready, full)
    var fresh2 = Topology(g2, allowed)
    assert_equal(topo.order, fresh2.order)
    assert_equal(topo.rank, fresh2.rank)
    assert_equal(topo.remaining, fresh2.remaining)


def test_simd_prefix_and_shift_helpers() raises:
    comptime for dt in [DType.int16, DType.int32]:
        comptime W = simd_width_of[dt]()
        var v = SIMD[dt, W](0)
        for i in range(W):
            v[i] = Scalar[dt]((i * 7 + 3) % 13 - 4)
        var pref = _prefix_max(v)
        var running = v[0]
        for i in range(W):
            running = max(running, v[i])
            assert_equal(pref[i], running)
        var last_v = _last(v)
        for i in range(W):
            assert_equal(last_v[i], v[W - 1])
        var fill = SIMD[dt, W](_neg[dt]())
        var shifted = _up[1](v, fill)
        assert_equal(shifted[0], _neg[dt]())
        for i in range(1, W):
            assert_equal(shifted[i], v[i - 1])


def test_simd_first_last_equal_and_rebase_straddling_lengths() raises:
    for n in [1, 7, 15, 16, 17, 31, 32, 33, 65]:
        var row = List[Int16](length=n, fill=Int16(10))
        row[0] = Int16(42)
        row[n - 1] = Int16(42)
        assert_equal(_first_equal[DType.int16](row.unsafe_ptr(), n, Int16(42)), 0)
        assert_equal(_last_equal[DType.int16](row.unsafe_ptr(), n, Int16(42)), n - 1)
        if n >= 5:
            row[0] = Int16(5)
            row[n - 1] = Int16(5)
            row[2] = Int16(99)
            row[n - 2] = Int16(99)
            assert_equal(_first_equal[DType.int16](row.unsafe_ptr(), n, Int16(99)), 2)
            assert_equal(_last_equal[DType.int16](row.unsafe_ptr(), n, Int16(99)), n - 2)
        var moved = List[Int16](length=n, fill=0)
        row[0] = _neg[DType.int16]()
        _moved[DType.int16](row.unsafe_ptr(), moved.unsafe_ptr(), n, -5)
        assert_equal(moved[0], _neg[DType.int16]())
        if n > 1:
            assert_equal(moved[1], row[1] - Int16(5))
        _rebased[DType.int16](moved.unsafe_ptr(), n, Int16(3))
        assert_equal(moved[0], _neg[DType.int16]())
        assert_true(len(moved) == n)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
