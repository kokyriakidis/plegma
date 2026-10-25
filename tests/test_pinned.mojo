from std.random import random_ui64, seed
from std.testing import assert_equal, assert_true, TestSuite

from karan import distance

from plegma.graph import Graph
from plegma.kernel import Scoring, Topology, Workspace
from plegma.pinned import WindowScratch, align_pinned


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


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
