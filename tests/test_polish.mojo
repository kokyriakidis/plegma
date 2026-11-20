from std.random import random_ui64, seed
from std.testing import assert_equal, assert_raises, assert_true, TestSuite

from karan import distance

from plegma.consensus import ConsensusOptions, consensus
from plegma.kernel import BOTH
from plegma.polish import polish


def _random(length: Int) -> List[UInt8]:
    var out = List[UInt8](capacity=length)
    var letters = "ACGT".as_bytes()
    for _ in range(length):
        out.append(letters[Int(random_ui64(0, 3))])
    return out^


def _mutated(seq: List[UInt8], rate: Float64) -> List[UInt8]:
    """`seq` with substitutions, deletions and insertions, each at `rate` / 3 per letter."""
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


def _text(seq: List[UInt8]) -> String:
    return String(unsafe_from_utf8=seq.copy())


def test_polish_never_raises_the_summed_distance() raises:
    seed(5)
    for trial in range(40):
        var truth = _random(Int(random_ui64(1, 400)))
        var reads = List[List[UInt8]]()
        var members = List[Int]()
        for k in range(Int(random_ui64(1, 12))):
            reads.append(_mutated(truth, 0.15))
            members.append(k)
        var start = _mutated(truth, 0.1)
        var polished = polish(start, reads, members, 30)
        var before = 0
        var after = 0
        for r in reads:
            before += distance(_text(start), _text(r))
            after += distance(_text(polished), _text(r))
        assert_true(after <= before, String(trial))


def test_polish_mends_one_error_in_one_round() raises:
    seed(3)
    var truth = _random(300)
    var reads = List[List[UInt8]]()
    var members = List[Int]()
    for k in range(9):
        reads.append(_mutated(truth, 0.03))
        members.append(k)
    var start = List[UInt8](truth[0:120])
    start.extend(List[UInt8](truth[121:]))
    var polished = polish(start, reads, members, 1)
    assert_equal(distance(_text(polished), _text(truth)), 0)


def test_polish_mends_a_tandem_repeat_consensus() raises:
    # A tandem array; the consensus starts one unit short and with a letter missing from a run.
    seed(11)
    var truth = _random(150)
    for _ in range(20):
        for c in "TGGTGA".as_bytes():
            truth.append(c)
    truth.extend(_random(150))
    var start = List[UInt8](truth[0:200])
    start.extend(List[UInt8](truth[206:240]))
    start.extend(List[UInt8](truth[241:]))
    var reads = List[List[UInt8]]()
    var members = List[Int]()
    for k in range(15):
        reads.append(_mutated(truth, 0.05))
        members.append(k)
    var polished = polish(start, reads, members, 20)
    assert_equal(distance(_text(polished), _text(truth)), 0)
    # Never worse in the reads' summed distance than where it started.
    var before = 0
    var after = 0
    for r in reads:
        before += distance(_text(start), _text(r))
        after += distance(_text(polished), _text(r))
    assert_true(after < before)


def test_likelihood_alone_mends_a_letter_and_a_unit() raises:
    seed(13)
    var truth = _random(150)
    for _ in range(12):
        for c in "TGGTGA".as_bytes():
            truth.append(c)
    truth.extend(_random(150))
    var reads = List[List[UInt8]]()
    var members = List[Int]()
    for k in range(12):
        reads.append(_mutated(truth, 0.04))
        members.append(k)
    # One letter missing from the flank, one unit from the array.
    var start = List[UInt8](truth[0:60])
    start.extend(List[UInt8](truth[61:180]))
    start.extend(List[UInt8](truth[186:]))
    var polished = polish(start, reads, members, 0, likely_rounds=20)
    assert_equal(distance(_text(polished), _text(truth)), 0)


def test_an_empty_read_leaves_the_consensus_alone() raises:
    """An empty read between the anchors among real ones: polishing leaves it out (all its alignment is edits), so
    the consensus is the one the others make; its deletions once set the likelihood model's rates (334 letters, not
    303)."""
    seed(23)
    var truth = _random(300)
    var reads = List[List[UInt8]]()
    for _ in range(3):
        reads.append(_mutated(truth, 0.06))
    var options = ConsensusOptions()
    options.poa.seeded = False
    var without = consensus(reads.copy(), options=options)
    reads.append(List[UInt8]())
    var with_empty = consensus(reads.copy(), options=options)
    assert_equal(_text(with_empty[0]), _text(without[0]))


def test_consensus_checks_its_arguments() raises:
    var reads = List[List[UInt8]]()
    reads.append(_random(50))
    reads.append(_random(50))
    with assert_raises(contains="anchors for"):
        _ = consensus(reads.copy(), [BOTH])
    with assert_raises(contains="anchor 7"):
        _ = consensus(reads.copy(), [BOTH, 7])
    with assert_raises(contains="strand 9"):
        _ = consensus(reads.copy(), strands=[0, 9])
    assert_equal(len(consensus(List[List[UInt8]]())), 1)


def test_unit_deletion_out_of_a_reads_band_does_not_wrap() raises:
    """Reads whose band cannot reach both ends of a repeat unit gain nothing by deleting it: ten of them (an
    insertion of 60 letters, past both bands, inside the repeat) wrapped the unit deletion's gain to the top of 32
    bits, the round was reverted, and the flank's error was never mended."""
    seed(1)
    var f1 = _random(150)
    var f2 = _random(150)
    var x = _random(60)
    var unit = "TGGTGA".as_bytes()
    var short = f1.copy()
    for _ in range(20):
        short.extend(List[UInt8](unit))
    short.extend(f2.copy())
    var long = f1.copy()
    for _ in range(10):
        long.extend(List[UInt8](unit))
    long.extend(x.copy())
    for _ in range(10):
        long.extend(List[UInt8](unit))
    long.extend(f2.copy())
    var consensus = short.copy()
    var at = len(short) - 75
    consensus[at] = UInt8(ord("A")) if consensus[at] != UInt8(ord("A")) else UInt8(ord("C"))
    var reads = List[List[UInt8]]()
    var members = List[Int]()
    for k in range(19):
        reads.append(long.copy() if k < 10 else short.copy())
        members.append(k)
    var polished = polish(consensus, reads, members, 10, likely_rounds=0)
    for t in range(1, 151):
        assert_equal(polished[len(polished) - t], short[len(short) - t])


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
