from std.random import random_ui64, seed
from std.testing import assert_equal, assert_true, TestSuite

from karan import distance

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


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
