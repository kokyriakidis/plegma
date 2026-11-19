from std.math import log10
from std.random import random_ui64, seed
from std.testing import assert_equal, assert_true, TestSuite

from plegma.homopolymer import FLANK_UNKNOWN, STRAND_UNKNOWN, HomopolymerModel, reads_say_one
from plegma.polish import polish

comptime PATH = "build/test_homopolymer.csv"


def _rows(mut text: String, base: String, strand: Int, left: String, right: String, n: Int, peak: Int, count: Int):
    """One (base, strand, flanks, n) row of the model, m 0 to 4: 0.9 at `peak`, 0.025 elsewhere."""
    for m in range(5):
        var p = 0.9 if m == peak else 0.025
        text.write(base, ",", strand, ",", left, ",", right, ",", n, ",", m, ",")
        text.write(count if m == peak else 0, ",", p, ",", 10.0 * log10(p), "\n")


def _model() raises -> HomopolymerModel:
    """shasta2's own test model: true lengths 2 and 3, observed 0 to 4, every base alike; in the * rows strand 0 reads
    a run right 90% of the time, strand 1 one letter short 90%; for A between C and G strand 0 reads it short too."""
    var text = String("base,strand,left,right,n,m,count,probability,logProbabilityDb\n")
    for base in [String("A"), String("C"), String("G"), String("T")]:
        for strand in range(2):
            for n in range(2, 4):
                _rows(text, base, strand, "*", "*", n, n if strand == 0 else n - 1, 100)
    for n in range(2, 4):
        _rows(text, "A", 0, "C", "G", n, n - 1, 0)
    with open(PATH, "w") as f:
        f.write(text)
    return HomopolymerModel.load(PATH)


def test_strands_and_flanks_as_shasta2_tests_them() raises:
    var model = _model()
    var any = FLANK_UNKNOWN
    # A strand 0 read of 3 and a strand 1 read of 2 both say 3; as one strand they tie, and a tie goes to 2.
    assert_equal(model.most_likely(0, any, any, [3, 2], [0, 1]), 3)
    assert_equal(model.most_likely(0, any, any, [3, 2], [0, 0]), 2)
    assert_equal(model.most_likely(0, any, any, [2, 2], [1, 1]), 3)
    # Between C and G a strand 0 read of 2 says 3; between other flanks, with no row of their own, 2.
    assert_equal(model.most_likely(0, 1, 2, [2, 2], [0, 0]), 3)
    assert_equal(model.most_likely(0, 2, 1, [2, 2], [0, 0]), 2)
    assert_equal(model.most_likely(0, 3, 2, [2, 2], [0, 0]), 2)
    # Strand unknown: the two strands half each; strand 0's read of 3 then outweighs nothing else.
    assert_equal(model.most_likely(0, any, any, [3], [STRAND_UNKNOWN]), 3)


def test_polish_takes_the_runs_most_likely_length() raises:
    # A truth with one run (TTT) and no other: three strand 0 reads report it whole, five strand 1 reads one short,
    # as the model says strand 1 reads it. Their vote makes it TT; the model, TTT.
    seed(7)
    var model = _model()
    # Random A, C and G, no letter twice in a row, and TTT in the middle.
    var letters = "ACG".as_bytes()
    var truth = List[UInt8]()
    var short = List[UInt8]()
    for i in range(120):
        if i == 60:
            for t in range(3):
                truth.append(UInt8(ord("T")))
                if t > 0:
                    short.append(UInt8(ord("T")))
            continue
        var c = letters[Int(random_ui64(0, 2))]
        while len(truth) > 0 and c == truth[len(truth) - 1]:
            c = letters[Int(random_ui64(0, 2))]
        truth.append(c)
        short.append(c)
    var reads = List[List[UInt8]]()
    var strands = List[Int]()
    var members = List[Int]()
    for k in range(8):
        reads.append(truth.copy() if k < 3 else short.copy())
        strands.append(0 if k < 3 else 1)
        members.append(k)
    var plain = polish(short, reads, members, 50, likely_rounds=0)
    assert_true(plain == short)
    var modelled = polish(short, reads, members, 50, likely_rounds=0, homopolymers=model, strands=strands)
    assert_true(modelled == truth)


def test_reads_say_one_counts_only_reads_with_a_length() raises:
    # Reads with no length there (-1) do not vote: counted, they took the lower median to 1.
    assert_true(not reads_say_one(_lengths("-1,-1,1,1,3,3,3")))
    assert_true(reads_say_one(_lengths("1,1,1,-1,-1,-1,-1")))
    assert_true(not reads_say_one(_lengths("1,1,-1,-1")))


def _lengths(text: String) raises -> List[Int]:
    var out = List[Int]()
    for word in text.split(","):
        out.append(Int(word))
    return out^


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
