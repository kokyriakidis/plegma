from std.testing import assert_equal, assert_true, TestSuite

from karan import Costs, Mode, align
from sporos import read_fasta


def test_karan_align() raises:
    var aligned = align("ACGTACGTTTGCA", "ACGTCGTTTTGCA")
    assert_equal(aligned.cost, 2)
    assert_true(aligned.cigar == "4=1D2=1I6=")
    var affine = align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", Costs.affine(4, 6, 2))
    assert_equal(affine.cost, 12)


def test_sporos_reads_fasta_and_fastq() raises:
    var path = "build/test_deps.fq"
    with open(path, "w") as f:
        f.write("@a first\nACGT\n+\nIIII\n@b\nGGCCA\n+\nIIIII\n")
    var records = read_fasta(path)
    assert_equal(len(records), 2)
    assert_true(records.names[0] == "a")
    assert_equal(len(records.seq(1)), 5)
    assert_equal(records.seq(1)[4], UInt8(ord("A")))


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
