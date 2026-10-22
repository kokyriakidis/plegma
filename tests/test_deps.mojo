from std.testing import assert_equal, assert_true, TestSuite

from karan import Costs, Mode, align
from sporos import Reference


def test_karan_align() raises:
    var aligned = align("ACGTACGTTTGCA", "ACGTCGTTTTGCA")
    assert_equal(aligned.cost, 2)
    assert_true(aligned.cigar == "4=1D2=1I6=")
    var affine = align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", Costs.affine(4, 6, 2))
    assert_equal(affine.cost, 12)


def test_sporos_mums() raises:
    var ref_bytes = List[UInt8]("GATTACAGGCCTTAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC".as_bytes())
    var query = List[UInt8]("TTTTGCATGCCTGCAGGTCGACTCTAGAGGAAAA".as_bytes())
    var seqs = List[Span[UInt8, origin_of(ref_bytes)]]()
    seqs.append(Span(ref_bytes))
    var reference = Reference(seqs)
    var mums = reference.mems(Span(query), min_len=20, mode="mum")
    assert_equal(len(mums), 1)
    assert_equal(mums[0].ref_pos, 17)
    assert_equal(mums[0].query_pos, 2)
    assert_equal(mums[0].length, 29)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
