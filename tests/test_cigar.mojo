"""Unit tests for `plegma/cigar.mojo`: CIGAR parsing, pairing extraction, and coordinate walkers."""

from std.testing import TestSuite, assert_equal, assert_raises

from karan import Alignment

from plegma.cigar import (
    DELETION,
    INSERTION,
    MATCH,
    MISMATCH,
    append_pairing_from_cigar,
    cigar_runs,
    haplotype_at,
    pairing_from_cigar,
    query_at,
)


def test_cigar_runs_parses_and_rejects_missing_length() raises:
    var runs = cigar_runs("3=1X2I1D4=")
    assert_equal(len(runs.ops), 5)
    assert_equal(runs.ops[0], MATCH)
    assert_equal(runs.counts[0], 3)
    assert_equal(runs.ops[1], MISMATCH)
    assert_equal(runs.counts[1], 1)
    assert_equal(runs.ops[2], INSERTION)
    assert_equal(runs.counts[2], 2)
    assert_equal(runs.ops[3], DELETION)
    assert_equal(runs.counts[3], 1)
    assert_equal(runs.ops[4], MATCH)
    assert_equal(runs.counts[4], 4)
    with assert_raises():
        _ = cigar_runs("=3X")


def test_pairing_from_cigar_with_offset() raises:
    var p = pairing_from_cigar("2=1I1D1X1=", 5, ref_offset=10)
    assert_equal(len(p), 5)
    assert_equal(p[0], 10)
    assert_equal(p[1], 11)
    assert_equal(p[2], -1)
    assert_equal(p[3], 13)
    assert_equal(p[4], 14)

    var buf = List[Int]()
    append_pairing_from_cigar("1I2=", 4, buf)
    assert_equal(len(buf), 3)
    assert_equal(buf[0], -1)
    assert_equal(buf[1], 4)
    assert_equal(buf[2], 5)


def test_query_at_and_haplotype_at() raises:
    # Ref span [2, 7) (5 bases), query span [1, 6) (5 bases): 2= 1I 1D 2=
    var aln = Alignment(2, 0, "2=1I1D2=", 2, 7, 1, 6)
    var qa = query_at(aln)
    assert_equal(len(qa), 6)
    assert_equal(qa[0], 1)
    assert_equal(qa[1], 2)
    assert_equal(qa[2], 4)  # insertion advances query before ref[2]
    assert_equal(qa[3], 4)  # deletion keeps query at 4
    assert_equal(qa[4], 5)
    assert_equal(qa[5], 6)

    # Haplotype inside read: query length 5, ref starts at 10: 2= 1D 1I 2=
    var hap = Alignment(2, 0, "2=1D1I2=", 10, 15, 0, 5)
    var ha = haplotype_at(hap)
    assert_equal(len(ha), 6)
    assert_equal(ha[0], 10)
    assert_equal(ha[1], 11)
    assert_equal(ha[2], 13)  # deletion in read advances read before next hap base
    assert_equal(ha[3], 13)  # insertion in hap keeps read at 13
    assert_equal(ha[4], 14)
    assert_equal(ha[5], 15)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
