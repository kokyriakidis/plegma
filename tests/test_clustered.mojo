from std.math import log10
from std.random import random_ui64, seed
from std.testing import assert_equal, assert_raises, assert_true, TestSuite

from karan import distance

from karan import Aligner, Costs, Mode, Ties
from plegma.clustered import clustered_consensus
from plegma.phased import _haplotype_at, phased_consensus
from plegma.runs import _length_in, _run_pairs
from plegma.consensus import ConsensusOptions, consensus
from plegma.homopolymer import HomopolymerModel
from plegma.kernel import BOTH, INSIDE, LEFT, RIGHT


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


def _b(text: StringSlice) -> List[UInt8]:
    return List[UInt8](text.as_bytes())


def _text(seq: List[UInt8]) -> String:
    return String(unsafe_from_utf8=seq)


def _other_haplotype(seq: List[UInt8]) -> List[UInt8]:
    """`seq` with a SNP every 60 letters and a 3-letter deletion and insertion: another allele of the place."""
    var out = List[UInt8]()
    var letters = "ACGT".as_bytes()
    for i in range(len(seq)):
        if i == 200 or i == 201 or i == 202:
            continue
        if i % 60 == 30:
            var code = 0
            for b in range(4):
                if letters[b] == seq[i]:
                    code = b
            out.append(letters[(code + 1) % 4])
            continue
        out.append(seq[i])
        if i == 400:
            out.append(letters[1])
            out.append(letters[2])
            out.append(letters[3])
    return out^


def _unpinned() -> ConsensusOptions:
    var options = ConsensusOptions()
    options.poa.seeded = False
    return options^


def _spells(cigar: String, consensus: List[UInt8], letters: List[UInt8], c0: Int, r0: Int) -> Bool:
    """Whether `cigar` from consensus `c0` and read `r0` pairs `=` letters equal and `X` letters different."""
    var i = c0
    var j = r0
    var count = 0
    for c in cigar.as_bytes():
        if c >= UInt8(ord("0")) and c <= UInt8(ord("9")):
            count = count * 10 + Int(c - UInt8(ord("0")))
            continue
        for _ in range(count):
            if c == UInt8(ord("=")) or c == UInt8(ord("X")):
                if (consensus[i] == letters[j]) != (c == UInt8(ord("="))):
                    return False
                i += 1
                j += 1
            elif c == UInt8(ord("I")):
                j += 1
            else:
                i += 1
        count = 0
    return True


def test_two_haplotypes_split_and_each_consensus_is_its_truth() raises:
    seed(41)
    var one = _random(800)
    var two = _other_haplotype(one)
    var reads = List[List[UInt8]]()
    var truth_of = List[Int]()
    for k in range(20):
        if k % 2 == 0:
            reads.append(_mutated(one, 0.04))
        else:
            reads.append(_mutated(two, 0.04))
        truth_of.append(k % 2)
    var got = clustered_consensus(reads.copy(), max_clusters=2, min_freq=0.2, options=_unpinned())
    assert_equal(len(got.consensus), 2)
    for k in range(2):
        var a = distance(_text(got.consensus[k]), _text(one))
        var b = distance(_text(got.consensus[k]), _text(two))
        assert_true(min(a, b) <= 1, String(t"cluster {k}: {a} and {b} from the two haplotypes"))
    # Every read in its haplotype's cluster, and its CIGAR spelling it against that consensus.
    for r in range(len(reads)):
        ref placed = got.placements[r]
        ref mine = got.consensus[placed.cluster]
        ref other = got.consensus[1 - placed.cluster]
        var truth = one.copy() if truth_of[r] == 0 else two.copy()
        assert_true(distance(_text(mine), _text(truth)) < distance(_text(other), _text(truth)), String(t"read {r}"))
        assert_true(
            _spells(placed.cigar, mine, reads[r], placed.consensus_start, placed.read_start), String(t"read {r}")
        )
        assert_equal(placed.read_end - placed.read_start, len(reads[r]))
    assert_equal(len(got.members[0]) + len(got.members[1]), len(reads))


def test_one_haplotype_is_one_cluster() raises:
    seed(43)
    var one = _random(600)
    var reads = List[List[UInt8]]()
    for _ in range(12):
        reads.append(_mutated(one, 0.05))
    var got = clustered_consensus(reads.copy(), max_clusters=2, min_freq=0.2, options=_unpinned())
    assert_equal(len(got.consensus), 1)
    assert_true(distance(_text(got.consensus[0]), _text(one)) <= 1)
    assert_equal(len(got.members[0]), 12)


def test_a_cluster_under_min_freq_joins_the_other() raises:
    seed(47)
    var one = _random(600)
    var two = _other_haplotype(one)
    var reads = List[List[UInt8]]()
    for k in range(10):
        if k == 0:
            reads.append(_mutated(two, 0.03))
        else:
            reads.append(_mutated(one, 0.03))
    var got = clustered_consensus(reads.copy(), max_clusters=2, min_freq=0.25, options=_unpinned())
    assert_equal(len(got.consensus), 1)
    assert_true(distance(_text(got.consensus[0]), _text(one)) <= 1)


def test_a_refinement_down_to_one_cluster() raises:
    """A refinement round that drops a cluster under the floor leaves one: one consensus of every read comes out
    (pgphase's chr20:22,644,019-22,644,096 region, 57 reads of a GAG repeat; it read past the clusters before)."""
    var reads = List[List[UInt8]]()
    var anchors = List[Int]()
    var strands = List[Int]()
    with open("tests/data/collapsed_clusters/reads.txt", "r") as f:
        for line in f.read().split("\n"):
            if line.byte_length() == 0:
                continue
            var parts = line.split(":")
            var seq = List[UInt8]()
            for b in parts[0].as_bytes():
                seq.append(b)
            reads.append(seq^)
            anchors.append(Int(parts[1]))
            strands.append(Int(parts[2]))
    var n = len(reads)
    var options = _unpinned()
    options.poa.scoring.match_score = 2
    options.poa.scoring.mismatch = 6
    options.poa.scoring.open1 = 6
    options.poa.scoring.extend1 = 2
    options.poa.scoring.open2 = 24
    options.poa.scoring.extend2 = 1
    options.two_lengths_db = 8.0
    var got = clustered_consensus(reads^, anchors, strands, 2, Float64(Float32(0.12)), options)
    assert_equal(len(got.consensus), 1)
    assert_equal(len(got.placements), n)
    for k in range(n):
        assert_equal(got.placements[k].cluster, 0)


def test_one_read_third_letter_does_not_hide_a_het() raises:
    # HG002 chr18:25,307,645, an A deleted from CATAGGGG: one MSA column tells the haplotypes apart, and one read
    # (another A after the G run) has a G there. Seeded at the furthest row alone, that read was the second mode and
    # the het one cluster.
    var with_a = List[UInt8]("AAATCGGCATAGGGGAAAGGGCCTG".as_bytes())
    var without = List[UInt8]("AAATCGGCATGGGGAAAGGGCCTG".as_bytes())
    var reads = List[List[UInt8]]()
    reads.append(List[UInt8]("AAATCGGCATAGGGGGAAAGGGCCTG".as_bytes()))
    reads.append(List[UInt8]("AAATCGGCATGGGGAAAAGGGCCTG".as_bytes()))
    for _ in range(8):
        reads.append(with_a.copy())
    for _ in range(8):
        reads.append(without.copy())
    var got = clustered_consensus(reads^, max_clusters=2, min_freq=0.2, options=_unpinned())
    assert_equal(len(got.consensus), 2)
    assert_true(
        (got.consensus[0] == with_a and got.consensus[1] == without)
        or (got.consensus[0] == without and got.consensus[1] == with_a)
    )


def _snp_and_run(letter: String, run: Int) -> List[UInt8]:
    """A haplotype `letter` apart from the other at one place, with a run of T `run` long further on."""
    var s = String("GATCCAGTACGGATCA") + letter + String("TGCAAC")
    for _ in range(run):
        s += "T"
    s += "GACCTAGGCTAA"
    return List[UInt8](s.as_bytes())


def test_a_letter_apart_outweighs_a_run_slip() raises:
    # HG002 chr20:11,477,892: the parents a letter apart, and each parent's reads split alike by a run's length (the
    # reads slip). Every MSA column counted the same, the reads split by the run, and neither consensus kept the
    # letter; a column of two letters counts SUBSTITUTION_WEIGHT gap columns.
    var reads = List[List[UInt8]]()
    for letter in ["C", "G"]:
        for k in range(20):
            reads.append(_snp_and_run(letter, 9 if k % 5 < 3 else 10))
    var got = clustered_consensus(reads^, max_clusters=2, min_freq=0.2, options=_unpinned())
    assert_equal(len(got.consensus), 2)
    var c = _snp_and_run("C", 9)
    var g = _snp_and_run("G", 9)
    assert_true(
        (got.consensus[0] == c and got.consensus[1] == g) or (got.consensus[0] == g and got.consensus[1] == c),
        String(t"{_text(got.consensus[0])} {_text(got.consensus[1])}"),
    )
    for k in range(2):
        assert_equal(len(got.members[k]), 20)


def test_anchored_and_inside_reads() raises:
    """Reads fixed at one anchor build the consensus with the spanning ones; a read inside builds none and is aligned
    within its haplotype's consensus."""
    seed(53)
    var one = _random(900)
    var two = _other_haplotype(one)
    var reads = List[List[UInt8]]()
    var anchors = List[Int]()
    for k in range(16):
        if k % 2 == 0:
            reads.append(_mutated(one, 0.04))
        else:
            reads.append(_mutated(two, 0.04))
        anchors.append(BOTH)
    reads.append(List[UInt8](_mutated(one, 0.04)[0:500]))
    anchors.append(LEFT)
    var tail = _mutated(two, 0.04)
    reads.append(List[UInt8](tail[len(tail) - 500 : len(tail)]))
    anchors.append(RIGHT)
    reads.append(List[UInt8](one[300:600]))
    anchors.append(INSIDE)
    var got = clustered_consensus(reads.copy(), anchors, max_clusters=2, min_freq=0.2, options=_unpinned())
    assert_equal(len(got.consensus), 2)
    ref inside = got.placements[len(reads) - 1]
    assert_true(distance(_text(got.consensus[inside.cluster]), _text(one)) <= 1)
    assert_equal(inside.cost, 0)
    assert_equal(inside.consensus_end - inside.consensus_start, 300)


def _wide_model() raises -> HomopolymerModel:
    """Every base, strand and true length n 2 to 40 read m long: n with 0.8, one off 0.07 each side, two off 0.02,
    any other m 1e-4."""
    var text = String("base,strand,left,right,n,m,count,probability,logProbabilityDb\n")
    for base in [String("A"), String("C"), String("G"), String("T")]:
        for strand in range(2):
            for n in range(2, 41):
                for m in range(46):
                    var off = abs(m - n)
                    var p = 0.8 if off == 0 else (0.07 if off == 1 else (0.02 if off == 2 else 1e-4))
                    text.write(base, ",", strand, ",*,*,", n, ",", m, ",", 100, ",", p, ",", 10.0 * log10(p), "\n")
    comptime path = "build/test_clustered_model.csv"
    with open(path, "w") as f:
        f.write(text)
    return HomopolymerModel.load(path)


def _with_run(left: List[UInt8], length: Int, right: List[UInt8]) -> List[UInt8]:
    var out = left.copy()
    for _ in range(length):
        out.append(UInt8(ord("A")))
    out.extend(right.copy())
    return out^


def _no_a(length: Int) -> List[UInt8]:
    var out = List[UInt8](capacity=length)
    var letters = "CGT".as_bytes()
    for _ in range(length):
        out.append(letters[Int(random_ui64(0, 2))])
    return out^


def test_homopolymer_lengths_far_apart_stay_two() raises:
    """Haplotypes apart only by a run of A 20 and 27 long, the reads each a letter off now and then: with the model
    both the clustered and the phased consensus keep the two lengths (the search spans both, not 4 about one)."""
    seed(59)
    var left = _no_a(150)
    var right = _no_a(150)
    var reads = List[List[UInt8]]()
    var strands = List[Int]()
    var haplotypes = List[Int]()
    var offsets = [0, 0, -1, 0, 1, 0, 0, -1]
    for h in range(2):
        for k in range(len(offsets)):
            reads.append(_with_run(left, (20 if h == 0 else 27) + offsets[k], right))
            strands.append(k % 2)
            haplotypes.append(h)
    var model = _wide_model()
    var options = _unpinned()
    var one = _with_run(left, 20, right)
    var two = _with_run(left, 27, right)
    var clustered = clustered_consensus(
        reads.copy(), strands=strands, max_clusters=2, min_freq=0.2, options=options, homopolymers=model
    )
    assert_equal(len(clustered.consensus), 2)
    var phased = phased_consensus(reads.copy(), haplotypes, strands=strands, options=options, homopolymers=model)
    assert_true(phased.consensus[0] == one, _text(phased.consensus[0]))
    assert_true(phased.consensus[1] == two, _text(phased.consensus[1]))
    for got in [clustered.consensus.copy(), phased.consensus.copy()]:
        assert_true((got[0] == one and got[1] == two) or (got[0] == two and got[1] == one))


def test_runs_longer_than_the_model_keep_their_length() raises:
    """Haplotypes apart by a run of A 20 and 45 long, the model's rows ending at 40: the model cannot judge the long
    run, so each haplotype keeps its own consensus's length there (it was cut to the model's longest row, 40)."""
    seed(61)
    var left = _no_a(150)
    var right = _no_a(150)
    var reads = List[List[UInt8]]()
    var strands = List[Int]()
    var haplotypes = List[Int]()
    var offsets = [0, 0, -1, 0, 1, 0, 0, -1]
    for h in range(2):
        for k in range(len(offsets)):
            reads.append(_with_run(left, (20 if h == 0 else 45) + offsets[k], right))
            strands.append(k % 2)
            haplotypes.append(h)
    var model = _wide_model()
    var one = _with_run(left, 20, right)
    var two = _with_run(left, 45, right)
    var phased = phased_consensus(reads.copy(), haplotypes, strands=strands, options=_unpinned(), homopolymers=model)
    assert_true(phased.consensus[0] == one, _text(phased.consensus[0]))
    assert_true(phased.consensus[1] == two, _text(phased.consensus[1]))
    var clustered = clustered_consensus(
        reads.copy(), strands=strands, max_clusters=2, min_freq=0.2, options=_unpinned(), homopolymers=model
    )
    assert_equal(len(clustered.consensus), 2)
    var got = clustered.consensus.copy()
    assert_true((got[0] == one and got[1] == two) or (got[0] == two and got[1] == one))


def _short_reads_model() raises -> HomopolymerModel:
    """Every base, strand and true length n 2 to 40 read m long, a run read short likelier than long (as HiFi reads
    one): n with 0.8, one short 0.1, one long 0.05, more short 0.01, more long 1e-6."""
    var text = String("base,strand,left,right,n,m,count,probability,logProbabilityDb\n")
    for base in [String("A"), String("C"), String("G"), String("T")]:
        for strand in range(2):
            for n in range(2, 41):
                for m in range(46):
                    var d = m - n
                    var p = 0.8 if d == 0 else (0.1 if d == -1 else (0.05 if d == 1 else (0.01 if d < 0 else 1e-6)))
                    text.write(base, ",", strand, ",*,*,", n, ",", m, ",", 100, ",", p, ",", 10.0 * log10(p), "\n")
    comptime path = "build/test_clustered_short_model.csv"
    with open(path, "w") as f:
        f.write(text)
    return HomopolymerModel.load(path)


def test_run_repair_of_two_haplotypes_lengthens_no_run_twice() raises:
    """HG002 chr18:46,302,804, both haplotypes' reads in one consensus: one has T4 C T11 where the other has T16.
    The other's reads, without that C, reported their T16 for both runs and the model made both 16, a consensus of
    86 letters from reads of 69; a read with no letter aligned to a run's flank now does not vote there."""
    var t16 = List[UInt8]("GTTCCACATTTCTCTGTGGCTCTCTTCATTTTCTTTTCTTTTCTTTTTTTTTTTTTTTTATGAGATGGA".as_bytes())
    var t4c = List[UInt8]("GTTCCACATTTCTCTGTGGCTCTCTTCATTTTCTTTTCTTTTCTTTTCTTTTTTTTTTTATGAGATGGA".as_bytes())
    var reads = List[List[UInt8]]()
    var strands = List[Int]()
    for k in range(55):
        reads.append(t4c.copy() if k % 2 == 0 else t16.copy())
        strands.append(k % 3 % 2)
    var options = ConsensusOptions()
    options.poa.seeded = False
    var parts = consensus(reads^, List[Int](length=55, fill=BOTH), strands, options, _short_reads_model())
    assert_equal(len(parts), 1)
    assert_true(parts[0] == t16 or parts[0] == t4c, _text(parts[0]))


def test_one_letter_where_the_other_haplotype_has_a_run() raises:
    """HG002 chr20:8,254,476: one haplotype's 42 reads have T where the other's 26 have TTT. The model holds no rows
    under 2 and the haplotypes' run search started at 2, so the first came out TT (the truth's 2 bp deletion written
    as 1 bp); its reads, all T, now keep it one letter."""
    var one = List[UInt8]("ACATGTCTCTGAGTCACTGTTACCATTTATAAAATGAAA".as_bytes())
    var three = List[UInt8]("ACATGTCTCTGAGTCACTTTGTTACCATTTATAAAATGAAA".as_bytes())
    var reads = List[List[UInt8]]()
    var strands = List[Int]()
    var haplotypes = List[Int]()
    for k in range(68):
        reads.append(one.copy() if k < 42 else three.copy())
        strands.append(k % 2)
        haplotypes.append(0 if k < 42 else 1)
    var options = ConsensusOptions()
    options.poa.seeded = False
    var phased = phased_consensus(
        reads.copy(), haplotypes, strands=strands, options=options, homopolymers=_short_reads_model()
    )
    assert_true(phased.consensus[0] == one, _text(phased.consensus[0]))
    assert_true(phased.consensus[1] == three, _text(phased.consensus[1]))


def test_one_repeat_unit_splits_a_haplotype_only_with_enough_reads() raises:
    """Clusters apart only by one TG of a (TG)10 / (TG)11 repeat: one cluster with the smaller under 30% of the
    reads (HiFi reads' slips there), two at half."""
    var left = List[UInt8]("ACCTGAGCAATCGGATCCAGTTACGGA".as_bytes())
    var right = List[UInt8]("CCATTGAGCGTACAGGTCAATGCCTAG".as_bytes())
    var ten = left.copy()
    var eleven = left.copy()
    for k in range(11):
        if k < 10:
            ten.extend(List[UInt8]("TG".as_bytes()))
        eleven.extend(List[UInt8]("TG".as_bytes()))
    ten.extend(right.copy())
    eleven.extend(right.copy())
    for minority in [4, 8]:
        var reads = List[List[UInt8]]()
        for k in range(16 - minority):
            reads.append(ten.copy())
        for k in range(minority):
            reads.append(eleven.copy())
        var got = clustered_consensus(reads^, max_clusters=2, min_freq=0.2, options=_unpinned())
        assert_equal(len(got.consensus), 1 if minority == 4 else 2, String(t"{minority} of 16 reads one TG longer"))


def test_a_minority_haplotype_backed_by_its_labels_keeps_its_letter() raises:
    """NA18941's CYP1A1 in panStar: 8 reads of haplotype 0 and 3 of haplotype 1 over a block, apart at one SNV every
    labelled read agrees on, and 4 more reads of another phase set (unlabelled) carrying haplotype 0's letter. Three
    of fifteen is under a quarter of the reads, but each spelling is carried by three of its own haplotype's
    labelled reads and none of the other's: a difference, not an error. Counted as one, the unlabelled reads go to
    haplotype 0 alone; uncounted, they went to both and outvoted haplotype 1's own three reads."""
    seed(41)
    var first = _random(700)
    var second = first.copy()
    var letters = "ACGT".as_bytes()
    second[350] = letters[0] if first[350] != letters[0] else letters[1]
    var reads = List[List[UInt8]]()
    var haplotypes = List[Int]()
    for _ in range(8):
        reads.append(first.copy())
        haplotypes.append(0)
    for _ in range(3):
        reads.append(second.copy())
        haplotypes.append(1)
    for _ in range(4):
        reads.append(first.copy())
        haplotypes.append(-1)
    var phased = phased_consensus(reads^, haplotypes)
    assert_equal(_text(phased.consensus[0]), _text(first))
    assert_equal(_text(phased.consensus[1]), _text(second), "haplotype 1 keeps its own reads' letter")


def test_arguments_are_checked() raises:
    var reads = List[List[UInt8]]()
    reads.append(_random(40))
    with assert_raises(contains="max_clusters"):
        _ = clustered_consensus(reads.copy(), max_clusters=0)
    with assert_raises(contains="min_freq"):
        _ = clustered_consensus(reads.copy(), min_freq=1.5)
    with assert_raises(contains="anchors for"):
        _ = clustered_consensus(reads.copy(), [BOTH, BOTH])
    var none = clustered_consensus(List[List[UInt8]]())
    assert_equal(len(none.consensus), 1)


def _swapped(seq: List[UInt8], at: List[Int]) -> List[UInt8]:
    """`seq` with the letter at each of `at` turned to another."""
    var out = seq.copy()
    var letters = "ACGT".as_bytes()
    for i in at:
        out[i] = letters[(Int(String("ACGT").find(chr(Int(out[i])))) + 1) % 4]
    return out^


def _slice(seq: List[UInt8], start: Int, end: Int) -> List[UInt8]:
    var out = List[UInt8](capacity=end - start)
    for i in range(start, end):
        out.append(seq[i])
    return out^


def _labelled_right_unlabelled_left(one: List[UInt8], two: List[UInt8]) raises -> List[List[UInt8]]:
    """`phased_consensus` of six reads of each haplotype labelled and anchored right, over letters 480 on, and six of
    each unlabelled and anchored left, over letters 0 to 720: no read spans, and only unlabelled reads reach the
    left end, as where pgphase's phase set ends inside a window (HG002's HLA-DPA1 thinned to 10x and 15x)."""
    var reads = List[List[UInt8]]()
    var anchors = List[Int]()
    var haplotypes = List[Int]()
    for h in range(2):
        var seq = one.copy() if h == 0 else two.copy()
        for _ in range(6):
            reads.append(_slice(seq, 480, len(seq)))
            anchors.append(RIGHT)
            haplotypes.append(h)
        for _ in range(6):
            reads.append(_slice(seq, 0, 720))
            anchors.append(LEFT)
            haplotypes.append(-1)
    return phased_consensus(reads^, haplotypes, anchors=anchors).consensus.copy()


def test_unlabelled_reads_build_the_haplotype_whose_letters_they_carry() raises:
    """The haplotypes differ at letters 240, 600 and 960; the unlabelled left reads carry 600, where the labelled
    consensus sequences differ too, and so each joins its own haplotype, which it alone carries to the left end and
    to letter 240. Built from the labelled reads alone, each haplotype began at letter 480."""
    seed(91)
    var one = _random(1200)
    var two = _swapped(one, [240, 600, 960])
    var got = _labelled_right_unlabelled_left(one, two)
    assert_true(got[0] == one, _text(got[0]))
    assert_true(got[1] == two, _text(got[1]))


def test_unlabelled_reads_where_the_haplotypes_read_alike_build_both() raises:
    """The haplotypes differ at letter 960 alone, which no unlabelled read reaches: each such read covers no place
    where the two differ and builds both, each then whole from its left end."""
    seed(92)
    var one = _random(1200)
    var two = _swapped(one, [960])
    var got = _labelled_right_unlabelled_left(one, two)
    assert_true(got[0] == one, _text(got[0]))
    assert_true(got[1] == two, _text(got[1]))


def test_unlabelled_reads_join_haplotypes_of_different_lengths() raises:
    """The second haplotype lacks 30 letters (590 to 620) and differs at 240, 900 and 960: the two consensus sequences
    differ in length, each labelled read is placed on its own haplotype (the shorter's at each place's span there)
    and the unlabelled left reads, which reach the deletion, still join theirs."""
    seed(93)
    var one = _random(1200)
    var swapped = _swapped(one, [240, 900, 960])
    var two = _slice(swapped, 0, 590)
    two.extend(_slice(swapped, 620, len(swapped)))
    var got = _labelled_right_unlabelled_left(one, two)
    assert_true(got[0] == one, _text(got[0]))
    assert_true(got[1] == two, _text(got[1]))


def test_a_haplotype_inside_a_read_has_each_letters_read_letter() raises:
    # `_haplotype_at`: a haplotype placed inside a longer read gives each haplotype letter the read letter it stands
    # before; a read letter the haplotype lacks counts before the haplotype letter after it, as `_query_at` counts.
    var haplotype = List[UInt8]("GATTACAGGTCCATTAGC".as_bytes())
    var flanked = List[UInt8]("CCCC".as_bytes())
    flanked.extend(haplotype.copy())
    flanked.extend(List[UInt8]("TTTTT".as_bytes()))
    var aligner = Aligner()
    var at = _haplotype_at(aligner.align(Span(flanked), Span(haplotype), Costs.edit(), Mode.INFIX, ties=Ties.LEFT))
    assert_equal(len(at), len(haplotype) + 1)
    for p in range(len(haplotype) + 1):
        assert_equal(at[p], 4 + p)
    # One letter more in the read after the haplotype's eighth: from there on, one further along.
    var longer = List[UInt8]("CCCC".as_bytes())
    longer.extend(List[UInt8](Span(haplotype)[:8]))
    longer.append(UInt8(ord("A")))
    longer.extend(List[UInt8](Span(haplotype)[8:]))
    longer.extend(List[UInt8]("TTTTT".as_bytes()))
    at = _haplotype_at(aligner.align(Span(longer), Span(haplotype), Costs.edit(), Mode.INFIX, ties=Ties.LEFT))
    for p in range(len(haplotype) + 1):
        assert_equal(at[p], 4 + p + (1 if p >= 8 else 0), String(t"letter {p}"))


def _pgphase_options() -> ConsensusOptions:
    """pgphase's plegma options for its MSA regions (msa.mojo `_options`): unpinned, its scoring."""
    var options = _unpinned()
    options.poa.scoring.match_score = 2
    options.poa.scoring.mismatch = 6
    options.poa.scoring.open1 = 6
    options.poa.scoring.extend1 = 2
    options.poa.scoring.open2 = 24
    options.poa.scoring.extend2 = 1
    options.two_lengths_db = 8.0
    return options^


def test_a_snp_in_a_run_keeps_its_reads() raises:
    """HG002 chr20:56,892,617-56,892,658 (a pgphase region no read's phase set held): a G at a poly-T's start on 28 of
    62 reads, the run a letter or two short or long on many. Moved by alignment cost at the scoring's costs (a gap
    letter dearer than a mismatch), reads went to the consensus of their run length: the G cluster's consensus lost
    its G, the two were apart by a run alone and merged, the G written on both. Moved as the first split weighs
    differences, the 28 stay together."""
    var reads = List[List[UInt8]]()
    var anchors = List[Int]()
    var strands = List[Int]()
    with open("tests/data/snp_in_run/reads.txt", "r") as f:
        for line in f.read().split("\n"):
            if line.byte_length() == 0:
                continue
            var parts = line.split(":")
            var seq = List[UInt8]()
            for b in parts[0].as_bytes():
                seq.append(b)
            reads.append(seq^)
            anchors.append(Int(parts[1]))
            strands.append(Int(parts[2]))
    var with_g = 0
    for r in reads:
        with_g += Int("CCTTG" in _text(r))
    assert_equal(with_g, 28)
    var got = clustered_consensus(reads^, anchors, strands, 2, Float64(Float32(0.12)), _pgphase_options())
    assert_equal(len(got.consensus), 2)
    var g = 0 if "CCTTG" in _text(got.consensus[0]) else 1
    assert_true("CCTTG" in _text(got.consensus[g]), _text(got.consensus[g]))
    assert_true(not ("CCTTG" in _text(got.consensus[1 - g])), _text(got.consensus[1 - g]))
    assert_true(len(got.members[g]) >= 26, String(t"{len(got.members[g])} reads with the G"))


def _pairs(a: String, b: String) raises -> List[Tuple[Int, Int]]:
    var aligner = Aligner()
    var runs_a = List[Int]()
    var runs_b = List[Int]()
    _run_pairs(aligner, List[UInt8](a.as_bytes()), List[UInt8](b.as_bytes()), 3, runs_a, runs_b)
    var out = List[Tuple[Int, Int]]()
    for k in range(len(runs_a)):
        out.append((runs_a[k], runs_b[k]))
    return out^


def test_run_pairs_finds_each_run_difference_and_only_those() raises:
    # A run a letter apart: the run's start in each.
    var one = _pairs("ACGTTTTTTACG", "ACGTTTTTACG")
    assert_equal(len(one), 1)
    assert_equal(one[0][0], 3)
    assert_equal(one[0][1], 3)
    # Beside a SNP elsewhere, the run is still found (the two apart by more than runs).
    var with_snp = _pairs("ACCGTTTTTTAGC", "ACAGTTTTTAGC")
    assert_equal(len(with_snp), 1)
    assert_equal(with_snp[0][0], 4)
    # One unit of a (TG)n apart: no run (HG002 chr20:21,208,762; at unit costs the unit split into a T of the TTT
    # before the repeat and a G, and that T read as a run's letter).
    var unit = _pairs("CTCCTGGAATTTGTGTGTGTGTGTGTGTGTGGCAGCAG", "CTCCTGGAATTTGTGTGTGTGTGTGTGTGGCAGCAG")
    assert_equal(len(unit), 0)
    # A substitution splitting a run: no pair of runs whose lengths differ by more than the gap.
    var split = _pairs("CTTTTGTTTTC", "CTTTTTTTTTTC")
    assert_equal(len(split), 0)
    # Either way round: a one-letter run against three is no run difference (the shorter under noise_run - 1).
    assert_equal(len(_pairs("ACTGA", "ACTTTGA")), 0)
    assert_equal(len(_pairs("ACTTTGA", "ACTGA")), 0)
    # One gap at the meeting of two runs (AAAA TTTT against AAA TTT: the alignment's one gap AT) is a letter of each.
    var junction = _pairs("CAAAATTTTG", "CAAATTTG")
    assert_equal(len(junction), 2)
    assert_equal(junction[0][0], 1)
    assert_equal(junction[1][0], 5)
    assert_equal(junction[1][1], 4)
    var aligner = Aligner()
    var ra = List[Int]()
    var rb = List[Int]()
    assert_true(
        _run_pairs(aligner, List[UInt8]("CAAAATTTTG".as_bytes()), List[UInt8]("CAAATTTG".as_bytes()), 3, ra, rb)
    )
    assert_true(
        not _run_pairs(
            aligner, List[UInt8]("CCGTTTTTTAGC".as_bytes()), List[UInt8]("CAGTTTTTAGC".as_bytes()), 3, ra, rb
        )
    )


def test_run_pairs_of_sequences_too_far_apart_is_none() raises:
    # Two haplotypes kilobases apart (a VNTR expansion, alleles hundreds of edits apart): an alignment past the
    # aligner's memory is no runs, not an error out of the consensus.
    seed(7)
    var a = _random(20000)
    var b = List[UInt8]()
    for k in range(len(a)):
        if k >= 6000 and k < 11000:
            continue
        b.append(a[k] if k % 50 != 0 else (UInt8(ord("A")) if a[k] != UInt8(ord("A")) else UInt8(ord("C"))))
    var aligner = Aligner()
    var ra = List[Int]()
    var rb = List[Int]()
    _ = _run_pairs(aligner, a, b, 3, ra, rb)


def _measured(cons: String, rd: String, s: Int, e: Int) raises -> Int:
    var c = List[UInt8](cons.as_bytes())
    var r = List[UInt8](rd.as_bytes())
    var aligner = Aligner()
    var found = aligner.align(Span(c), Span(r), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
    return _length_in(found.cigar, found.reference_start, r, s, e, c[s], c)


def test_a_read_run_is_its_own_between_the_flank_letters() raises:
    assert_equal(_measured("CGTAAAATGC", "CGTAAAAATGC", 3, 7), 5)
    assert_equal(_measured("CGTAAAATGC", "CGTAAATGC", 3, 7), 3)
    # A repeat unit's gap placed leftmost takes the run's last letter: the read's own run is still A4.
    assert_equal(_measured("CAAAATATATAG", "CAAAATATAG", 1, 5), 4)
    assert_equal(_measured("CAAAATATATAG", "CAAATATAG", 1, 5), 3)
    # Another letter inside the run: no vote.
    assert_equal(_measured("CGTAAAATGC", "CGTAAGAATGC", 3, 7), -1)
    # A read without the letter between two runs of one letter: its run ends at another letter than the flank's.
    assert_equal(_measured("GATTTTCTTTTTTGA", "GATTTTTTTTTTGA", 2, 6), -1)
    assert_equal(_measured("GATTTTCTTTTTTGA", "GATTTTTTTTTTGA", 7, 13), -1)


def test_phased_with_no_labelled_read_builds_from_the_unlabelled() raises:
    """Reads of neither haplotype alone make both haplotypes their one consensus, never two empty ones."""
    seed(211)
    var truth = _random(300)
    var reads = List[List[UInt8]]()
    for _ in range(6):
        reads.append(_mutated(truth, 0.01))
    var got = phased_consensus(reads.copy(), List[Int](length=6, fill=-1))
    assert_equal(len(got.consensus), 2)
    assert_true(distance(_text(got.consensus[0]), _text(truth)) <= 1)
    assert_true(got.consensus[0] == got.consensus[1])


def test_phased_haplotype_without_labels_is_its_unlabelled_reads() raises:
    """A phase set covering one haplotype only: 8 reads labelled 0, the other allele's 5 reads unlabelled. Haplotype 1
    is made from the unlabelled reads, so the heterozygote stays (taking haplotype 0's consensus for it would lose
    it)."""
    seed(223)
    var a = _random(400)
    var b = a.copy()
    var letters = "ACGT".as_bytes()
    for i in [60, 170, 290]:
        b[i] = letters[0] if b[i] != letters[0] else letters[1]
    var reads = List[List[UInt8]]()
    var labels = List[Int]()
    for k in range(13):
        reads.append(a.copy() if k < 8 else b.copy())
        labels.append(0 if k < 8 else -1)
    var options = ConsensusOptions()
    options.poa.seeded = False
    var got = phased_consensus(reads^, labels, options=options)
    assert_true(got.consensus[0] == a)
    assert_true(got.consensus[1] == b)


def test_phased_reads_all_inside_build_as_spanning() raises:
    """As `clustered_consensus` has it: with every read inside, all of them build as if spanning."""
    seed(227)
    var truth = _random(150)
    var reads = List[List[UInt8]]()
    for _ in range(6):
        reads.append(truth.copy())
    var labels: List[Int] = [0, 0, 0, 1, 1, 1]
    var got = phased_consensus(reads^, labels, List[Int](length=6, fill=INSIDE))
    assert_true(got.consensus[0] == truth)
    assert_true(got.consensus[1] == truth)


def test_phased_swapping_the_labels_swaps_the_haplotypes() raises:
    """Labels 0 and 1 swapped give the two haplotypes swapped: with labelled-only consensus sequences of one length,
    the longer's role (where `_differences` puts the places, which haplotype an unlabelled read is placed on) went
    to haplotype 0 by its label, and a LEFT read built the other haplotype once swapped."""
    var reads: List[List[UInt8]] = [
        _b(
            "TACTGTGGATTGGTACAGGATAATAATAATAATAAAAAAAAAAAATCATGGCAGTCAATCCTGTCCCGGGACTTACTTACTTACTTACTTACTTACTTACTTGGGGCTAACATCCTCACGAAGTTATTTATAATGACGC"
        ),
        _b(
            "GATAATAATAACAATAAAAAAAAAAAAAAATCATGGCAGTCAATCCTGTCCCGGGACTTACTTACTTACTCACTTACTTACTTGGGGCTAACATCCTCACGAAGTTATTTA"
        ),
        _b(
            "TGGTAAACAATCACTCCTATACGAAGAAGCTACTGTGGATTGGTACAGGATAATAATAATAATAAAAAAAAAAAAAAATCATGGCAGTCAATCCTGTCCCGGGACTTACTTACTTACTCACTTACTTACTTGGGGCTAACATCCTTCACGAAGTTATTTATAATGACGC"
        ),
        _b(
            "GGTACAGGATAATAATAATAATAAAAAAAAAAAATCCTGGCAGTCAATCCCTGTCCCGGGACTTACTTACTTACTTACTTACTTACTTACTTGGGGCTAACATCCTCACGAAGTTATTTATAATGACGC"
        ),
        _b("TGTAAACAATCACTCCTATACGAAGAAGCTACTGTGGATTGGTACAGGATAATAATAATAATAAAAAAAAAAAATCA"),
        _b(
            "TAATAATAAAAAAAAAAAAAAATCATGGCAGTCAATCCTGGCCCGGGACTTACTTACTTACTCACTTACTTACTTGGGGCTAACATCCTCACGAAGTTATTTATAATGACGC"
        ),
        _b("TACTAATAATAAAAAAAAAAAAAAATCATGGCAGTCAATCCTGTCCCGGGACTTACTTACTTACTCACTT"),
        _b(
            "TGTAAACAATCACTCCTATACGAAGAAGCTACTGTGGATTGGTACAGGATAATAATAATAATAAAAAAAAAAAATCATGGCAGTCAATCCTGTCCCGGGACTTACTTACTTACTTACTTACTTACTTACTTGGGGCTAACATCCTCACGAAGTTATTTATAATGACGC"
        ),
        _b(
            "TGTAAACAATCACTCCTATACGAAGAAGTACTGTGGATTGGTACAGGATAATAATAATAATAAAAAAAAAAAAAAATCATGGCAGTCAGATCCTGTCCCGGGACTTACTTACTTACTCACTTACTTACTTGGGGCTAACATCCTCACGAAGTTATTTATAAATGACGC"
        ),
        _b("AAAATCATGGCAGTCAATCCTGTCCCGGGACTTACTTACTTACTTACTTACTTACTTACTTGGGGCTAACATCCTCACGAAGGTTATTTATAATGACGC"),
    ]
    var labels: List[Int] = [0, 1, 1, 0, -1, 1, -1, 0, 1, 0]
    var swapped = List[Int]()
    for h in labels:
        swapped.append(1 - h if h >= 0 else -1)
    var anchors: List[Int] = [2, 3, 0, 2, 1, 2, 3, 0, 0, 2]
    var a = phased_consensus(reads.copy(), labels, anchors)
    var b = phased_consensus(reads.copy(), swapped, anchors)
    assert_true(a.consensus[0] == b.consensus[1])
    assert_true(a.consensus[1] == b.consensus[0])


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
