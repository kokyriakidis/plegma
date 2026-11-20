from std.random import random_ui64, seed
from std.testing import assert_equal, assert_true, TestSuite

from karan import distance

from plegma.graph import Graph
from plegma.kernel import BOTH, LEFT, RIGHT
from plegma.consensus import consensus
from plegma.poa import INPUT_ORDER, PoaOptions, poa
from plegma.polish import polish


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


def _b(text: StringSlice) -> List[UInt8]:
    return List[UInt8](text.as_bytes())


def _slice(seq: List[UInt8], start: Int, end: Int) -> List[UInt8]:
    var out = List[UInt8](capacity=end - start)
    for i in range(start, end):
        out.append(seq[i])
    return out^


def _assembled(
    truth: List[UInt8], n_both: Int, n_left: Int, n_right: Int, rate: Float64, overhang: Int = 0, pinned: Bool = True
) raises -> List[UInt8]:
    """The polished consensus of reads of `truth`: `n_both` spanning it, `n_left` fixed at its start, `n_right` at
    its end; with `overhang`, the one-sided reads run that many letters past the far end instead (into sequence
    beyond it), as a local assembler cuts them."""
    var seqs = List[List[UInt8]]()
    var names = List[String]()
    var anchors = List[Int]()
    var length = len(truth)
    var before = _random(overhang)
    var after = _random(overhang)
    for k in range(n_both + n_left + n_right):
        var anchor = BOTH if k < n_both else (LEFT if k < n_both + n_left else RIGHT)
        var start = 0
        var end = length
        if anchor == LEFT:
            end = Int(random_ui64(UInt64(length * 3 // 10), UInt64(length * 95 // 100)))
        elif anchor == RIGHT:
            start = Int(random_ui64(UInt64(length // 20), UInt64(length * 7 // 10)))
        var piece = _slice(truth, start, end)
        if overhang > 0 and anchor == LEFT:
            piece = truth + after.copy()
        elif overhang > 0 and anchor == RIGHT:
            piece = before + truth.copy()
        seqs.append(_mutated(piece, rate))
        names.append(String(t"r{k}"))
        anchors.append(anchor)
    var options = PoaOptions()
    options.seeded = pinned
    var made = poa(names, seqs, anchors, options, Graph())
    # Every read is in the graph and its row spells it.
    var rows = made.graph.msa()
    for k in range(len(seqs)):
        var spelled = List[UInt8]()
        for c in rows[made.shown[k]].as_bytes():
            if c != UInt8(ord("-")):
                spelled.append(c)
        assert_true(spelled == seqs[k])
    assert_equal(len(made.consensus), 1)
    var everyone = List[Int]()
    for k in range(len(seqs)):
        everyone.append(k)
    return polish(made.consensus[0], seqs, everyone, 50, anchors=anchors)


def test_few_spanning_reads_assemble_the_truth() raises:
    seed(89)
    var truth = _random(3000)
    var consensus = _assembled(truth, 1, 8, 8, 0.05)
    assert_true(distance(String(unsafe_from_utf8=consensus), String(unsafe_from_utf8=truth)) <= 3)


def test_no_spanning_read_joins_the_two_sides() raises:
    seed(97)
    var truth = _random(3000)
    var consensus = _assembled(truth, 0, 10, 10, 0.05)
    assert_true(distance(String(unsafe_from_utf8=consensus), String(unsafe_from_utf8=truth)) <= 3)


def test_reads_past_the_far_anchor_do_not_extend_the_consensus() raises:
    # More one-sided reads than spanning ones, each running past the far anchor: the consensus still ends where the
    # spanning reads do.
    seed(103)
    var truth = _random(3000)
    var consensus = _assembled(truth, 3, 8, 8, 0.05, overhang=330)
    assert_true(distance(String(unsafe_from_utf8=consensus), String(unsafe_from_utf8=truth)) <= 3)


def test_reads_past_the_far_anchor_unpinned() raises:
    # The same, each read aligned to the whole graph: the left reads' letters past the end join the graph after it,
    # and the right reads still end at the right anchor, not at the end of those letters.
    seed(103)
    var truth = _random(3000)
    var consensus = _assembled(truth, 3, 8, 8, 0.05, overhang=330, pinned=False)
    assert_true(distance(String(unsafe_from_utf8=consensus), String(unsafe_from_utf8=truth)) <= 3)


def test_another_copys_one_sided_reads_get_no_vote() raises:
    # Five spanning reads of the truth, two left reads of it, and eight left reads of another copy (2% of its letters
    # substituted), outvoting them: the consensus is the truth's all the same.
    seed(107)
    var truth = _random(2000)
    var other = truth.copy()
    var letters = "ACGT".as_bytes()
    for i in range(0, len(other), 50):
        other[i] = letters[(Int(random_ui64(1, 3)) + Int(_code(other[i]))) % 4]
    var seqs = List[List[UInt8]]()
    var names = List[String]()
    var anchors = List[Int]()
    for k in range(15):
        if k >= 5 and k < 13:
            seqs.append(_mutated(other, 0.01))
        else:
            seqs.append(_mutated(truth, 0.01))
        names.append(String(t"r{k}"))
        anchors.append(BOTH if k < 5 else LEFT)
    var made = poa(names, seqs, anchors, PoaOptions(), Graph())
    assert_equal(len(made.foreign), 8)
    for r in made.foreign:
        assert_true(r >= 5 and r < 13)
    assert_true(distance(String(unsafe_from_utf8=made.consensus[0]), String(unsafe_from_utf8=truth)) <= 3)


def test_a_spanning_read_wins_a_tie() raises:
    # One read on both anchors, one fixed at the right anchor alone that differs from it at three places (a paralog's,
    # missing the left anchor): one vote each, and the read on both anchors is the surer.
    seed(109)
    var truth = _random(400)
    var other = truth.copy()
    var letters = "ACGT".as_bytes()
    for i in [250, 300, 350]:
        other[i] = letters[(_code(other[i]) + 1) % 4]
    var seqs: List[List[UInt8]] = [truth.copy(), _slice(other, 150, 400)]
    var names: List[String] = ["both", "right"]
    var anchors: List[Int] = [BOTH, RIGHT]
    var made = poa(names, seqs, anchors, PoaOptions(), Graph())
    assert_true(made.consensus[0] == truth)


def _code(letter: UInt8) -> Int:
    if letter == UInt8(ord("A")):
        return 0
    if letter == UInt8(ord("C")):
        return 1
    return 2 if letter == UInt8(ord("G")) else 3


def test_one_side_assembles_a_prefix_where_reads_agree() raises:
    seed(101)
    var truth = _random(3000)
    var consensus = _assembled(truth, 0, 12, 0, 0.05)
    # The consensus stops where fewer than three reads go on: a prefix of the truth, most of it.
    assert_true(len(consensus) >= 1000)
    var prefix = _slice(truth, 0, len(consensus))
    assert_true(distance(String(unsafe_from_utf8=consensus), String(unsafe_from_utf8=prefix)) <= 10)


def _fasta(path: String) raises -> List[Tuple[String, List[UInt8]]]:
    """A FASTA file's records: each header line (without its `>`) and its letters."""
    var out = List[Tuple[String, List[UInt8]]]()
    var header = String()
    var letters = List[UInt8]()
    var any = False
    with open(path, "r") as f:
        for line in f.read().splitlines():
            if line.startswith(">"):
                if any:
                    out.append((header, letters^))
                    letters = List[UInt8]()
                header = String(line.removeprefix(">"))
                any = True
            else:
                for byte in line.as_bytes():
                    letters.append(byte)
    if any:
        out.append((header, letters^))
    return out^


def test_hg00512_drb1_joins_where_both_sides_have_reads() raises:
    """Real reads: HG00512's HLA-DRB1 haplotype 1 (DRB1*08:03) as panStar windows it, 9 reads anchored left and 18
    right, none spanning; the right side reaches the left end on one read alone. Its consensus must be the haplotype
    HG00512's own assembly (HGSVC3 verkko) spells there. Joined at the middle of the pins it was 13 letters from it,
    9 insertions and 2 deletions of that one read over 1.1 kb where 9 left reads disagreed, and the allele was read
    *08:01:05."""
    var seqs = List[List[UInt8]]()
    var anchors = List[Int]()
    var strands = List[Int]()
    for record in _fasta("tests/data/hg00512_drb1/reads.fa"):
        var fields = record[0].split(" ")
        anchors.append(LEFT if fields[1] == "anchor=left" else RIGHT)
        strands.append(Int(fields[2].removeprefix("strand=")))
        seqs.append(record[1].copy())
    assert_equal(len(seqs), 27)
    var truth = _fasta("tests/data/hg00512_drb1/truth.fa")[0][1].copy()
    var made = consensus(seqs^, anchors, strands)
    assert_equal(len(made), 1)
    var edits = distance(String(unsafe_from_utf8=made[0]), String(unsafe_from_utf8=truth))
    assert_equal(edits, 0, String(t"{edits} edits from HG00512's assembly"))


def test_no_reads_make_an_empty_poa() raises:
    var names = List[String]()
    var seqs = List[List[UInt8]]()
    var made = poa(names, seqs, List[Int](), PoaOptions(), Graph())
    assert_equal(len(made.graph), 0)
    assert_equal(len(made.consensus), 0)


def test_a_right_read_does_not_end_on_a_left_reads_overhang() raises:
    """Pinned or not, a right read ends among the nodes the spanning reads' end is reached from, never on the
    letters a left read carries past the right anchor."""
    for pinned in [True, False]:
        var names: List[String] = ["s1", "s2", "s3", "l", "r"]
        var seqs: List[List[UInt8]] = [_b("AC"), _b("AC"), _b("AC"), _b("ACGTACGTAA"), _b("TTGCA")]
        var anchors: List[Int] = [BOTH, BOTH, BOTH, LEFT, RIGHT]
        var options = PoaOptions()
        options.seeded = pinned
        var made = poa(names, seqs, anchors, options, Graph())
        ref graph = made.graph
        ref left = graph.paths[made.shown[3]]
        var spanning_end = graph.paths[made.shown[0]][1]
        var past = False
        var overhang = List[Int32]()
        for node in left:
            if past:
                overhang.append(node)
            past = past or node == spanning_end
        for node in graph.paths[made.shown[4]]:
            assert_true(node not in overhang)


def test_unjoined_right_reads_align_to_the_right_side() raises:
    """Two sides overlapping by 150 letters, too few to join (`JOIN_LEAST`), stay apart; every right read is aligned
    to the right side's part of the graph, never fused onto the left side's path where the two overlap."""
    seed(131)
    var truth = _random(900)
    var names = List[String]()
    var seqs = List[List[UInt8]]()
    var anchors = List[Int]()
    for k in range(12):
        var left = k < 6
        seqs.append(_mutated(_slice(truth, 0, 500) if left else _slice(truth, 350, 900), 0.02))
        names.append(String(t"r{k}"))
        anchors.append(LEFT if left else RIGHT)
    for pinned in [True, False]:
        var options = PoaOptions()
        options.seeded = pinned
        var copy = seqs.copy()
        var made = poa(names, copy, anchors, options, Graph())
        assert_equal(len(made.consensus), 2)
        ref graph = made.graph
        var left_nodes = List[Int32]()
        for k in range(6):
            for node in graph.paths[made.shown[k]]:
                left_nodes.append(node)
        for k in range(6, 12):
            for node in graph.paths[made.shown[k]]:
                assert_true(node not in left_nodes)


def test_scores_scaled_up_keep_the_msa() raises:
    """Every score a hundredfold is the same optimum, so the same MSA: the band's 16-bit rows then cannot hold a
    row's spread, and the kernel must widen its lanes rather than wrap."""
    seed(137)
    for _ in range(3):
        var truth = _random(3000)
        var names = List[String]()
        var seqs = List[List[UInt8]]()
        var anchors = List[Int]()
        for k in range(6):
            seqs.append(_mutated(truth, 0.1))
            names.append(String(t"r{k}"))
            anchors.append(BOTH)
        var options = PoaOptions()
        options.seeded = False
        options.order = INPUT_ORDER
        var scaled = options.copy()
        scaled.scoring.match_score *= 100
        scaled.scoring.mismatch *= 100
        scaled.scoring.open1 *= 100
        scaled.scoring.open2 *= 100
        scaled.scoring.extend1 *= 100
        scaled.scoring.extend2 *= 100
        var a = poa(names, seqs, anchors, options, Graph()).graph.msa()
        var b = poa(names, seqs, anchors, scaled, Graph()).graph.msa()
        for k in range(len(a)):
            assert_true(a[k] == b[k])


def test_sides_overlapping_wholly_are_joined() raises:
    """With no spanning read, left and right reads that all reach the far anchor (the two sides' consensus sequences
    alike end to end), or a right read reaching a few letters past the left anchor, are one assembly: pins with right
    starting before left, or ending where left does, were refused and the two sides came out apart."""
    seed(149)
    var truth = _random(3000)
    for overhang in [0, 8]:
        var names = List[String]()
        var seqs = List[List[UInt8]]()
        var anchors = List[Int]()
        for k in range(12):
            var left = k < 6
            var piece = truth.copy()
            if not left and k == 6 and overhang > 0:
                piece = _random(overhang) + truth.copy()
            seqs.append(_mutated(piece, 0.01))
            names.append(String(t"r{k}"))
            anchors.append(LEFT if left else RIGHT)
        var made = poa(names, seqs, anchors, PoaOptions(), Graph())
        assert_equal(len(made.consensus), 1)
        assert_true(distance(String(unsafe_from_utf8=made.consensus[0]), String(unsafe_from_utf8=truth)) <= 10)


def test_a_join_inside_a_tandem_array_is_never_a_wrong_copy() raises:
    """A 5 kb local assembly whose two sides overlap by 500 letters wholly inside a 25-copy array of a 40-letter unit
    (copies 1% apart), no read spanning, reads at 1% errors: one 19-mer pin unique in both sides was another copy's,
    and the join there dropped 240 letters. Every copy number is now a candidate, scored over the whole overlap; this
    one's true placement and the one a copy off differ by two letters, too few to tell. Either the assembly is the
    truth's length or its two sides stay apart: never a copy dropped or doubled."""
    var seqs = List[List[UInt8]]()
    var anchors = List[Int]()
    for record in _fasta("tests/data/tandem_join/reads.fa"):
        anchors.append(LEFT if "anchor=left" in record[0] else (RIGHT if "anchor=right" in record[0] else BOTH))
        seqs.append(record[1].copy())
    var truth = _fasta("tests/data/tandem_join/truth.fa")[0][1].copy()
    var made = consensus(seqs^, anchors)
    if len(made) == 1:
        assert_true(distance(String(unsafe_from_utf8=made[0]), String(unsafe_from_utf8=truth)) <= 10)
    else:
        assert_equal(len(made), 2)


def test_a_read_ending_on_another_letter_still_reaches_there() raises:
    """Three one-letter left reads, C, A and C: A is in C's column, not on the consensus path, and counting only the
    path's nodes it reached nothing, so the one-sided cut (where the last of fewer than -c reads stops) emptied the
    consensus."""
    var names: List[String] = ["a", "b", "c"]
    var seqs: List[List[UInt8]] = [_b("C"), _b("A"), _b("C")]
    var anchors: List[Int] = [LEFT, LEFT, LEFT]
    var made = poa(names, seqs, anchors, PoaOptions(), Graph())
    assert_equal(len(made.consensus), 1)
    assert_equal(len(made.consensus[0]), 1)


def test_the_consensus_does_not_depend_on_read_order() raises:
    """The guide order goes by the reads, not their input places: a read set and the same reads reversed or rotated
    give one consensus (abPOA's order, the earlier read of the most similar pair first, did not: 3 of 150 sets of 3 to
    6 reads at 8% errors)."""
    seed(151)
    for _ in range(40):
        var truth = _random(400)
        var n = 3 + Int(random_ui64(0, 3))
        var reads = List[List[UInt8]]()
        for _ in range(n):
            reads.append(_mutated(truth, 0.08))
        var reversed = List[List[UInt8]]()
        var rotated = List[List[UInt8]]()
        for k in range(n):
            reversed.append(reads[n - 1 - k].copy())
            rotated.append(reads[(k + 1) % n].copy())
        var want = consensus(reads^)
        assert_true(consensus(reversed^) == want)
        assert_true(consensus(rotated^) == want)


def test_a_one_sided_read_out_of_band_is_aligned_again_wider() raises:
    """A 1 kb local assembly at 10% errors, its two sides overlapping by 200 letters: a right read aligned to the
    right side's graph lost its path in the band (its diagonal follows the graph's heaviest remaining path) and scored
    far under a third of a perfect score; never retried (an end is always reached), it left the right side 18 letters
    long and the two sides apart. Aligned again with the band doubled, the assembly is one sequence near the truth."""
    var seqs = List[List[UInt8]]()
    var anchors = List[Int]()
    for record in _fasta("tests/data/band_retry/reads.fa"):
        anchors.append(LEFT if "anchor=left" in record[0] else RIGHT)
        seqs.append(record[1].copy())
    var truth = _fasta("tests/data/band_retry/truth.fa")[0][1].copy()
    var made = consensus(seqs^, anchors)
    assert_equal(len(made), 1)
    assert_true(distance(String(unsafe_from_utf8=made[0]), String(unsafe_from_utf8=truth)) <= 5)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
