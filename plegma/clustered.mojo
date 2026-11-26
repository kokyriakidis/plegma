"""Several consensus sequences of one read set: the reads clustered into haplotypes (abPOA's `-d`, as a diploid
caller asks for), each cluster's consensus polished against its own reads, and every read aligned to its cluster's.

The clusters start as `cluster_rows` splits the spanning reads by the MSA columns where they disagree; then, as
k-means does, each read goes to the consensus it differs from least (`MOVING_COSTS`: a substitution weighs more than
a letter slipped in a run) and each cluster's consensus is made again
from its reads, until no read moves. Two clusters stay apart only while each holds at least `min_freq` of the reads
and their consensus sequences differ by more than homopolymer run lengths the reads fit as one (`_two_lengths`: the
reads' length errors there split one haplotype's reads as readily as two haplotypes'), nor only by one unit of a
short tandem repeat with the smaller under `MIN_UNIT_SPLIT_FREQ` of the reads. Every run two clusters differ by is
set so, whatever else they differ by: one length in both where the reads fit one. A read covering neither anchor
(`INSIDE`) builds no consensus; it is aligned inside the one it fits best.
"""

from std.math import log10

from karan import Aligner, Alignment, Costs, Mode, Ties

from .clusters import cluster_rows
from .consensus import ConsensusOptions, ConsensusScratch, consensus
from .graph import Graph
from .alphabet import encode_base
from .homopolymer import FLANK_UNKNOWN, STRAND_UNKNOWN, HomopolymerModel, reads_say_one
from .kernel import BOTH, INSIDE, LEFT, RIGHT, Scoring
from .poa import poa_with_scratch


comptime RUN_PAIR_COSTS = Costs(4, 6, 2, -1, 0)
"""`_run_pairs`'s alignment of two consensus sequences: gap-affine, a substitution 4, a gap of k letters 6 + 2k
(WFA2-lib's default), so an indel of several letters stays one gap. It only finds the runs the two differ by (the
records are written later, by the caller's alignment to its reference); minimap2's asm5 (40, 78 + 7k, 162 + 3k) kept
more than karan's 80 MB of wavefronts on long pairs (HG002 chr20's VNTR regions)."""

comptime MOVING_COSTS = Costs(5, 2, 1, -1, 0)
"""The refinement's costs of a read against a consensus: a one-letter gap 3, a substitution 5, under the 6 of the
insertion and deletion that would otherwise spell it (no costs weigh a substitution more than two one-letter gaps:
the alignment then writes it as those two)."""

comptime MAX_REFINEMENTS = 4
"""Rounds of moving reads to their nearest consensus, at most."""

comptime MIN_UNIT_SPLIT_FREQ = 0.3
"""Two clusters apart by one unit of a tandem repeat of 2 to 6 letters stay apart only while the smaller holds this
share of the reads: reads' slips in short tandem repeats split one haplotype's reads by a unit. On HG002 HiFi chr12,
18 and 20 (pgphase's regions no phase set held), such splits with the smaller under 30% of the reads were 1 true het
to 25 false; at 30% or more, 67 to 11."""


@fieldwise_init
struct Placement(Copyable, Movable):
    """A read's alignment to its cluster's consensus: the `cluster`; the CIGAR (`=`, `X`, `I` a letter the read has
    and the consensus not, `D` one the consensus has and the read not), its `cost` under the scoring's costs, and
    the spans it covers, consensus letters `consensus_start ..< consensus_end` and read letters `read_start ..<
    read_end` (a read anchored at one end, or inside, covers part of the consensus)."""

    var cluster: Int
    var cigar: String
    var cost: Int
    var consensus_start: Int
    var consensus_end: Int
    var read_start: Int
    var read_end: Int


struct Clustered(Movable):
    """`clustered_consensus`'s outcome: per cluster its `consensus` and its `members` (read indices, in input
    order); per read its `placements` entry."""

    var consensus: List[List[UInt8]]
    var members: List[List[Int]]
    var placements: List[Placement]

    def __init__(out self):
        self.consensus = List[List[UInt8]]()
        self.members = List[List[Int]]()
        self.placements = List[Placement]()


def alignment_costs(scoring: Scoring) raises -> Costs:
    """`scoring` as costs, twice its scores so an odd match converts exactly: a mismatch costs 2 (mismatch +
    match), a gap letter 2 extend + match over 2 opens. Globally the least costly alignment is the best-scoring one."""
    var per_match = Int(scoring.match_score)
    var mismatch = 2 * (Int(scoring.mismatch) + per_match)
    if scoring.open2 > 0 and scoring.extend2 > 0:
        return Costs.two_piece(
            mismatch,
            2 * Int(scoring.open1),
            2 * Int(scoring.extend1) + per_match,
            2 * Int(scoring.open2),
            2 * Int(scoring.extend2) + per_match,
        )
    return Costs.affine(mismatch, 2 * Int(scoring.open1), 2 * Int(scoring.extend1) + per_match)


def _mode(anchor: Int) -> Mode:
    """How a read anchored so aligns to a consensus: whole (`BOTH`), from its start (`LEFT`, the consensus's end
    free), to its end (`RIGHT`), or anywhere inside it (`INSIDE`)."""
    if anchor == LEFT:
        return Mode.PREFIX
    if anchor == RIGHT:
        return Mode.SUFFIX
    if anchor == INSIDE:
        return Mode.INFIX
    return Mode.GLOBAL


def _aligned(
    mut aligner: Aligner, consensus: List[UInt8], seq: List[UInt8], anchor: Int, costs: Costs
) raises -> Alignment:
    return aligner.align(Span(consensus), Span(seq), costs, _mode(anchor), ties=Ties.LEFT)


def _consensus_of(
    mut scratch: ConsensusScratch,
    seqs: List[List[UInt8]],
    reads: List[Int],
    anchors: List[Int],
    strands: List[Int],
    options: ConsensusOptions,
    homopolymers: HomopolymerModel,
) raises -> List[UInt8]:
    """The consensus of `reads` alone; of two sides not joined (no read spanning), the longer."""
    var mine = List[List[UInt8]](capacity=len(reads))
    var mine_anchors = List[Int](capacity=len(reads))
    var mine_strands = List[Int](capacity=len(reads))
    for r in reads:
        mine.append(seqs[r].copy())
        mine_anchors.append(anchors[r])
        mine_strands.append(strands[r])
    var parts = consensus(scratch, mine^, mine_anchors, mine_strands, options, homopolymers)
    if len(parts) > 1 and len(parts[1]) > len(parts[0]):
        return parts[1].copy()
    return parts[0].copy()


def _apart_by_runs_alone(
    mut aligner: Aligner, a: List[UInt8], b: List[UInt8], noise_run: Int, mut runs: List[Int], mut single: Bool
) raises -> Bool:
    """Whether `a` and `b` differ only by letters added to or dropped from runs of one letter at least `noise_run`
    long (never with `noise_run` 0), each such run's start in `a` then in `runs`, and `single` whether each by one
    letter alone: what reads' length errors in homopolymers make, the most common split of one haplotype's reads in
    two (on HG002 HiFi, 446 of 466 such splits where the truth is homozygous)."""
    runs.clear()
    single = True
    if noise_run <= 0:
        return False
    var found = aligner.align(Span(a), Span(b), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
    if found.cost == 0:
        return False
    var i = 0  # in a
    var j = 0  # in b
    var count = 0
    for c in found.cigar.as_bytes():
        if c >= UInt8(ord("0")) and c <= UInt8(ord("9")):
            count = count * 10 + Int(c - UInt8(ord("0")))
            continue
        if c == UInt8(ord("=")):
            i += count
            j += count
        elif c == UInt8(ord("X")):
            return False
        elif c == UInt8(ord("D")):
            # Letters of `a`'s run `b` lacks: all of the run's letter, the run long enough even in `b`.
            for k in range(count):
                if a[i + k] != a[i]:
                    return False
            if _run_at(a, i) - count < noise_run - 1 or _run_at(a, i) < noise_run:
                return False
            runs.append(_run_start(a, i))
            single = single and count == 1
            i += count
        else:
            for k in range(count):
                if b[j + k] != b[j]:
                    return False
            # As a deletion's: the run long enough even in `a` (else `a` against `b` and `b` against `a` disagree).
            if _run_at(b, j) < noise_run or _run_at(b, j) - count < noise_run - 1:
                return False
            # The letters `b` has more: `a`'s run of them beside the gap.
            if i < len(a) and a[i] == b[j]:
                runs.append(_run_start(a, i))
            elif i > 0 and a[i - 1] == b[j]:
                runs.append(_run_start(a, i - 1))
            else:
                return False
            single = single and count == 1
            j += count
        count = 0
    return True


def _run_pairs(
    mut aligner: Aligner,
    a: List[UInt8],
    b: List[UInt8],
    noise_run: Int,
    mut runs_a: List[Int],
    mut runs_b: List[Int],
) -> Bool:
    """Every place `a` and `b` differ by letters added to or dropped from a run of one letter at least `noise_run`
    long, whatever else they differ by: the run's start in `a` in `runs_a` and the same run's start in `b` in
    `runs_b`, ascending, each run once; True when that is every difference (the two apart by runs alone, as one
    alignment says, not two). A gap at the meeting of two runs (AAAA TTTT against AAA TTT: one gap AT) is a letter of
    each. A run's length is decided apart from what else the two differ by: a run length is the reads' commonest
    error wherever it lies, and two consensus sequences apart by a short tandem repeat's units (or a SNP) beside a run
    kept each its own reads' run length untested against both haplotypes' reads. Two sequences too far apart to align
    within the aligner's memory: no runs, and not apart by runs alone."""
    runs_a.clear()
    runs_b.clear()
    if noise_run <= 0 or a == b:
        return False
    # Gaps opened at a cost, so an indel of several letters stays one gap: at unit costs a repeat unit the two differ
    # by can split into a letter dropped from a run beside the repeat and the rest further on (HG002 chr20:21,208,762,
    # one TG of (TG)n: read as one T of the TTT before it and a G, both haplotypes then cut to TT).
    var cigar: String
    try:
        cigar = aligner.align(Span(a), Span(b), RUN_PAIR_COSTS, Mode.GLOBAL, ties=Ties.LEFT).cigar
    except:
        return False
    var all_runs = True
    var i = 0  # in a
    var j = 0  # in b
    var count = 0
    for c in cigar.as_bytes():
        if c >= UInt8(ord("0")) and c <= UInt8(ord("9")):
            count = count * 10 + Int(c - UInt8(ord("0")))
            continue
        if c == UInt8(ord("=")) or c == UInt8(ord("M")):
            i += count
            j += count
        elif c == UInt8(ord("X")):
            all_runs = False
            i += count
            j += count
        elif c == UInt8(ord("D")):
            all_runs = _gap_runs(a, i, count, b, j, noise_run, runs_a, runs_b, False) and all_runs
            i += count
        else:
            all_runs = _gap_runs(b, j, count, a, i, noise_run, runs_b, runs_a, True) and all_runs
            j += count
        count = 0
    return all_runs and len(runs_a) > 0


def _gap_runs(
    longer: List[UInt8],
    at: Int,
    count: Int,
    shorter: List[UInt8],
    at_short: Int,
    noise_run: Int,
    mut runs_long: List[Int],
    mut runs_short: List[Int],
    swapped: Bool,
) -> Bool:
    """`_run_pairs`'s reading of one gap: `longer`'s letters `at ..< at + count` that `shorter` lacks before its
    letter `at_short`. One letter's block lengthens the run beside it; two (the meeting of two runs) lengthen the run
    ending at the gap's start and the one starting at its end; each pair in `runs_long` and `runs_short` (`a`'s and
    `b`'s lists the other way round when `swapped`). True when the gap is run letters alone."""
    var cut = at
    while cut < at + count and longer[cut] == longer[at]:
        cut += 1
    for k in range(cut, at + count):
        if longer[k] != longer[cut]:
            return False
    if cut == at + count:
        # One letter: the shorter's run beside the gap, after it or before it.
        var x = longer[at]
        var side = at_short if at_short < len(shorter) and shorter[at_short] == x else (
            at_short - 1 if at_short > 0 and shorter[at_short - 1] == x else -1
        )
        if side < 0:
            return False
        return _pair(
            longer,
            shorter,
            _run_start(longer, at),
            _run_start(shorter, side),
            count,
            noise_run,
            runs_long,
            runs_short,
            swapped,
        )
    # Two letters: the first block ends the run before the gap, the second starts the run after it.
    if at_short == 0 or at_short >= len(shorter):
        return False
    if shorter[at_short - 1] != longer[at] or shorter[at_short] != longer[cut]:
        return False
    var first = _pair(
        longer,
        shorter,
        _run_start(longer, at),
        _run_start(shorter, at_short - 1),
        cut - at,
        noise_run,
        runs_long,
        runs_short,
        swapped,
    )
    var second = _pair(longer, shorter, cut, at_short, at + count - cut, noise_run, runs_long, runs_short, swapped)
    return first and second


def _pair(
    longer: List[UInt8],
    shorter: List[UInt8],
    s_long: Int,
    s_short: Int,
    count: Int,
    noise_run: Int,
    mut runs_long: List[Int],
    mut runs_short: List[Int],
    swapped: Bool,
) -> Bool:
    if swapped:
        return _pair_runs(shorter, longer, s_short, s_long, count, noise_run, runs_short, runs_long)
    return _pair_runs(longer, shorter, s_long, s_short, count, noise_run, runs_long, runs_short)


def _pair_runs(
    a: List[UInt8],
    b: List[UInt8],
    sa: Int,
    sb: Int,
    count: Int,
    noise_run: Int,
    mut runs_a: List[Int],
    mut runs_b: List[Int],
) -> Bool:
    """`_run_pairs`'s test of one block of `count` letters between `a`'s run at `sa` and `b`'s at `sb`: a pair only
    where the two runs' lengths differ by the gap alone (a substitution in or beside a run pairs runs of different
    extent: C TTTT G TTTT C against C T10 C would pair a T4 with the T10) and the longer is at least `noise_run`, the
    shorter `noise_run - 1`, whichever of the two holds it (`a` against `b` as `b` against `a`); each run once."""
    var la = _run_at(a, sa)
    var lb = _run_at(b, sb)
    if abs(la - lb) != count or max(la, lb) < noise_run or min(la, lb) < noise_run - 1:
        return False
    if len(runs_a) > 0 and (runs_a[len(runs_a) - 1] == sa or runs_b[len(runs_b) - 1] == sb):
        return True
    runs_a.append(sa)
    runs_b.append(sb)
    return True


def _one_unit_apart(a: List[UInt8], b: List[UInt8]) -> Bool:
    """Whether `a` and `b` differ only by one unit, 2 to 6 letters not all one, of a tandem repeat: the longer one's
    extra letters, where the two first differ, repeat the unit beside them."""
    var longer = a.copy() if len(a) > len(b) else b.copy()
    var shorter = b.copy() if len(a) > len(b) else a.copy()
    var p = len(longer) - len(shorter)
    if p < 2 or p > 6:
        return False
    var i = 0
    while i < len(shorter) and longer[i] == shorter[i]:
        i += 1
    for k in range(i, len(shorter)):
        if longer[k + p] != shorter[k]:
            return False
    var one_letter = True
    for k in range(i + 1, i + p):
        one_letter = one_letter and longer[k] == longer[i]
    if one_letter:
        return False
    var after = i + 2 * p <= len(longer)
    var before = i >= p
    for k in range(p):
        after = after and longer[i + p + k] == longer[i + k]
        before = before and longer[i - p + k] == longer[i + k]
    return after or before


def _run_start(seq: List[UInt8], i: Int) -> Int:
    var b = i
    while b > 0 and seq[b - 1] == seq[i]:
        b -= 1
    return b


def _length_in(
    cigar: String, consensus_start: Int, letters: List[UInt8], s: Int, e: Int, letter: UInt8, consensus: List[UInt8]
) -> Int:
    """The length a read reports for `consensus`'s run `s ..< e` of `letter`, from its CIGAR against it (from
    `consensus_start`): the read's own run of the letter where the alignment puts the consensus's (its letters aligned
    within the run or inserted inside it or at its ends, and those of the letter beside them), bounded in the read by
    the consensus's two flank letters (`s - 1` and `e`, each another letter); -1 when the read's run is not so
    bounded (a run at the consensus's end has no flank), or holds another letter (no vote, as shasta2's). The read's
    run, not the letters the alignment puts on the consensus's: a gap placed leftmost can take the run's last letter
    with the start of a repeat unit beside it (consensus C AAAA TATATA G, read C AAAA TATA G: the deletion AT took the
    fourth A and the read's A4 counted 3; at unit costs the unit split into letters dropped either side,
    1=1D3=1D6=). Bounded by the flank letters themselves: a read without the letter between two runs of one letter
    (T4 C T11 read as T15) has its run end at another letter than the flank's, and does not count its whole stretch
    for each (as polish's `_run_lengths` refuses it). A read's vote needs only its own letters there, not one
    alignment of them: requiring each flank aligned to a letter would leave out the reads whose gap a unit beside the
    run placed at the flank."""
    if s == 0 or e >= len(consensus):
        return -1
    var i = consensus_start
    var j = 0
    var count = 0
    var first = -1
    var last = -1
    var at_s = -1
    for c in cigar.as_bytes():
        if c >= UInt8(ord("0")) and c <= UInt8(ord("9")):
            count = count * 10 + Int(c - UInt8(ord("0")))
            continue
        for _ in range(count):
            if i == s and at_s < 0:
                at_s = j
            if c == UInt8(ord("I")):
                if i >= s and i <= e and letters[j] == letter:
                    first = j if first < 0 else first
                    last = j
                elif i > s and i < e:
                    return -1
                j += 1
            elif c == UInt8(ord("D")):
                i += 1
            else:
                if i >= s and i < e:
                    if letters[j] != letter:
                        return -1
                    first = j if first < 0 else first
                    last = j
                i += 1
                j += 1
        count = 0
    if first < 0:
        # Every letter of the run deleted: none in the read where the consensus has it, if its flanks say so.
        if at_s <= 0 or at_s >= len(letters):
            return -1
        if letters[at_s - 1] == consensus[s - 1] and letters[at_s] == consensus[e]:
            return 0
        return -1
    for q in range(first, last + 1):
        if letters[q] != letter:
            return -1
    while first > 0 and letters[first - 1] == letter:
        first -= 1
    while last + 1 < len(letters) and letters[last + 1] == letter:
        last += 1
    if first == 0 or last + 1 >= len(letters):
        return -1
    if letters[first - 1] != consensus[s - 1] or letters[last + 1] != consensus[e]:
        return -1
    return last - first + 1


def _flanks(consensus: List[UInt8], s: Int, e: Int) -> Tuple[Int, Int]:
    """The letters either side of `consensus`'s run `s ..< e` as the model codes them (`FLANK_UNKNOWN` at an end)."""
    var left = Int(encode_base(consensus[s - 1])) if s > 0 else FLANK_UNKNOWN
    var right = Int(encode_base(consensus[e])) if e < len(consensus) else FLANK_UNKNOWN
    return (left if left <= 3 else FLANK_UNKNOWN, right if right <= 3 else FLANK_UNKNOWN)


def _two_lengths(
    mut aligner: Aligner,
    made: List[List[UInt8]],
    starts: List[Int],
    groups: List[List[Int]],
    seqs: List[List[UInt8]],
    anchors: List[Int],
    strands: List[Int],
    model: HomopolymerModel,
    mut shorter: Int,
    mut longer: Int,
    mut single: Int,
    mut said: List[List[Int]],
) raises -> Float64:
    """How much better two true lengths of one run explain the reads' lengths there than one, under `model`: the
    log likelihood ratio in dB, each read under two from either length with probability a half. The run starts at
    `starts[k]` in `made[k]` (k 0 and 1, the two clusters' consensus sequences), and each cluster's reads
    (`groups[k]`) are measured on their own consensus, with its flanks: the two may differ beside the run. The two
    lengths in `shorter` and `longer`, the one in `single` (-1 where the model cannot say), searched from 4 under the
    shorter of the two clusters' lengths to 4 over the longer; each cluster's reads' lengths in `said[k]` (-1 none).
    All the reads together, not each cluster's: clusters split by the run's length hold the reads biased apart.
    Where the model does not reach that far (a run longer than its rows), it cannot say: the two lengths the
    clusters have, apart (`Float64.MAX`) when they differ by more than a letter, as without a model, else one
    (`-Float64.MAX`)."""
    single = -1
    var base = Int(encode_base(made[0][starts[0]]))
    var run_length = List[Int]()
    var lengths = List[Int]()
    var on = List[Int]()
    var lefts = List[Int]()
    var rights = List[Int]()
    said.clear()
    for k in range(2):
        var s = starts[k]
        var e = s
        while e < len(made[k]) and made[k][e] == made[k][s]:
            e += 1
        run_length.append(e - s)
        var f = _flanks(made[k], s, e)
        var mine = List[Int]()
        for r in groups[k]:
            var found = aligner.align(Span(made[k]), Span(seqs[r]), Costs.edit(), _mode(anchors[r]), ties=Ties.LEFT)
            var m = _length_in(found.cigar, found.reference_start, seqs[r], s, e, made[k][s], made[k])
            mine.append(m)
            if m >= 0:
                lengths.append(m)
                on.append(strands[r])
                lefts.append(f[0])
                rights.append(f[1])
        said.append(mine^)
    if base > 3 or len(lengths) < 2:
        return -Float64.MAX
    var low = max(2, min(run_length[0], run_length[1]) - 4)
    var high = max(run_length[0], run_length[1]) + 4
    if not model.covers(base, high):
        shorter = min(run_length[0], run_length[1])
        longer = max(run_length[0], run_length[1])
        return Float64.MAX if longer - shorter > 1 else -Float64.MAX
    var one = -Float64.MAX
    for n in range(low, high + 1):
        var total = 0.0
        for k in range(len(lengths)):
            total += model.log_probability(base, lefts[k], rights[k], n, lengths[k], on[k])
        # Ties to the longer: reads lose letters in runs.
        if total >= one:
            one = total
            single = n
    var two = -Float64.MAX
    shorter = run_length[0]
    longer = run_length[0] + 1
    for n1 in range(low, high + 1):
        for n2 in range(n1 + 1, high + 1):
            var total = 0.0
            for k in range(len(lengths)):
                var x = model.log_probability(base, lefts[k], rights[k], n1, lengths[k], on[k])
                var y = model.log_probability(base, lefts[k], rights[k], n2, lengths[k], on[k])
                var top = max(x, y)
                total += top + 10.0 * log10(0.5 * 10.0 ** ((x - top) / 10.0) + 0.5 * 10.0 ** ((y - top) / 10.0))
            if total > two:
                two = total
                shorter = n1
                longer = n2
    return two - one


def _run_at(seq: List[UInt8], i: Int) -> Int:
    """The length of the run of one letter holding position `i` of `seq`."""
    var b = i
    var e = i
    while b > 0 and seq[b - 1] == seq[i]:
        b -= 1
    while e + 1 < len(seq) and seq[e + 1] == seq[i]:
        e += 1
    return e - b + 1


def _first_split(
    mut scratch: ConsensusScratch,
    seqs: List[List[UInt8]],
    spanning: List[Int],
    max_clusters: Int,
    min_freq: Float64,
    options: ConsensusOptions,
) raises -> List[List[Int]]:
    """The spanning reads split as `cluster_rows` splits their MSA's rows (one cluster when it finds no split)."""
    var names = List[String](capacity=len(spanning))
    var mine = List[List[UInt8]](capacity=len(spanning))
    for r in spanning:
        names.append(String(t"r{r}"))
        mine.append(seqs[r].copy())
    var made = poa_with_scratch(
        names, mine, List[Int](length=len(spanning), fill=BOTH), options.poa, Graph(), scratch.poa, scratch.aligner
    )
    var by_path = made.graph.msa()
    var rows = List[String](capacity=len(spanning))
    for k in range(len(spanning)):
        rows.append(by_path[made.shown[k]])
    var out = List[List[Int]]()
    for part in cluster_rows(rows, max_clusters, min_freq):
        var reads = List[Int](capacity=len(part))
        for k in part:
            reads.append(spanning[k])
        out.append(reads^)
    return out^


def clustered_consensus(
    var seqs: List[List[UInt8]],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    max_clusters: Int = 2,
    min_freq: Float64 = 0.25,
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> Clustered:
    """At most `max_clusters` consensus sequences of `seqs` (each anchored as `anchors` says: `BOTH`, `LEFT`,
    `RIGHT` or `INSIDE`; empty: all `BOTH`; on `strands`' strand for `homopolymers`), each from the reads nearest
    it and made as `consensus` makes one, a cluster kept apart only with at least `min_freq` of the reads that make
    consensus sequences; and every read's alignment to its cluster's consensus (at `options.poa.scoring`'s costs,
    gaps leftmost). With one cluster, or none found, the consensus of every read."""
    var scratch = ConsensusScratch()
    return clustered_consensus(scratch, seqs^, anchors, strands, max_clusters, min_freq, options, homopolymers)


def clustered_consensus(
    mut scratch: ConsensusScratch,
    var seqs: List[List[UInt8]],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    max_clusters: Int = 2,
    min_freq: Float64 = 0.25,
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> Clustered:
    """As `clustered_consensus`, reusing `scratch` across calls."""
    var n = len(seqs)
    if max_clusters < 1:
        raise Error(t"max_clusters {max_clusters}: must be 1 or more")
    if not (min_freq >= 0.0 and min_freq <= 1.0):
        raise Error(t"min_freq {min_freq}: must be within 0 and 1")
    if len(anchors) != 0 and len(anchors) != n:
        raise Error(t"{len(anchors)} anchors for {n} sequences")
    if len(strands) != 0 and len(strands) != n:
        raise Error(t"{len(strands)} strands for {n} sequences")
    var anchor_of = anchors.copy() if len(anchors) > 0 else List[Int](length=n, fill=BOTH)
    var strand_of = strands.copy() if len(strands) > 0 else List[Int](length=n, fill=STRAND_UNKNOWN)
    var out = Clustered()
    if n == 0:
        out.consensus.append(List[UInt8]())
        out.members.append(List[Int]())
        return out^
    # The reads that make consensus sequences; with every read inside, all of them as if spanning.
    var builders = List[Int]()
    for r in range(n):
        if anchor_of[r] != INSIDE:
            builders.append(r)
    var build_anchor = anchor_of.copy()
    if len(builders) == 0:
        for r in range(n):
            builders.append(r)
            build_anchor[r] = BOTH
    # `consensus` checks the anchors and strands as it is first called.
    var groups = List[List[Int]]()
    if max_clusters >= 2:
        var spanning = List[Int]()
        for r in builders:
            if build_anchor[r] == BOTH:
                spanning.append(r)
        if len(spanning) >= 2:
            groups = _first_split(scratch, seqs, spanning, max_clusters, min_freq, options)
    if len(groups) <= 1:
        groups = [builders.copy()]
    var costs = alignment_costs(options.poa.scoring)
    # A read moves to the consensus it differs from least, a substitution weighing more than a letter slipped in a run
    # (reads slip in runs far more often than they misread a letter), as the first split weighs them. At the
    # scoring's costs a gap letter cost more than a mismatch, so reads moved by their run lengths, not by a SNP's
    # letter: HG002 chr20:56,892,630, G in a poly-T on 28 of 62 reads: the G cluster came to hold one T-run length's
    # reads instead, its consensus lost the G, and the two merged as apart by a run alone (written G on both).
    var moving = MOVING_COSTS
    ref aligner = scratch.aligner
    var made = List[List[UInt8]]()
    for g in groups:
        made.append(_consensus_of(scratch, seqs, g, build_anchor, strand_of, options, homopolymers))
    if len(groups) >= 2:
        var floor = max(1, Int(min_freq * Float64(len(builders)) + 0.999999))
        var label = List[Int](length=n, fill=-1)
        for k in range(len(groups)):
            for r in groups[k]:
                label[r] = k
        for _ in range(MAX_REFINEMENTS):
            var moved = False
            var next = List[List[Int]](length=len(made), fill=List[Int]())
            for r in builders:
                var best = label[r] if label[r] >= 0 else 0
                var best_cost = _aligned(aligner, made[best], seqs[r], build_anchor[r], moving).cost
                for k in range(len(made)):
                    if k == best:
                        continue
                    var cost = _aligned(aligner, made[k], seqs[r], build_anchor[r], moving).cost
                    if cost < best_cost:
                        best = k
                        best_cost = cost
                if best != label[r]:
                    moved = True
                    label[r] = best
                next[best].append(r)
            var kept = List[List[Int]]()
            for g in next:
                if len(g) >= floor:
                    kept.append(g.copy())
            if len(kept) < len(next):
                # A cluster under the floor: its reads are another's, so one cluster fewer.
                groups = kept^
                if len(groups) <= 1:
                    break
                # Every label again: the dropped cluster's reads kept its number, past the clusters left (or another's).
                for r in builders:
                    label[r] = -1
                for k in range(len(groups)):
                    for r in groups[k]:
                        label[r] = k
                made.clear()
                for g in groups:
                    made.append(_consensus_of(scratch, seqs, g, build_anchor, strand_of, options, homopolymers))
                continue
            if not moved:
                break
            groups = next^
            made.clear()
            for g in groups:
                made.append(_consensus_of(scratch, seqs, g, build_anchor, strand_of, options, homopolymers))
        # A read of a cluster dropped in the last round joins the consensus it differs from least.
        if len(groups) >= 2:
            var placed = List[Bool](length=n, fill=False)
            for g in groups:
                for r in g:
                    placed[r] = True
            for r in builders:
                if placed[r]:
                    continue
                var best = 0
                var best_cost = _aligned(aligner, made[0], seqs[r], build_anchor[r], moving).cost
                for k in range(1, len(made)):
                    var cost = _aligned(aligner, made[k], seqs[r], build_anchor[r], moving).cost
                    if cost < best_cost:
                        best = k
                        best_cost = cost
                groups[best].append(r)
        # Clusters spelling the same consensus are one; so are clusters apart only by homopolymer lengths, unless a
        # homopolymer model finds two lengths at one of those runs (`options.two_lengths_db`).
        # Each pair that is one joins (`root`: union of clusters); with two clusters, as one flag.
        var root = List[Int]()
        for k in range(len(groups)):
            root.append(k)
        var runs = List[Int]()
        var single = True
        # Over the clusters still standing: a refinement that leaves one stops with `made` as it was (the two
        # consensus sequences of the round before), and one cluster has nothing to merge.
        for a in range(len(groups)):
            for b in range(a + 1, len(groups)):
                if made[a] == made[b]:
                    _join(root, a, b)
                elif _one_unit_apart(made[a], made[b]) and Float64(
                    min(len(groups[a]), len(groups[b]))
                ) < MIN_UNIT_SPLIT_FREQ * Float64(len(builders)):
                    _join(root, a, b)
                else:
                    if not homopolymers.loaded():
                        # Without a model: apart unless every difference is a single letter of a run.
                        if _apart_by_runs_alone(aligner, made[a], made[b], options.noise_run, runs, single) and single:
                            _join(root, a, b)
                        continue
                    # With a model, each run the two differ by, whatever else they differ by: two lengths where all
                    # the reads fit two better than one (`_two_lengths`), the cluster with the shorter run the
                    # shorter, else the one length in both. Clusters apart by runs alone stay apart only where a
                    # run has two. One run's two lengths no longer make every other run two (each was fitted to
                    # the pair it found), nor does another difference spare a run the test.
                    var runs_a = List[Int]()
                    var runs_b = List[Int]()
                    var only_runs = _run_pairs(aligner, made[a], made[b], options.noise_run, runs_a, runs_b)
                    var set_a = List[Int]()
                    var set_b = List[Int]()
                    var two = False
                    var pair = List[List[UInt8]]()
                    pair.append(made[a].copy())
                    pair.append(made[b].copy())
                    var pair_groups = List[List[Int]]()
                    pair_groups.append(groups[a].copy())
                    pair_groups.append(groups[b].copy())
                    for k in range(len(runs_a)):
                        var len_a = _run_at(made[a], runs_a[k])
                        var len_b = _run_at(made[b], runs_b[k])
                        var n1 = 0
                        var n2 = 0
                        var one = -1
                        var said = List[List[Int]]()
                        var starts = List[Int]()
                        starts.append(runs_a[k])
                        starts.append(runs_b[k])
                        var gain = _two_lengths(
                            aligner,
                            pair,
                            starts,
                            pair_groups,
                            seqs,
                            build_anchor,
                            strand_of,
                            homopolymers,
                            n1,
                            n2,
                            one,
                            said,
                        )
                        if gain > options.two_lengths_db:
                            two = True
                            set_a.append(n1 if len_a < len_b else n2)
                            set_b.append(n2 if len_a < len_b else n1)
                        elif one > 0:
                            set_a.append(one)
                            set_b.append(one)
                        else:
                            set_a.append(len_a)
                            set_b.append(len_b)
                        # A cluster whose reads report the letter once there has no run, as in `phased_consensus`:
                        # the search starts at 2.
                        if len(said) == 2 and reads_say_one(said[0]):
                            set_a[k] = 1
                        if len(said) == 2 and reads_say_one(said[1]):
                            set_b[k] = 1
                    made[a] = _set_runs(made[a].copy(), runs_a, set_a)
                    made[b] = _set_runs(made[b].copy(), runs_b, set_b)
                    if (only_runs and not two) or made[a] == made[b]:
                        _join(root, a, b)
        # The clusters joined made one again from their reads; with every cluster one, the consensus of every read.
        # One flag for all pairs made three clusters one wherever two of them were (a third haplotype's lost).
        var merged = List[List[Int]]()
        var made_merged = List[List[UInt8]]()
        for k in range(len(groups)):
            if _root(root, k) != k:
                continue
            var members = List[Int]()
            var parts = 0
            for j in range(len(groups)):
                if _root(root, j) == k:
                    members.extend(groups[j].copy())
                    parts += 1
            made_merged.append(
                made[k].copy() if parts
                == 1 else _consensus_of(scratch, seqs, members, build_anchor, strand_of, options, homopolymers)
            )
            merged.append(members^)
        if len(groups) <= 1 or len(merged) <= 1:
            groups = [builders.copy()]
            made = [_consensus_of(scratch, seqs, builders, build_anchor, strand_of, options, homopolymers)]
        else:
            groups = merged^
            made = made_merged^
    var cluster_of = List[Int](length=n, fill=-1)
    for k in range(len(groups)):
        for r in groups[k]:
            cluster_of[r] = k
    return _placed(aligner, made^, seqs, cluster_of, anchor_of, costs)


def _root(root: List[Int], k: Int) -> Int:
    var r = k
    while root[r] != r:
        r = root[r]
    return r


def _join(mut root: List[Int], a: Int, b: Int):
    """Clusters `a` and `b` one (`root`: each cluster's parent, a root its own)."""
    var ra = _root(root, a)
    var rb = _root(root, b)
    if ra != rb:
        root[max(ra, rb)] = min(ra, rb)


def _placed(
    mut aligner: Aligner,
    var made: List[List[UInt8]],
    seqs: List[List[UInt8]],
    cluster_of: List[Int],
    anchors: List[Int],
    costs: Costs,
) raises -> Clustered:
    """The clusters' consensus sequences `made` with every read placed: in its cluster (`cluster_of`), or with none
    (-1) in the one it aligns to most cheaply."""
    var out = Clustered()
    out.consensus = made^
    out.members = List[List[Int]](length=len(out.consensus), fill=List[Int]())
    for r in range(len(seqs)):
        var best = cluster_of[r]
        var found: Alignment
        if best >= 0:
            found = _aligned(aligner, out.consensus[best], seqs[r], anchors[r], costs)
        else:
            best = 0
            found = _aligned(aligner, out.consensus[0], seqs[r], anchors[r], costs)
            for k in range(1, len(out.consensus)):
                var other = _aligned(aligner, out.consensus[k], seqs[r], anchors[r], costs)
                if other.cost < found.cost:
                    best = k
                    found = other^
        out.members[best].append(r)
        out.placements.append(
            Placement(
                best,
                found.cigar,
                found.cost,
                found.reference_start,
                found.reference_end,
                found.query_start,
                found.query_end,
            )
        )
    return out^


def _set_runs(a: List[UInt8], runs: List[Int], lengths: List[Int]) -> List[UInt8]:
    """`a` with its runs starting at `runs` (ascending) `lengths` long."""
    var out = List[UInt8](capacity=len(a) + 8)
    var at = 0
    for k in range(len(runs)):
        var s = runs[k]
        var e = s
        while e < len(a) and a[e] == a[s]:
            e += 1
        for i in range(at, s):
            out.append(a[i])
        for _ in range(lengths[k]):
            out.append(a[s])
        at = e
    for i in range(at, len(a)):
        out.append(a[i])
    return out^


def _read_lengths(
    mut aligner: Aligner,
    consensus: List[UInt8],
    s: Int,
    e: Int,
    seqs: List[List[UInt8]],
    reads: List[Int],
    anchors: List[Int],
    mut lengths: List[Int],
) raises:
    """The length each of `reads` reports for `consensus`'s run `s ..< e` (-1 for none, see `_length_in`)."""
    lengths.clear()
    for r in reads:
        var found = aligner.align(Span(consensus), Span(seqs[r]), Costs.edit(), _mode(anchors[r]), ties=Ties.LEFT)
        lengths.append(_length_in(found.cigar, found.reference_start, seqs[r], s, e, consensus[s], consensus))


def _best_length(
    model: HomopolymerModel,
    base: Int,
    left: Int,
    right: Int,
    lengths: List[Int],
    strands: List[Int],
    reads: List[Int],
    low: Int,
    high: Int,
    mut score: Float64,
) -> Int:
    """The run length in `low ..= high` the reads' lengths make most likely under `model` (ties to the longer: reads
    lose letters in runs), its log likelihood in dB in `score`; reads with no length (-1) do not vote."""
    var best = -1
    score = -Float64.MAX
    for n in range(low, high + 1):
        var total = 0.0
        var any = False
        for k in range(len(reads)):
            if lengths[k] < 0:
                continue
            any = True
            total += model.log_probability(base, left, right, n, lengths[k], strands[reads[k]])
        if any and total >= score:
            score = total
            best = n
    return best


def _best_shared_length(
    model: HomopolymerModel,
    base: Int,
    flanks0: Tuple[Int, Int],
    flanks1: Tuple[Int, Int],
    lengths0: List[Int],
    lengths1: List[Int],
    strands: List[Int],
    reads0: List[Int],
    reads1: List[Int],
    low: Int,
    high: Int,
    mut score: Float64,
) -> Int:
    """`_best_length` over two haplotypes' reads together, each read scored with its own haplotype's flanks: the one
    length both share most likely, its log likelihood in dB in `score` (ties to the longer); -1 with no read."""
    var best = -1
    score = -Float64.MAX
    for n in range(low, high + 1):
        var total = 0.0
        var any = False
        for k in range(len(reads0)):
            if lengths0[k] >= 0:
                any = True
                total += model.log_probability(base, flanks0[0], flanks0[1], n, lengths0[k], strands[reads0[k]])
        for k in range(len(reads1)):
            if lengths1[k] >= 0:
                any = True
                total += model.log_probability(base, flanks1[0], flanks1[1], n, lengths1[k], strands[reads1[k]])
        if any and total >= score:
            score = total
            best = n
    return best


comptime TO_BOTH = 2
"""`_sides`: a read for both haplotypes, covering no place where the two differ."""

comptime MIN_CARRIERS = 1
"""`_sides`: a place where the two haplotypes' consensus sequences differ decides reads only where each spelling is
carried by at least this many reads of the window, labelled or not, and by `MIN_CARRIER_SHARE` of those spelling it
either way. The share does the work: a spelling few of many reads carry is a read's error made the consensus's, not a
second haplotype. HG002's HLA-DPA1 thinned to 15x had haplotype 2's left half built from 1 or 2 labelled reads, whose
36 errors read as 36 places apart from haplotype 1, and every unlabelled read, carrying the truth there, went to
haplotype 1 alone. Counting instead the labelled reads behind each consensus (3 of each) set aside the true
differences low depth leaves under 3 labelled reads, and the unlabelled reads there, given to both haplotypes, made
heterozygotes homozygous (CYP4A11 among them, unflagged); a floor of 2 carriers did the same at 10x and below (NAT2,
MTHFR, HTR2A, GRIK1, XPC), a second allele there being one read in four."""

comptime MIN_CARRIER_SHARE = 0.25
"""`_sides`: the least share of the reads spelling a place either way that each spelling must have to decide reads
there, a second haplotype being about half. On HG002 thinned to 30x down to 5x (two read sets) this kept every
HLA-DPA1 consensus exact at 15x and 10x and left no wrong call unflagged in panStar; it still set aside a second allele
one read in five carried at 10x and below (NAT2, MTHFR thinned), as the reads alone cannot tell it from an error
there: at 0.15 those were kept, but so were DPA1's errors where one read in six carried them (33 edits at 15x)."""

comptime PHASED_CARRIERS = 3
"""`_sides`: a place also decides reads where pgphase's labels back it: each haplotype's spelling carried by at least
this many of its own labelled reads, and by `PHASED_AGREEMENT` of its labelled reads spelling it either way. The
share of all reads cannot tell a minority haplotype from an error (NA18941's CYP1A1: 8 reads of one haplotype, 3 of
the other, 4 of another phase set carrying the first's letter; at 3 of 15 the place was set aside, the 4 went to
both and outvoted the second's own 3), but an error does not follow the labels: one read's error made a thin
consensus's (HG002's HLA-DPA1 at 15x, a haplotype built from 1 or 2 reads) has under 3 carriers of its own."""

comptime PHASED_AGREEMENT = 0.75

comptime SITE_FLANK = 8
"""Letters of the longer haplotype either side of a place where the two differ, compared with a read's letters there
along with it: two alignments may put the same gap a letter or two apart, which compared alone would read as the
other haplotype's spelling; with the flank both are one string. Places within twice this of each other are one."""


@fieldwise_init
struct _Site(Copyable, Movable):
    """A place where two haplotypes differ, on the longer's letters `start ..< end` (`SITE_FLANK` either side
    included) and the shorter's `short_start ..< short_end`: the longer spells it `longer`, the shorter `shorter`."""

    var start: Int
    var end: Int
    var short_start: Int
    var short_end: Int
    var longer: List[UInt8]
    var shorter: List[UInt8]


def _ops(cigar: String) -> List[Tuple[Int, UInt8]]:
    """karan's run-length CIGAR as (length, letter) pairs."""
    var out = List[Tuple[Int, UInt8]]()
    var number = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            number = number * 10 + Int(byte - UInt8(ord("0")))
        else:
            out.append((number, byte))
            number = 0
    return out^


def _query_at(alignment: Alignment) -> List[Int]:
    """For each reference letter `reference_start ..= reference_end` of `alignment`, the query letter it stands
    before (a gap's letters count before the reference letter after them): the query's letters over reference
    `s ..< e` are `query_at[s - reference_start] ..< query_at[e - reference_start]`."""
    var out = List[Int](capacity=alignment.reference_end - alignment.reference_start + 1)
    var q = alignment.query_start
    var ops = _ops(alignment.cigar)
    for ref op in ops:
        var length = op[0]
        var letter = op[1]
        if letter == UInt8(ord("I")):
            q += length
            continue
        for _ in range(length):
            out.append(q)
            if letter != UInt8(ord("D")):
                q += 1
    out.append(q)
    return out^


def _differences(mut aligner: Aligner, longer: List[UInt8], shorter: List[UInt8], costs: Costs) raises -> List[_Site]:
    """Where `shorter`, placed inside `longer` (`Mode.INFIX`: a consensus cut short where its reads stop is not
    thereby another haplotype), spells it otherwise, each place widened by `SITE_FLANK` and those that then touch
    joined."""
    var out = List[_Site]()
    if len(longer) == 0 or len(shorter) == 0 or longer == shorter:
        return out^
    var placed = aligner.align(Span(longer), Span(shorter), costs, Mode.INFIX, ties=Ties.LEFT)
    var at = _query_at(placed)
    var spans = List[Tuple[Int, Int]]()
    var r = placed.reference_start
    var ops = _ops(placed.cigar)
    for ref op in ops:
        var length = op[0]
        var letter = op[1]
        if letter == UInt8(ord("=")) or letter == UInt8(ord("M")):
            r += length
            continue
        var consumed = 0 if letter == UInt8(ord("I")) else length
        var start = max(placed.reference_start, r - SITE_FLANK)
        var end = min(placed.reference_end, r + consumed + SITE_FLANK)
        if len(spans) > 0 and start <= spans[len(spans) - 1][1]:
            spans[len(spans) - 1] = (spans[len(spans) - 1][0], max(end, spans[len(spans) - 1][1]))
        else:
            spans.append((start, end))
        r += consumed
    for ref span in spans:
        var long_letters = List[UInt8]()
        for i in range(span[0], span[1]):
            long_letters.append(longer[i])
        var short_start = at[span[0] - placed.reference_start]
        var short_end = at[span[1] - placed.reference_start]
        var short_letters = List[UInt8]()
        for i in range(short_start, short_end):
            short_letters.append(shorter[i])
        out.append(_Site(span[0], span[1], short_start, short_end, long_letters^, short_letters^))
    return out^


def _haplotype_at(alignment: Alignment) -> List[Int]:
    """`_query_at` the other way round: for a haplotype placed inside a read (`alignment`'s query the haplotype, its
    reference the read; `Mode.INFIX`, so every haplotype letter is placed), for each haplotype letter `0 ..= len` the
    read letter it stands before, a gap's letters counting before the haplotype letter after them."""
    var out = List[Int](capacity=alignment.query_end - alignment.query_start + 1)
    var r = alignment.reference_start
    var ops = _ops(alignment.cigar)
    for ref op in ops:
        var length = op[0]
        var letter = op[1]
        if letter == UInt8(ord("D")):
            r += length
            continue
        for _ in range(length):
            out.append(r)
            if letter != UInt8(ord("I")):
                r += 1
    out.append(r)
    return out^


def _spelling(seq: List[UInt8], start: Int, end: Int, at: List[Int], lo: Int, hi: Int, site: _Site) -> Int:
    """Which spelling of `site` a read carries, placed over a haplotype's letters `start ..< end` (`at`: the read
    letter each haplotype letter from `start` stands before), the site lying on that haplotype's `lo ..< hi`: 0 the
    longer's spelling, 1 the shorter's, -1 neither, or a place the read does not span wholly short of its
    placement's ends, where a read running past the haplotype's end piles its overhang."""
    if lo <= start or hi >= end:
        return -1
    var letters = List[UInt8]()
    for q in range(at[lo - start], at[hi - start]):
        letters.append(seq[q])
    if letters == site.longer:
        return 0
    if letters == site.shorter:
        return 1
    return -1


def _sides(
    mut aligner: Aligner,
    made: List[List[UInt8]],
    seqs: List[List[UInt8]],
    builders: List[Int],
    reads: List[Int],
    costs: Costs,
    labels: List[Int] = List[Int](),
) raises -> List[Int]:
    """For each of `reads`, the haplotype of `made` whose letters it carries where the two differ: 0 or 1, `TO_BOTH`
    when it covers no such place (a stretch where the two read alike, or one neither reaches), -1 when it carries as
    many of one's as of the other's. Every read of `builders` (labelled or not) is placed on a haplotype
    (`Mode.INFIX`) and its letters over each place it spans wholly compared with the two spellings there -- a read
    labelled with the shorter's haplotype on the shorter, at each place's span there, every other read on the
    longer; a place
    decides only where each spelling is carried by `MIN_CARRIERS` of them and `MIN_CARRIER_SHARE` of those spelling
    it either way, or where `labels` (a haplotype per read of `seqs`, -1 none) back it (`PHASED_CARRIERS`)."""
    var out = List[Int](length=len(reads), fill=TO_BOTH)
    var long_index = 0 if len(made[0]) >= len(made[1]) else 1
    var short_index = 1 - long_index
    ref longer = made[long_index]
    var sites = _differences(aligner, longer, made[short_index], costs)
    if len(sites) == 0:
        return out^
    # Each read's spelling at each place: row per read of `builders`.
    var spelt = List[List[Int]]()
    var row_of = List[Int](length=len(seqs), fill=-1)
    for b in range(len(builders)):
        ref seq = seqs[builders[b]]
        # A read labelled with the shorter's haplotype is placed on the shorter, the places read at their spans there
        # (`_Site.short_start`): it carries that haplotype's letters, and placed on the longer it would differ from it
        # wherever the two do. NA19007's HLA-DRB1 (*09:01 and *08:03, hundreds of edits apart through the introns):
        # 77 reads placed on the other haplotype at half a second each, 39 of the window's 47 s, for 2 unlabelled
        # reads. Every other read is placed on the longer.
        var h = labels[builders[b]] if len(labels) > 0 else -1
        var on_shorter = h == short_index
        ref haplotype = made[short_index if on_shorter else long_index]
        # A read is placed inside the haplotype, unless it is the longer of the two: then the haplotype inside it.
        # A haplotype built from a few labelled reads can be a fraction of the window (HG02106's HLA-DPB1: 2,716 and
        # 3,732 letters from 7 of 71 reads), and a 15.7 kb read forced inside 3.7 kb is a 12 kb gap, which an
        # alignment's cost grows with: 1.8 s a read, 57 of the window's 59 s. The other way round the two differ by
        # what the read's errors and the haplotypes' differences make, and every haplotype letter has its read
        # letter.
        var start: Int
        var end: Int
        var at: List[Int]
        if len(seq) > len(haplotype):
            var inside = aligner.align(Span(seq), Span(haplotype), costs, Mode.INFIX, ties=Ties.LEFT)
            at = _haplotype_at(inside)
            start = 0
            end = len(haplotype)
        else:
            var placed = aligner.align(Span(haplotype), Span(seq), costs, Mode.INFIX, ties=Ties.LEFT)
            at = _query_at(placed)
            start = placed.reference_start
            end = placed.reference_end
        var row = List[Int](length=len(sites), fill=-1)
        for k in range(len(sites)):
            ref site = sites[k]
            if on_shorter:
                row[k] = _spelling(seq, start, end, at, site.short_start, site.short_end, site)
            else:
                row[k] = _spelling(seq, start, end, at, site.start, site.end, site)
        spelt.append(row^)
        row_of[builders[b]] = b
    var counted = List[Bool](length=len(sites), fill=False)
    for k in range(len(sites)):
        var carriers = [0, 0]
        # own[2 * h + c]: reads labelled h spelling the place c (0 the longer haplotype's spelling).
        var own = List[Int](length=4, fill=0)
        for b in range(len(spelt)):
            var c = spelt[b][k]
            if c < 0:
                continue
            carriers[c] += 1
            var h = labels[builders[b]] if len(labels) > 0 else -1
            if h == 0 or h == 1:
                own[2 * h + c] += 1
        var fewer = min(carriers[0], carriers[1])
        counted[k] = fewer >= MIN_CARRIERS and Float64(fewer) >= MIN_CARRIER_SHARE * Float64(carriers[0] + carriers[1])
        if not counted[k]:
            # Labels back it: the longer haplotype's reads spell it 0, the other's 1.
            var mine = own[2 * long_index]
            var theirs = own[2 * short_index + 1]
            var mine_all = own[2 * long_index] + own[2 * long_index + 1]
            var theirs_all = own[2 * short_index] + own[2 * short_index + 1]
            counted[k] = (
                mine >= PHASED_CARRIERS
                and theirs >= PHASED_CARRIERS
                and Float64(mine) >= PHASED_AGREEMENT * Float64(mine_all)
                and Float64(theirs) >= PHASED_AGREEMENT * Float64(theirs_all)
            )
    for i in range(len(reads)):
        var votes = [0, 0]
        ref row = spelt[row_of[reads[i]]]
        for k in range(len(sites)):
            if counted[k] and row[k] >= 0:
                # Spelling 0 is the longer haplotype's.
                votes[long_index if row[k] == 0 else 1 - long_index] += 1
        if votes[0] > votes[1]:
            out[i] = 0
        elif votes[1] > votes[0]:
            out[i] = 1
        elif votes[0] > 0:
            out[i] = -1
    return out^


def phased_consensus(
    var seqs: List[List[UInt8]],
    haplotypes: List[Int],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> Clustered:
    """The two haplotypes' consensus sequences of reads already phased (`haplotypes`: 0 or 1 per read, -1 for a read
    of neither), each made as `consensus` makes one from its reads, and every read's alignment to its haplotype's
    (see `clustered_consensus`). A read of neither builds as well: made first from the labelled reads, the two are
    made again with each such read on the haplotype whose letters it carries where they differ, on both where it
    covers no such place, on neither where it carries as many of each (`_sides`); it is then placed where it fits
    best. Where the two differ by the length of a run of one letter and a homopolymer model is loaded, each such run
    is decided on both haplotypes' reads together, whatever else the two differ by: their own lengths for the two only
    when that explains the reads `options.two_lengths_db` better than one length for both, else that one length in
    both (a run is far more often the same on the two haplotypes than not). Two clusters always, the
    second's consensus the first's when no read is on haplotype 1."""
    var scratch = ConsensusScratch()
    return phased_consensus(scratch, seqs^, haplotypes, anchors, strands, options, homopolymers)


def phased_consensus(
    mut scratch: ConsensusScratch,
    var seqs: List[List[UInt8]],
    haplotypes: List[Int],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> Clustered:
    """As `phased_consensus`, reusing `scratch` across calls."""
    var n = len(seqs)
    if len(haplotypes) != n:
        raise Error(t"{len(haplotypes)} haplotypes for {n} sequences")
    if len(anchors) != 0 and len(anchors) != n:
        raise Error(t"{len(anchors)} anchors for {n} sequences")
    if len(strands) != 0 and len(strands) != n:
        raise Error(t"{len(strands)} strands for {n} sequences")
    var anchor_of = anchors.copy() if len(anchors) > 0 else List[Int](length=n, fill=BOTH)
    var strand_of = strands.copy() if len(strands) > 0 else List[Int](length=n, fill=STRAND_UNKNOWN)
    var groups = List[List[Int]](length=2, fill=List[Int]())
    var cluster_of = List[Int](length=n, fill=-1)
    for r in range(n):
        var h = haplotypes[r]
        if h != 0 and h != 1 and h != -1:
            raise Error(t"read {r}: haplotype {h} is not 0, 1 or -1")
        if h >= 0:
            cluster_of[r] = h
            if anchor_of[r] != INSIDE:
                groups[h].append(r)
    var made = List[List[UInt8]]()
    for k in range(2):
        made.append(
            _consensus_of(scratch, seqs, groups[k], anchor_of, strand_of, options, homopolymers) if len(groups[k])
            > 0 else (List[UInt8]())
        )
    if len(made[1]) == 0 and len(groups[1]) == 0:
        made[1] = made[0].copy()
    if len(made[0]) == 0 and len(groups[0]) == 0:
        made[0] = made[1].copy()
    ref aligner = scratch.aligner
    # Reads of neither haplotype build too: each joins the haplotype whose letters it carries where the two differ,
    # or both where it covers no such place, and the two are made again. Left out of building, they left a haplotype
    # pieced together from its few labelled reads wherever pgphase's phase set ends inside the window: HG002's
    # HLA-DPA1 thinned to 15x, its 14 reads of the window's left end all unlabelled, had haplotype 2's left half from
    # its one spanning read, 36 edits from the truth; at 10x, with no labelled read reaching that end, both
    # haplotypes stopped short, at 3,915 and 8,058 letters of 18,180. With them: 0 and 1 edits, then 0 and 0.
    var unlabelled = List[Int]()
    for r in range(n):
        if haplotypes[r] < 0 and anchor_of[r] != INSIDE:
            unlabelled.append(r)
    if len(unlabelled) > 0 and len(groups[0]) + len(groups[1]) > 0:
        var costs = alignment_costs(options.poa.scoring)
        var builders = groups[0].copy()
        builders.extend(groups[1].copy())
        builders.extend(unlabelled.copy())
        var sides = _sides(aligner, made, seqs, builders, unlabelled, costs, haplotypes)
        var grew = [False, False]
        for i in range(len(unlabelled)):
            for k in range(2):
                if sides[i] == k or sides[i] == TO_BOTH:
                    groups[k].append(unlabelled[i])
                    grew[k] = True
        for k in range(2):
            if grew[k]:
                made[k] = _consensus_of(scratch, seqs, groups[k], anchor_of, strand_of, options, homopolymers)
    # Each run the two differ by is decided on both haplotypes' reads together, whatever else the two differ by.
    var runs0 = List[Int]()
    var runs1 = List[Int]()
    if homopolymers.loaded() and made[0] != made[1]:
        _ = _run_pairs(aligner, made[0], made[1], max(options.noise_run, 2), runs0, runs1)
    if len(runs0) > 0:
        var first = List[Int]()
        var second = List[Int]()
        var lengths0 = List[Int]()
        var lengths1 = List[Int]()
        for k in range(len(runs0)):
            var s = runs0[k]
            var e = s
            while e < len(made[0]) and made[0][e] == made[0][s]:
                e += 1
            var s1 = runs1[k]
            var e1 = s1
            while e1 < len(made[1]) and made[1][e1] == made[1][s1]:
                e1 += 1
            var base = Int(encode_base(made[0][s]))
            # Each haplotype's reads scored with its own consensus's flanks: the two may differ beside the run.
            var flanks0 = _flanks(made[0], s, e)
            var flanks1 = _flanks(made[1], s1, e1)
            # From 4 under the shorter of the two haplotypes' lengths to 4 over the longer: they may be far apart.
            var low = max(2, min(e - s, e1 - s1) - 4)
            var high = max(e - s, e1 - s1) + 4
            if not homopolymers.covers(base, high):
                # A run longer than the model's rows: it cannot say, so each haplotype keeps its own consensus's
                # length (searched within the rows alone, both came out at the longest row: T59 as T50).
                first.append(e - s)
                second.append(e1 - s1)
                continue
            # Each haplotype's reads measured on its own consensus: the two may differ beside the run.
            _read_lengths(aligner, made[0], s, e, seqs, groups[0], anchor_of, lengths0)
            _read_lengths(aligner, made[1], s1, e1, seqs, groups[1], anchor_of, lengths1)
            var shared_score = 0.0
            var score0 = 0.0
            var score1 = 0.0
            var shared = _best_shared_length(
                homopolymers,
                base,
                flanks0,
                flanks1,
                lengths0,
                lengths1,
                strand_of,
                groups[0],
                groups[1],
                low,
                high,
                shared_score,
            )
            var own0 = _best_length(
                homopolymers, base, flanks0[0], flanks0[1], lengths0, strand_of, groups[0], low, high, score0
            )
            var own1 = _best_length(
                homopolymers, base, flanks1[0], flanks1[1], lengths1, strand_of, groups[1], low, high, score1
            )
            if shared < 0 or own0 < 0 or own1 < 0:
                # Too few reads to say: as the two consensus sequences have it.
                first.append(e - s)
                second.append(e1 - s1)
                continue
            if score0 + score1 - shared_score > options.two_lengths_db:
                first.append(own0)
                second.append(own1)
            else:
                first.append(shared)
                second.append(shared)
            # A haplotype whose reads report the letter once there has no run (`reads_say_one`): the search above
            # starts at 2 (chr20:8,254,476, 42 reads T against the other haplotype's 26 TTT: written TT, the truth T).
            if reads_say_one(lengths0):
                first[k] = 1
            if reads_say_one(lengths1):
                second[k] = 1
        var a = made[0].copy()
        var b = made[1].copy()
        made[0] = _set_runs(a, runs0, first)
        made[1] = _set_runs(b, runs1, second)
    return _placed(aligner, made^, seqs, cluster_of, anchor_of, alignment_costs(options.poa.scoring))
