"""Homopolymer runs two consensus sequences differ by, and the reads' lengths there: the runs found by an affine
alignment of the two (`_run_pairs`), each read's length of a run between the letters it aligns to the run's flanks
(`_length_in`), and the length or lengths a homopolymer model finds most likely for them (`_two_lengths`,
`_best_length`, `_best_shared_length`). Read length errors in runs split one haplotype's reads far more often than
a true difference does, so `clustered_consensus` and `phased_consensus` both decide every run they differ by here."""

from std.math import log10

from karan import Aligner, Costs, Mode, Ties

from .alphabet import encode_base
from .consensus import consensus
from .homopolymer import FLANK_UNKNOWN, HomopolymerModel
from .kernel import LEFT
from .placement import _mode


comptime RUN_PAIR_COSTS = Costs(4, 6, 2, -1, 0)
"""`_run_pairs`'s alignment of two consensus sequences: gap-affine, a substitution 4, a gap of k letters 6 + 2k
(WFA2-lib's default), so an indel of several letters stays one gap. It only finds the runs the two differ by (the
records are written later, by the caller's alignment to its reference); minimap2's asm5 (40, 78 + 7k, 162 + 3k) kept
more than karan's 80 MB of wavefronts on long pairs (HG002 chr20's VNTR regions)."""


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
