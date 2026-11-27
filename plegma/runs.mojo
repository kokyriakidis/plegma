"""Homopolymer runs two consensus sequences differ by, and the reads' lengths there: the runs found by an affine
alignment of the two (`_run_pairs`), each read's length of a run between the letters it aligns to the run's flanks
(`_length_in`), and the length or lengths a homopolymer model finds most likely for them (`_two_lengths`,
`_best_length`, `_best_shared_length`).

Read alignments against a cluster's consensus can be computed once (`_align_group_reads`, in parallel when
`threads > 1`) and reused across every differing run (`_two_lengths_from`, `_read_lengths_from`), avoiding
re-aligning all reads for every run."""

from std.atomic import Atomic
from std.math import log10

from karan import Aligner, Alignment, Costs, Mode, Ties
from karan.common import spread, thread_count

from .alphabet import encode_base
from .homopolymer import FLANK_UNKNOWN, HomopolymerModel
from .placement import _mode

comptime RUN_PAIR_COSTS = Costs(4, 6, 2, -1, 0)
"""`_run_pairs`'s alignment of two consensus sequences: gap-affine, a substitution 4, a gap of k letters 6 + 2k."""


def _apart_by_runs_alone(
    mut aligner: Aligner, a: List[UInt8], b: List[UInt8], noise_run: Int, mut runs: List[Int], mut single: Bool
) raises -> Bool:
    """Whether `a` and `b` differ only by letters added to or dropped from runs of one letter at least `noise_run`
    long (never with `noise_run` 0), each such run's start in `a` then in `runs`, and `single` whether each by one
    letter alone."""
    runs.clear()
    single = True
    if noise_run <= 0:
        return False
    var found = aligner.align(Span(a), Span(b), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
    if found.cost == 0:
        return False
    var i = 0
    var j = 0
    var count = 0
    for c in found.cigar.as_bytes():
        if c >= 48 and c <= 57:
            count = count * 10 + Int(c - 48)
            continue
        if c == UInt8(ord("=")):
            i += count
            j += count
        elif c == UInt8(ord("X")):
            return False
        elif c == UInt8(ord("D")):
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
            if _run_at(b, j) < noise_run or _run_at(b, j) - count < noise_run - 1:
                return False
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
    long, whatever else they differ by."""
    runs_a.clear()
    runs_b.clear()
    if noise_run <= 0 or a == b:
        return False
    var cigar: String
    try:
        cigar = aligner.align(Span(a), Span(b), RUN_PAIR_COSTS, Mode.GLOBAL, ties=Ties.LEFT).cigar
    except:
        return False
    var all_runs = True
    var i = 0
    var j = 0
    var count = 0
    for c in cigar.as_bytes():
        if c >= 48 and c <= 57:
            count = count * 10 + Int(c - 48)
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
    var cut = at
    while cut < at + count and longer[cut] == longer[at]:
        cut += 1
    for k in range(cut, at + count):
        if longer[k] != longer[cut]:
            return False
    if cut == at + count:
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
    """The length a read reports for `consensus`'s run `s ..< e` of `letter`, from its CIGAR against it."""
    if s == 0 or e >= len(consensus):
        return -1
    var i = consensus_start
    var j = 0
    var count = 0
    var first = -1
    var last = -1
    var at_s = -1
    for c in cigar.as_bytes():
        if c >= 48 and c <= 57:
            count = count * 10 + Int(c - 48)
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


def _align_group_reads(
    mut aligner: Aligner,
    consensus: List[UInt8],
    seqs: List[List[UInt8]],
    reads: List[Int],
    anchors: List[Int],
    threads: Int = 1,
) raises -> List[Alignment]:
    """Aligns every read in `reads` to `consensus` at edit costs (in parallel when `threads > 1`), so multiple
    differing runs reuse the same alignments."""
    var n = len(reads)
    var out = List[Alignment](capacity=n)
    var workers = thread_count(threads, n)
    if workers > 1 and n >= 2:
        var order = List[UInt64](capacity=n)
        for k in range(n):
            order.append((UInt64(len(seqs[reads[k]])) << 32) | UInt64(k))
        sort(order)
        order.reverse()
        var results = List[Optional[Alignment]](length=n, fill=None)
        var failed = List[Bool](length=n, fill=False)
        var results_p = results.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def align_worker(
            slot: Int,
        ) {mut taken, imm order, imm consensus, imm seqs, imm reads, imm anchors, imm results_p, imm failed_p}:
            var local_aligner = Aligner(warm=False)
            while True:
                var pos = Int(taken.fetch_add(1))
                if pos >= len(order):
                    return
                var k = Int(order[pos] & 0xFFFFFFFF)
                var r = reads[k]
                try:
                    results_p[unsafe_offset=k] = local_aligner.align(
                        Span(consensus), Span(seqs[r]), Costs.edit(), _mode(anchors[r]), ties=Ties.LEFT
                    )
                except:
                    failed_p[unsafe_offset=k] = True

        spread(align_worker, workers, workers)
        for k in range(n):
            if failed[k]:
                var r = reads[k]
                out.append(
                    aligner.align(Span(consensus), Span(seqs[r]), Costs.edit(), _mode(anchors[r]), ties=Ties.LEFT)
                )
            else:
                out.append(results[k].take())
        return out^

    for r in reads:
        out.append(aligner.align(Span(consensus), Span(seqs[r]), Costs.edit(), _mode(anchors[r]), ties=Ties.LEFT))
    return out^


def _two_lengths_from(
    made: List[List[UInt8]],
    starts: List[Int],
    groups: List[List[Int]],
    aligned0: List[Alignment],
    aligned1: List[Alignment],
    seqs: List[List[UInt8]],
    strands: List[Int],
    model: HomopolymerModel,
    mut shorter: Int,
    mut longer: Int,
    mut single: Int,
    mut said: List[List[Int]],
) -> Float64:
    """`_two_lengths` using pre-computed read alignments `aligned0` and `aligned1`."""
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
        var mine = List[Int](capacity=len(groups[k]))
        ref aligned_k = aligned0 if k == 0 else aligned1
        for idx in range(len(groups[k])):
            var r = groups[k][idx]
            ref found = aligned_k[idx]
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
    var n_len = len(lengths)
    var n_span = high - low + 1
    var lp = List[Float64](length=n_span * n_len, fill=0.0)
    var one = -Float64.MAX
    for n in range(low, high + 1):
        var row = (n - low) * n_len
        var total = 0.0
        for k in range(n_len):
            var v = model.log_probability(base, lefts[k], rights[k], n, lengths[k], on[k])
            lp[row + k] = v
            total += v
        if total >= one:
            one = total
            single = n
    var two = -Float64.MAX
    shorter = run_length[0]
    longer = run_length[0] + 1
    for n1 in range(low, high + 1):
        var row1 = (n1 - low) * n_len
        for n2 in range(n1 + 1, high + 1):
            var row2 = (n2 - low) * n_len
            var total = 0.0
            for k in range(n_len):
                var x = lp[row1 + k]
                var y = lp[row2 + k]
                var top = max(x, y)
                total += top + 10.0 * log10(0.5 * 10.0 ** ((x - top) / 10.0) + 0.5 * 10.0 ** ((y - top) / 10.0))
            if total > two:
                two = total
                shorter = n1
                longer = n2
    return two - one


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
    threads: Int = 1,
) raises -> Float64:
    """How much better two true lengths of one run explain the reads' lengths there than one, under `model`."""
    var aligned0 = _align_group_reads(aligner, made[0], seqs, groups[0], anchors, threads)
    var aligned1 = _align_group_reads(aligner, made[1], seqs, groups[1], anchors, threads)
    return _two_lengths_from(
        made, starts, groups, aligned0, aligned1, seqs, strands, model, shorter, longer, single, said
    )


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


def _read_lengths_from(
    aligned: List[Alignment],
    consensus: List[UInt8],
    s: Int,
    e: Int,
    seqs: List[List[UInt8]],
    reads: List[Int],
    mut lengths: List[Int],
):
    """The length each of `reads` reports for `consensus`'s run `s ..< e`, from pre-computed `aligned`."""
    lengths.clear()
    for k in range(len(reads)):
        var r = reads[k]
        ref found = aligned[k]
        lengths.append(_length_in(found.cigar, found.reference_start, seqs[r], s, e, consensus[s], consensus))


def _read_lengths(
    mut aligner: Aligner,
    consensus: List[UInt8],
    s: Int,
    e: Int,
    seqs: List[List[UInt8]],
    reads: List[Int],
    anchors: List[Int],
    mut lengths: List[Int],
    threads: Int = 1,
) raises:
    """The length each of `reads` reports for `consensus`'s run `s ..< e` (-1 for none, see `_length_in`)."""
    var aligned = _align_group_reads(aligner, consensus, seqs, reads, anchors, threads)
    _read_lengths_from(aligned, consensus, s, e, seqs, reads, lengths)


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
    """The run length in `low ..= high` the reads' lengths make most likely under `model` (ties to the longer)."""
    var best = -1
    score = -Float64.MAX
    for n in range(low, high + 1):
        if not model.covers(base, n):
            continue
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
    """`_best_length` over two haplotypes' reads together, each read scored with its own haplotype's flanks."""
    var best = -1
    score = -Float64.MAX
    for n in range(low, high + 1):
        if not model.covers(base, n):
            continue
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
