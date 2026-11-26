"""What polishing aligns: each read aligned to the consensus (`_Realigned`, `_Paths`: karan at unit costs, a read
anchored at one end as it is anchored, `past_the_end`), its forward and backward edit distance matrices in a
band around that alignment (`_Bands`, `_forward_column`, `_backward_column`), reusable scratch buffers
(`PolishScratch`), the consensus cut into segments outside tandem repeats (`RepeatMask`, `_segments`), and the
repeat periods a whole unit is edited by (`_periods`)."""

from std.atomic import Atomic
from std.math import iota

from karan import Aligner, Alignment, Costs, Mode, Ties
from karan.common import spread, thread_count

from .alphabet import base_code, encode_base
from .kernel import LEFT, RIGHT

comptime BAND = 24
"""Read letters either side of a sequence's alignment path that its forward and backward matrices span, and so the
longest repeat unit polishing edits whole (`MAX_PERIOD`)."""
comptime SEGMENT = 64
"""Consensus letters a segment spans at least."""
comptime APART = 64
"""Edits one round takes in a segment are at least this many letters apart."""
comptime REPEAT_K = 12
comptime REPEAT_DISTANCE = 500
"""No segment boundary within a k-mer of a tandem repeat (a 12-mer recurring within 500 letters)."""
comptime INF = Int32(1 << 28)
comptime SLOTS = 9
"""Per consensus position: insert A, C, G or T before it, delete it, substitute A, C, G or T for it."""
comptime OUTLIER = 0.3
"""A sequence more than this share edits (over the longer of the two spans its alignment covers) takes no part in
polishing."""
comptime TIE_RUN = 4
"""Runs of one letter at least this long take the longer length of a tie (`_lengthen_ties`)."""
comptime MAX_PERIOD = BAND
"""Repeat units up to this long are deleted or copied whole."""


@inline(.always)
def _code(letter: UInt8) -> Int:
    """A, C, G and T as 0..3 (upper or lower case); 4 for anything else."""
    return base_code(letter)


struct _Realigned(Movable):
    """The sequences' alignments to one consensus, as `_Paths.fill` found them. A fill against the same consensus
    takes them as they are (an alignment depends on its two sequences alone), and a sequence with the letters and
    anchor of one before it takes that one's (`same_as`): the likelihood rounds mostly start on the consensus the
    edit rounds end on, and the runs on theirs, two thirds of the letters polishing aligned."""

    var consensus: List[UInt8]
    var found: List[Optional[Alignment]]
    """Per sequence, its alignment to `consensus`, once found."""
    var same_as: List[Int]
    """Per sequence, the first member with its letters and anchor (itself if none before it)."""

    def __init__(out self, var same_as: List[Int]):
        self.consensus = List[UInt8]()
        self.found = List[Optional[Alignment]]()
        self.same_as = same_as^

    def on(mut self, consensus: List[UInt8]):
        """The alignments kept to `consensus`: those to another forgotten."""
        if len(self.found) == len(self.same_as) and self.consensus == consensus:
            return
        self.consensus = consensus.copy()
        self.found = List[Optional[Alignment]](length=len(self.same_as), fill=None)


def _align_one(mut aligner: Aligner, consensus: List[UInt8], seq: List[UInt8], anchor: Int) raises -> Alignment:
    """Aligns one read to `consensus` at unit costs according to its `anchor` (`BOTH`, `LEFT` or `RIGHT`)."""
    var mode = Mode.GLOBAL
    if anchor == LEFT:
        mode = Mode.PREFIX
    elif anchor == RIGHT:
        mode = Mode.SUFFIX
    var aligned = aligner.align(Span(consensus), Span(seq), Costs.edit(), mode, ties=Ties.LEFT)
    # Every consensus letter the read lacks costs the other way at least one: tried only where that could still win.
    if (anchor == LEFT or anchor == RIGHT) and len(consensus) - len(seq) < aligned.cost:
        var other = past_the_end(aligner, consensus, seq, anchor == LEFT)
        if other.cost < aligned.cost:
            aligned = other^
    return aligned^


struct _Paths(Movable):
    """Every sequence's alignment to the consensus: at consensus boundary i (before letter i), sequence k's path
    enters at read offset `enter[k * stride + i]` and leaves at `leave[k * stride + i]`, the letters between
    inserted before letter i."""

    var enter: List[Int32]
    var leave: List[Int32]
    var discord: List[Int32]
    """Per consensus boundary i: how many sequences insert before letter i, or delete or substitute it."""
    var cover: List[Int32]
    """Per consensus boundary i: how many sequences' alignments span it (a read anchored at one end spans a part)."""
    var span_start: List[Int]
    var span_end: List[Int]
    """Per sequence: the consensus boundaries its alignment spans, span_start ..= span_end."""
    var spanned: Int
    """The consensus letters the sequences' alignments span, summed."""
    var total: Int
    var stride: Int
    var matched: Int
    var substituted: Int
    var inserted: Int
    var deleted: Int
    """The alignments' letters of each kind, summed over the sequences."""
    var copying: Int
    """Inserted letters the same as a consensus letter beside them."""
    var worker_aligners: List[Aligner]

    def __init__(out self):
        self.enter = List[Int32]()
        self.leave = List[Int32]()
        self.discord = List[Int32]()
        self.cover = List[Int32]()
        self.span_start = List[Int]()
        self.span_end = List[Int]()
        self.spanned = 0
        self.total = 0
        self.stride = 0
        self.matched = 0
        self.substituted = 0
        self.inserted = 0
        self.deleted = 0
        self.copying = 0
        self.worker_aligners = List[Aligner]()

    def fill(
        mut self,
        mut aligner: Aligner,
        mut realigned: _Realigned,
        consensus: List[UInt8],
        sequences: List[List[UInt8]],
        members: List[Int],
        anchors: List[Int],
        weights: List[Int] = List[Int](),
        threads: Int = 1,
    ) raises:
        """Every member aligned to the consensus: whole (`anchors` empty, or `BOTH`), or, a read anchored at one end,
        against a prefix (`LEFT`) or a suffix (`RIGHT`) of it. A path outside its alignment's span stays at its
        read's start before it and its end after it. The alignments come from `realigned` where it holds them."""
        realigned.on(consensus)
        var length = len(consensus)
        self.stride = length + 1
        var count = len(members)
        self.enter.resize(count * self.stride, 0)
        self.leave.resize(count * self.stride, 0)
        self.discord.resize(self.stride, 0)
        self.cover.resize(self.stride, 0)
        self.span_start.resize(count, 0)
        self.span_end.resize(count, 0)
        self.spanned = 0
        for i in range(self.stride):
            self.discord[i] = 0
            self.cover[i] = 0
        self.total = 0
        self.matched = 0
        self.substituted = 0
        self.inserted = 0
        self.deleted = 0
        self.copying = 0

        # Batch any uncached distinct alignments in parallel when multiple threads are requested.
        if threads > 1 and count >= 2:
            var missing = List[Int]()
            var queued = List[Bool](length=len(sequences), fill=False)
            for k in range(count):
                var first = realigned.same_as[members[k]]
                if not realigned.found[first] and not queued[first]:
                    queued[first] = True
                    missing.append(first)
            var workers = thread_count(threads, len(missing))
            if workers > 1 and len(missing) >= 2:
                while len(self.worker_aligners) < workers:
                    self.worker_aligners.append(Aligner(warm=False))
                var order = List[UInt64](capacity=len(missing))
                for idx in range(len(missing)):
                    order.append((UInt64(len(sequences[missing[idx]])) << 32) | UInt64(idx))
                sort(order)
                order.reverse()
                var results = List[Optional[Alignment]](length=len(missing), fill=None)
                var failed = List[Bool](length=len(missing), fill=False)
                var aligners_p = self.worker_aligners.unsafe_ptr()
                var results_p = results.unsafe_ptr()
                var failed_p = failed.unsafe_ptr()
                var taken = Atomic[Int64](0)

                def align_missing(
                    slot: Int,
                ) {
                    mut taken,
                    imm order,
                    imm missing,
                    imm consensus,
                    imm sequences,
                    imm anchors,
                    imm aligners_p,
                    imm results_p,
                    imm failed_p,
                }:
                    ref local_aligner = aligners_p[unsafe_offset=slot]
                    while True:
                        var pos = Int(taken.fetch_add(1))
                        if pos >= len(order):
                            return
                        var idx = Int(order[pos] & 0xFFFFFFFF)
                        var r = missing[idx]
                        var anchor = anchors[r] if len(anchors) > 0 else 0
                        try:
                            results_p[unsafe_offset=idx] = _align_one(local_aligner, consensus, sequences[r], anchor)
                        except:
                            failed_p[unsafe_offset=idx] = True

                spread(align_missing, workers, workers)
                for idx in range(len(missing)):
                    var r = missing[idx]
                    if failed[idx]:
                        var anchor = anchors[r] if len(anchors) > 0 else 0
                        realigned.found[r] = _align_one(aligner, consensus, sequences[r], anchor)
                    else:
                        realigned.found[r] = results[idx].take()

        for k in range(count):
            # A member standing for `weights[k]` reads alike counts that many times.
            var w = weights[k] if len(weights) > 0 else 1
            ref seq = sequences[members[k]]
            var anchor = anchors[members[k]] if len(anchors) > 0 else 0
            var first = realigned.same_as[members[k]]
            if not realigned.found[first]:
                realigned.found[first] = _align_one(aligner, consensus, seq, anchor)
            ref found = realigned.found[first].value()
            self.total += w * found.cost
            var base = k * self.stride
            var i = found.reference_start
            var j = found.query_start
            self.span_start[k] = i
            self.span_end[k] = found.reference_end
            self.spanned += w * (found.reference_end - found.reference_start)
            for b in range(i + 1):
                self.enter[base + b] = Int32(j)
                self.leave[base + b] = Int32(j)
            var cigar_bytes = found.cigar.as_bytes()
            var run_n = 0
            for ci in range(len(cigar_bytes)):
                var ch = cigar_bytes[ci]
                if ch >= 48 and ch <= 57:
                    run_n = run_n * 10 + Int(ch - 48)
                    continue
                var n = run_n if run_n > 0 else 1
                run_n = 0
                if ch == 73:  # 'I'
                    self.inserted += w * n
                    for _ in range(n):
                        if (i > 0 and seq[j] == consensus[i - 1]) or (i < length and seq[j] == consensus[i]):
                            self.copying += w
                        self.discord[i] += Int32(w)
                        j += 1
                elif ch == 68:  # 'D'
                    self.deleted += w * n
                    for _ in range(n):
                        self.leave[base + i] = Int32(j)
                        self.discord[i] += Int32(w)
                        i += 1
                        self.enter[base + i] = Int32(j)
                elif ch == 61 or ch == 77:  # '=' or 'M'
                    self.matched += w * n
                    for _ in range(n):
                        self.leave[base + i] = Int32(j)
                        j += 1
                        i += 1
                        self.enter[base + i] = Int32(j)
                elif ch == 88:  # 'X'
                    self.substituted += w * n
                    for _ in range(n):
                        self.leave[base + i] = Int32(j)
                        self.discord[i] += Int32(w)
                        j += 1
                        i += 1
                        self.enter[base + i] = Int32(j)
            self.leave[base + i] = Int32(j)
            for b in range(i + 1, self.stride):
                self.enter[base + b] = Int32(j)
                self.leave[base + b] = Int32(j)
            for b in range(self.span_start[k], i + 1):
                self.cover[b] += Int32(w)


struct _Bands(Movable):
    """One sequence's band in a segment: column c (consensus boundary s + c) spans read offsets `low[c]..=high[c]`,
    its cells from `start[c]` in a flat matrix."""

    var low: List[Int]
    var high: List[Int]
    var start: List[Int]

    def __init__(out self):
        self.low = List[Int]()
        self.high = List[Int]()
        self.start = List[Int]()

    def fill(mut self, paths: _Paths, k: Int, s: Int, e: Int, first: Int, last: Int):
        var columns = e - s + 1
        self.low.resize(columns, 0)
        self.high.resize(columns, 0)
        self.start.resize(columns + 1, 0)
        var base = k * paths.stride
        self.start[0] = 0
        for c in range(columns):
            self.low[c] = max(first, Int(paths.enter[base + s + c]) - BAND)
            self.high[c] = min(last, Int(paths.leave[base + s + c]) + BAND)
            self.start[c + 1] = self.start[c] + self.high[c] - self.low[c] + 1

    @inline(.always)
    def at(self, matrix: List[Int32], c: Int, j: Int) -> Int32:
        """Column c's cell at read offset j, INF outside the band."""
        if j < self.low[c] or j > self.high[c]:
            return INF
        return matrix[self.start[c] + j - self.low[c]]

    @inline(.always)
    def at_float(self, matrix: List[Float64], c: Int, j: Int) -> Float64:
        """Column c's probability at read offset j, 0 outside the band."""
        if j < self.low[c] or j > self.high[c]:
            return 0.0
        return matrix[self.start[c] + j - self.low[c]]


struct _EditWorker(Movable):
    """Per-worker reusable DP scratch for parallel edit-distance polishing."""

    var bands: _Bands
    var forward: List[Int32]
    var backward: List[Int32]
    var copied: List[Int32]
    var ahead: List[Int32]
    var gains: List[Int32]
    var periods: List[Int]
    var block_at: List[Int]
    var block_period: List[Int]
    var block_copy: List[Bool]
    var block_gain: List[Int32]

    def __init__(out self):
        self.bands = _Bands()
        self.forward = List[Int32]()
        self.backward = List[Int32]()
        self.copied = List[Int32]()
        self.ahead = List[Int32]()
        self.gains = List[Int32]()
        self.periods = List[Int]()
        self.block_at = List[Int]()
        self.block_period = List[Int]()
        self.block_copy = List[Bool]()
        self.block_gain = List[Int32]()


struct _LikelyWorker(Movable):
    """Per-worker reusable DP scratch for parallel likelihood polishing."""

    var bands: _Bands
    var forward: List[Float64]
    var slotted: List[Float64]
    var backward: List[Float64]
    var forward_scale: List[Float64]
    var backward_scale: List[Float64]
    var copied: List[Float64]
    var ahead: List[Float64]
    var gains: List[Float64]
    var periods: List[Int]
    var block_at: List[Int]
    var block_period: List[Int]
    var block_copy: List[Bool]
    var block_gain: List[Float64]

    def __init__(out self):
        self.bands = _Bands()
        self.forward = List[Float64]()
        self.slotted = List[Float64]()
        self.backward = List[Float64]()
        self.forward_scale = List[Float64]()
        self.backward_scale = List[Float64]()
        self.copied = List[Float64]()
        self.ahead = List[Float64]()
        self.gains = List[Float64]()
        self.periods = List[Int]()
        self.block_at = List[Int]()
        self.block_period = List[Int]()
        self.block_copy = List[Bool]()
        self.block_gain = List[Float64]()


struct PolishScratch(Movable):
    """Reusable workspace buffers for `polish`, `_lengthen_ties`, `_likely`, and `_runs` across calls."""

    var paths: _Paths
    var bands: _Bands
    var backward: List[Int32]
    var forward: List[Int32]
    var copied: List[Int32]
    var ahead: List[Int32]
    var gains: List[Int32]
    var site: List[Bool]
    var periods: List[Int]
    var block_at: List[Int]
    var block_period: List[Int]
    var block_copy: List[Bool]
    var block_gain: List[Int32]
    var likely_forward: List[Float64]
    var likely_slotted: List[Float64]
    var likely_backward: List[Float64]
    var likely_forward_scale: List[Float64]
    var likely_backward_scale: List[Float64]
    var likely_copied: List[Float64]
    var likely_ahead: List[Float64]
    var likely_gains: List[Float64]
    var likely_block_gain: List[Float64]
    var lengths: List[Int]
    var on: List[Int]
    var which: List[Int]
    var edit_workers: List[_EditWorker]
    var likely_workers: List[_LikelyWorker]

    def __init__(out self):
        self.paths = _Paths()
        self.bands = _Bands()
        self.backward = List[Int32]()
        self.forward = List[Int32]()
        self.copied = List[Int32]()
        self.ahead = List[Int32]()
        self.gains = List[Int32]()
        self.site = List[Bool]()
        self.periods = List[Int]()
        self.block_at = List[Int]()
        self.block_period = List[Int]()
        self.block_copy = List[Bool]()
        self.block_gain = List[Int32]()
        self.likely_forward = List[Float64]()
        self.likely_slotted = List[Float64]()
        self.likely_backward = List[Float64]()
        self.likely_forward_scale = List[Float64]()
        self.likely_backward_scale = List[Float64]()
        self.likely_copied = List[Float64]()
        self.likely_ahead = List[Float64]()
        self.likely_gains = List[Float64]()
        self.likely_block_gain = List[Float64]()
        self.lengths = List[Int]()
        self.on = List[Int]()
        self.which = List[Int]()
        self.edit_workers = List[_EditWorker]()
        self.likely_workers = List[_LikelyWorker]()


struct RepeatMask(Movable):
    """The stretches of a sequence inside tandem repeats: a k-mer start is masked when the same k-mer starts
    again within `distance` letters. An anchor any of whose k-mers is masked sits in a repeat, where an exact
    match, unique or not, may pair different copies."""

    var k: Int
    var masked_before: List[Int32]
    """How many masked k-mer starts precede each position."""

    def __init__(out self, sequence: ImmSpan[UInt8, _], k: Int, distance: Int):
        self.k = k
        var n = len(sequence)
        self.masked_before = List[Int32](length=n + 1, fill=0)
        if k <= 0 or n < k:
            return
        var marked = List[Bool](length=n, fill=False)
        var last = Dict[UInt64, Int]()
        var mask = (UInt64(1) << UInt64(2 * k)) - 1
        var key = UInt64(0)
        var valid = 0
        for i in range(n):
            var code = encode_base(sequence[i])
            valid = valid + 1 if code < 4 else 0
            key = ((key << 2) | UInt64(code & 3)) & mask
            if valid < k:
                continue
            var start = i - k + 1
            var seen = last.get(key)
            if seen and start - seen.value() <= distance:
                marked[start] = True
                marked[seen.value()] = True
            last[key] = start
        for i in range(n):
            self.masked_before[i + 1] = self.masked_before[i] + Int32(1 if marked[i] else 0)

    def covers(self, start: Int, end: Int) -> Bool:
        """Whether any k-mer of `start ..< end` is masked."""
        if self.k <= 0 or end - start < self.k:
            return False
        return self.masked_before[end - self.k + 1] > self.masked_before[start]


def _segments(consensus: List[UInt8], discord: List[Int32], mask: RepeatMask) -> List[Int]:
    """Segment boundaries 0 = b0 < b1 < ... = len: each at least SEGMENT letters on from the last and outside
    tandem repeats, where the fewest sequences disagree nearby among the next 32 such positions."""
    var length = len(consensus)
    var out: List[Int] = [0]
    var s = 0
    while s < length:
        var e = length
        if s + 2 * SEGMENT < length:
            var best = -1
            var best_discord = Int.MAX
            var seen = 0
            var p = s + SEGMENT
            while p < length - SEGMENT and seen < 32:
                if not mask.covers(max(0, p - REPEAT_K), min(length, p + REPEAT_K)):
                    var near = 0
                    for q in range(max(0, p - 4), min(length + 1, p + 5)):
                        near += Int(discord[q])
                    if near < best_discord:
                        best = p
                        best_discord = near
                    seen += 1
                p += 1
            if best >= 0:
                e = best
        out.append(e)
        s = e
    return out^


def _periods(consensus: List[UInt8], i: Int, s: Int, mut out: List[Int]):
    """The two shortest periods p (2..MAX_PERIOD) the consensus repeats with at boundary i: the p letters before
    it, inside the segment from s, and the p after it differ at most at p / 6 places."""
    out.clear()
    for p in range(2, MAX_PERIOD + 1):
        if i - p < s or i + p > len(consensus):
            break
        var differ = 0
        for t in range(p):
            if consensus[i - p + t] != consensus[i + t]:
                differ += 1
                if differ > p // 6:
                    break
        if differ <= p // 6:
            out.append(p)
            if len(out) == 2:
                return


def past_the_end(mut aligner: Aligner, consensus: List[UInt8], seq: List[UInt8], left: Bool) raises -> Alignment:
    """`seq`, anchored at the consensus's start (`left`) or end, aligned past the consensus's other end: the whole
    consensus against a prefix (or a suffix) of `seq`, the rest of `seq` free."""
    var turned = aligner.align(
        Span(seq), Span(consensus), Costs.edit(), Mode.PREFIX if left else Mode.SUFFIX, ties=Ties.LEFT
    )
    var cigar = List[UInt8](capacity=turned.cigar.byte_length())
    comptime INSERTED = UInt8(ord("I"))
    comptime DELETED = UInt8(ord("D"))
    for c in turned.cigar.as_bytes():
        cigar.append(DELETED if c == INSERTED else (INSERTED if c == DELETED else c))
    return Alignment(
        turned.cost,
        turned.score,
        String(unsafe_from_utf8=cigar),
        turned.query_start,
        turned.query_end,
        turned.reference_start,
        turned.reference_end,
    )


comptime LANES = 8
"""The int32 lanes of a segment's forward and backward columns (a column is some 70 cells: wider vectors leave
longer scalar tails)."""
comptime FAR = Int32(1 << 30)
"""Above every cell's score, INF and what it gathers along a column included: a lane no running minimum takes."""


@inline(.always)
def _up[k: Int](v: SIMD[DType.int32, LANES]) -> SIMD[DType.int32, LANES]:
    """`v` moved `k` lanes up (lane i + k holds lane i), `FAR` below."""
    return SIMD[DType.int32, LANES](FAR).join(v).slice[LANES, offset=LANES - k]()


@inline(.always)
def _down[k: Int](v: SIMD[DType.int32, LANES]) -> SIMD[DType.int32, LANES]:
    """`v` moved `k` lanes down (lane i holds lane i + k), `FAR` above."""
    return v.join(SIMD[DType.int32, LANES](FAR)).slice[LANES, offset=k]()


def _forward_column(mut forward: List[Int32], bands: _Bands, c: Int, seq: List[UInt8], letter: UInt8, first: Int):
    """Forward column c + 1 from column c."""
    var lo = bands.low[c + 1]
    var hi = bands.high[c + 1]
    var at = bands.start[c + 1]
    var before_lo = bands.low[c]
    var before_hi = bands.high[c]
    var before_at = bands.start[c]
    var cells = forward.unsafe_ptr()
    var letters = seq.unsafe_ptr()
    var core_lo = max(lo, before_lo + 1, first + 1)
    var core_hi = min(hi, before_hi)
    var j = lo
    while j <= hi:
        if j >= core_lo and j + LANES - 1 <= core_hi:
            var beside = cells.unsafe_load[width=LANES](before_at + j - before_lo)
            var diagonal = cells.unsafe_load[width=LANES](before_at + j - 1 - before_lo)
            var differ = letters.unsafe_load[width=LANES](j - 1).ne(letter).cast[DType.int32]()
            cells.unsafe_store(at + j - lo, min(beside + 1, diagonal + differ))
            j += LANES
            continue
        var v = bands.at(forward, c, j) + 1
        if j > first:
            v = min(v, bands.at(forward, c, j - 1) + Int32(1 if seq[j - 1] != letter else 0))
        cells[unsafe_offset=at + j - lo] = v
        j += 1
    var lane = iota[DType.int32, LANES]()
    var above = FAR
    j = lo
    while j + LANES - 1 <= hi:
        var x = cells.unsafe_load[width=LANES](at + j - lo)
        x = min(x, _up[1](x) + 1)
        x = min(x, _up[2](x) + 2)
        x = min(x, _up[4](x) + 4)
        x = min(x, SIMD[DType.int32, LANES](above) + lane)
        cells.unsafe_store(at + j - lo, x)
        above = x[LANES - 1] + 1
        j += LANES
    while j <= hi:
        cells[unsafe_offset=at + j - lo] = min(cells[unsafe_offset=at + j - lo], above)
        above = cells[unsafe_offset=at + j - lo] + 1
        j += 1


def _backward_column(
    mut backward: List[Int32],
    bands: _Bands,
    c: Int,
    seq: List[UInt8],
    letter: UInt8,
    first: Int,
    last: Int,
    free_start: Bool,
    free_end: Bool,
):
    """Backward column c from column c + 1."""
    var lo = bands.low[c]
    var hi = bands.high[c]
    var at = bands.start[c]
    var after_lo = bands.low[c + 1]
    var after_hi = bands.high[c + 1]
    var after_at = bands.start[c + 1]
    var cells = backward.unsafe_ptr()
    var letters = seq.unsafe_ptr()
    var core_lo = max(lo, after_lo)
    var core_hi = min(hi, after_hi - 1, last - 1)
    var j = lo
    while j <= hi:
        if j >= core_lo and j + LANES - 1 <= core_hi:
            var beside = cells.unsafe_load[width=LANES](after_at + j - after_lo)
            var diagonal = cells.unsafe_load[width=LANES](after_at + j + 1 - after_lo)
            var differ = letters.unsafe_load[width=LANES](j).ne(letter).cast[DType.int32]()
            cells.unsafe_store(at + j - lo, min(beside + 1, diagonal + differ))
            j += LANES
            continue
        var v = bands.at(backward, c + 1, j) + 1
        if j < last:
            v = min(v, bands.at(backward, c + 1, j + 1) + Int32(1 if seq[j] != letter else 0))
        cells[unsafe_offset=at + j - lo] = v
        j += 1
    if free_end and hi == last:
        cells[unsafe_offset=at + hi - lo] = 0
    var lane = iota[DType.int32, LANES]()
    var below = FAR
    j = hi + 1
    while j - LANES >= lo:
        j -= LANES
        var x = cells.unsafe_load[width=LANES](at + j - lo)
        x = min(x, _down[1](x) + 1)
        x = min(x, _down[2](x) + 2)
        x = min(x, _down[4](x) + 4)
        x = min(x, SIMD[DType.int32, LANES](below) + (LANES - 1 - lane))
        cells.unsafe_store(at + j - lo, x)
        below = x[0] + 1
    while j > lo:
        j -= 1
        cells[unsafe_offset=at + j - lo] = min(cells[unsafe_offset=at + j - lo], below)
        below = cells[unsafe_offset=at + j - lo] + 1
    if free_start and lo == first:
        cells[unsafe_offset=at] = min(cells[unsafe_offset=at], bands.at(backward, c + 1, first))
