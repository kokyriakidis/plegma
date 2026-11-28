"""Polishing's second stage: one-letter and whole-unit edits kept while they raise the reads' likelihood under a pair
HMM (`_likely`), the matrices filled in the same bands as the edit distance's."""

from std.atomic import Atomic
from std.math import log
from std.sys import simd_width_of

from karan import Aligner
from karan.common import spread, thread_count

from .bands import (
    _Bands,
    _LikelyWorker,
    MAX_PERIOD,
    _Paths,
    PolishScratch,
    REPEAT_DISTANCE,
    REPEAT_K,
    _Realigned,
    RepeatMask,
    SEGMENT,
    SLOTS,
    _code,
    _periods,
    _segments,
)

comptime MAX_COPYING = 0.7
"""Past this share of inserted letters copying a neighbour, the sequences' indels lengthen and shorten runs (as
nanopore reads' do), which the pair HMM's uniform insertions misjudge: on such read sets its edits made the
consensus worse, and `_likely` leaves it as it is. Uniform errors copy a neighbour some 45% of the time."""
comptime MIN_INSERTED = 20
"""The fewest inserted letters `MAX_COPYING`'s share is judged on: on clean reads one insertion copying its
neighbour was all of them, and the stage was left out."""
comptime RUN_GAIN = 1.0
"""A likelihood edit lengthening a run (a letter inserted beside its like, or a one-letter unit copied) gains at
least this much."""
comptime MIN_GAIN = 0.01
"""A likelihood edit gains at least this much (natural log, summed over the sequences)."""


@inline(.always)
def _log(x: Float64) -> Float64:
    return log(x) if x > 0 else -1e300


def _sites(paths: _Paths, support: Int, share: Float64, mut site: List[Bool]):
    """Per consensus boundary: whether, within two letters, `support` sequences disagree and `share` of those whose
    alignments span it."""
    var stride = paths.stride
    site.resize(stride, False)
    var near = 0
    for i in range(min(2, stride)):
        near += Int(paths.discord[i])
    for i in range(stride):
        if i + 2 < stride:
            near += Int(paths.discord[i + 2])
        if i - 3 >= 0:
            near -= Int(paths.discord[i - 3])
        site[i] = near >= max(support, Int(Float64(paths.cover[i]) * share + 0.999))


def _lengthens_run(consensus: List[UInt8], i: Int, code: Int) -> Bool:
    """Whether inserting the letter of `code` before consensus letter `i` lengthens a run of it."""
    return (i > 0 and _code(consensus[i - 1]) == code) or (i < len(consensus) and _code(consensus[i]) == code)


@fieldwise_init
struct _HmmRates(ImplicitlyCopyable, TrivialRegisterPassable):
    """Pair-HMM transition and emission probabilities for one likelihood round."""

    var pd: Float64
    var stay: Float64
    var leave: Float64
    var kept: Float64
    var right: Float64
    var wrong: Float64


struct _LikelyEdit(Copyable, Movable):
    """One segment's chosen likelihood edit (empty `has_edit=False` when none beats `MIN_GAIN`)."""

    var has_edit: Bool
    var start: Int
    var end: Int
    var text: List[UInt8]

    def __init__(out self):
        self.has_edit = False
        self.start = 0
        self.end = 0
        self.text = List[UInt8]()


def _score_likely_read(
    paths: _Paths,
    sequences: List[List[UInt8]],
    members: List[Int],
    current: List[UInt8],
    s: Int,
    e: Int,
    top: Int,
    length: Int,
    site: List[Bool],
    block_at: List[Int],
    block_period: List[Int],
    block_copy: List[Bool],
    rates: _HmmRates,
    k: Int,
    mut bands: _Bands,
    mut forward: List[Float64],
    mut slotted: List[Float64],
    mut backward: List[Float64],
    mut forward_scale: List[Float64],
    mut backward_scale: List[Float64],
    mut copied: List[Float64],
    mut ahead: List[Float64],
    mut gains: List[Float64],
    mut block_gain: List[Float64],
):
    """Computes read `k`'s forward/backward pair-HMM DP over segment `s..e` and adds its log-likelihood gains into
    `gains` and `block_gain`."""
    if paths.span_start[k] > s or paths.span_end[k] < e:
        return
    var columns = e - s + 1
    var c_last = columns - 1
    var pd = rates.pd
    var stay = rates.stay
    var leave = rates.leave
    var kept = rates.kept
    var right = rates.right
    var wrong = rates.wrong
    ref seq = sequences[members[k]]
    var kb = k * paths.stride
    var first = Int(paths.enter[kb + s])
    var last = Int(paths.enter[kb + e]) if e < length else Int(paths.leave[kb + length])
    bands.fill(paths, k, s, e, first, last)
    var cells = bands.start[columns]
    forward.resize(cells, 0.0)
    slotted.resize(cells, 0.0)
    backward.resize(cells, 0.0)
    forward_scale.resize(columns, 0.0)
    backward_scale.resize(columns, 0.0)
    for j in range(bands.low[0], bands.high[0] + 1):
        forward[j - bands.low[0]] = 1.0 if j == first else 0.0
    forward_scale[0] = 0.0
    for c in range(columns):
        var lo = bands.low[c]
        var hi = bands.high[c]
        var at = bands.start[c]
        var h = 0.0
        for j in range(lo, hi + 1):
            h = forward[at + j - lo] + stay * h
            slotted[at + j - lo] = leave * h
        if c == c_last:
            break
        var letter = current[s + c]
        var nlo = bands.low[c + 1]
        var nhi = bands.high[c + 1]
        var nat = bands.start[c + 1]
        comptime W = simd_width_of[DType.float64]()
        var fwd_p = forward.unsafe_ptr()
        var sl_p = slotted.unsafe_ptr()
        var seq_p = seq.unsafe_ptr()
        var v_pd = SIMD[DType.float64, W](pd)
        var v_kept = SIMD[DType.float64, W](kept)
        var v_right = SIMD[DType.float64, W](right)
        var v_wrong = SIMD[DType.float64, W](wrong)
        var core_lo = max(nlo, lo + 1)
        var core_hi = min(nhi, hi)
        var most = 0.0
        var most_v = SIMD[DType.float64, W](0.0)
        var j = nlo
        while j <= nhi:
            if j >= core_lo and j + W - 1 <= core_hi:
                var match_v = seq_p.unsafe_load[width=W](j - 1).eq(letter).select(v_right, v_wrong)
                var v = (
                    sl_p.unsafe_load[width=W](at + j - lo) * v_pd
                    + (sl_p.unsafe_load[width=W](at + j - 1 - lo) * v_kept) * match_v
                )
                fwd_p.unsafe_store(nat + j - nlo, v)
                most_v = max(most_v, v)
                j += W
                continue
            var v = sl_p[unsafe_offset=at + j - lo] * pd if j >= lo and j <= hi else 0.0
            if j - 1 >= lo and j - 1 <= hi:
                v += (
                    sl_p[unsafe_offset=at + j - 1 - lo]
                    * kept
                    * (right if seq_p[unsafe_offset=j - 1] == letter else wrong)
                )
            fwd_p[unsafe_offset=nat + j - nlo] = v
            most = max(most, v)
            j += 1
        most = max(most, most_v.reduce_max())
        if most <= 0:
            most = 1e-300
        var v_most = SIMD[DType.float64, W](most)
        j = nlo
        while j + W - 1 <= nhi:
            fwd_p.unsafe_store(nat + j - nlo, fwd_p.unsafe_load[width=W](nat + j - nlo) / v_most)
            j += W
        while j <= nhi:
            fwd_p[unsafe_offset=nat + j - nlo] /= most
            j += 1
        forward_scale[c + 1] = forward_scale[c] + log(most)
    var lo_last = bands.low[c_last]
    var at_last = bands.start[c_last]
    for j in reversed(range(lo_last, bands.high[c_last] + 1)):
        if e == length:
            backward[at_last + j - lo_last] = leave if j == last else stay * backward[at_last + j + 1 - lo_last]
        else:
            backward[at_last + j - lo_last] = 1.0 if j == last else 0.0
    backward_scale[c_last] = 0.0
    for c in reversed(range(c_last)):
        var letter = current[s + c]
        var lo = bands.low[c]
        var hi = bands.high[c]
        var at = bands.start[c]
        var nlo = bands.low[c + 1]
        var nhi = bands.high[c + 1]
        var nat = bands.start[c + 1]
        var b = 0.0
        var most = 0.0
        for j in reversed(range(lo, hi + 1)):
            var through = pd * backward[nat + j - nlo] if j >= nlo and j <= nhi else 0.0
            if j < last and j + 1 >= nlo and j + 1 <= nhi:
                through += kept * (right if seq[j] == letter else wrong) * backward[nat + j + 1 - nlo]
            b = leave * through + stay * b
            backward[at + j - lo] = b
            most = max(most, b)
        if most <= 0:
            most = 1e-300
        comptime WB = simd_width_of[DType.float64]()
        var bwd_p = backward.unsafe_ptr()
        var v_bmost = SIMD[DType.float64, WB](most)
        var bj = lo
        while bj + WB - 1 <= hi:
            bwd_p.unsafe_store(at + bj - lo, bwd_p.unsafe_load[width=WB](at + bj - lo) / v_bmost)
            bj += WB
        while bj <= hi:
            bwd_p[unsafe_offset=at + bj - lo] /= most
            bj += 1
        backward_scale[c] = backward_scale[c + 1] + log(most)
    var base = _log(backward[bands.start[0] + first - bands.low[0]]) + backward_scale[0]
    for c in range(columns):
        var i = s + c
        if i >= top or not site[i]:
            continue
        var lo = bands.low[c]
        var hi = bands.high[c]
        var at = bands.start[c]
        var deleted = 0.0
        var skipped = 0.0
        var read = 0.0
        var read_by = SIMD[DType.float64, 4](0)
        var skipped_here = 0.0
        var read_here = 0.0
        var read_here_by = SIMD[DType.float64, 4](0)
        for j in range(lo, hi + 1):
            var after = slotted[at + j - lo]
            var code = _code(seq[j]) if j < last else 4
            skipped_here += after * backward[at + j - lo]
            if j < last and j + 1 <= hi:
                var u = after * backward[at + j + 1 - lo]
                read_here += u
                if code < 4:
                    read_here_by[code] += u
            if i < e:
                var next = bands.at_float(backward, c + 1, j)
                deleted += forward[at + j - lo] * next
                skipped += after * next
                if j < last:
                    var t = after * bands.at_float(backward, c + 1, j + 1)
                    read += t
                    if code < 4:
                        read_by[code] += t
        var row = c * SLOTS
        var here = forward_scale[c] + backward_scale[c] - base
        for x in range(4):
            var p = pd * skipped_here + kept * (wrong * read_here + (right - wrong) * read_here_by[x])
            gains[row + x] += _log(p) + here
        if i < e:
            var onward = forward_scale[c] + backward_scale[c + 1] - base
            gains[row + 4] += _log(deleted) + onward
            for x in range(4):
                var p = pd * skipped + kept * (wrong * read + (right - wrong) * read_by[x])
                gains[row + 5 + x] += _log(p) + onward
    for b in range(len(block_at)):
        var c = block_at[b] - s
        var p = block_period[b]
        var lo = bands.low[c]
        var hi = bands.high[c]
        var at = bands.start[c]
        var total = 0.0
        var scale: Float64
        if not block_copy[b]:
            for j in range(lo, hi + 1):
                total += bands.at_float(forward, c - p, j) * backward[at + j - lo]
            scale = forward_scale[c - p] + backward_scale[c]
        else:
            var wide = min(last, hi + p)
            copied.resize(wide - lo + 1, 0.0)
            ahead.resize(wide - lo + 1, 0.0)
            for j in range(lo, wide + 1):
                copied[j - lo] = slotted[at + j - lo] if j <= hi else 0.0
            scale = forward_scale[c]
            for t in range(p):
                if t > 0:
                    var h = 0.0
                    for j in range(lo, wide + 1):
                        h = copied[j - lo] + stay * h
                        copied[j - lo] = leave * h
                var letter = current[s + c - p + t]
                var most = 0.0
                for j in range(lo, wide + 1):
                    var v = copied[j - lo] * pd
                    if j > lo:
                        v += copied[j - 1 - lo] * kept * (right if seq[j - 1] == letter else wrong)
                    ahead[j - lo] = v
                    most = max(most, v)
                if most <= 0:
                    most = 1e-300
                for j in range(lo, wide + 1):
                    ahead[j - lo] /= most
                scale += log(most)
                swap(copied, ahead)
            for j in range(lo, hi + 1):
                total += copied[j - lo] * backward[at + j - lo]
            scale += backward_scale[c]
        if total > 0:
            block_gain[b] += log(total) + scale - base
        else:
            block_gain[b] += Float64(p) * log(min(pd, stay * leave))


def _pick_likely_edit(
    current: List[UInt8],
    s: Int,
    length: Int,
    run_gain: Float64,
    gains: List[Float64],
    block_at: List[Int],
    block_period: List[Int],
    block_copy: List[Bool],
    block_gain: List[Float64],
) -> _LikelyEdit:
    """Chooses the single best likelihood edit in a segment (empty if none exceeds `MIN_GAIN`)."""
    comptime LETTERS = "ACGT"
    var best = MIN_GAIN
    var best_x = -1
    for x in range(len(gains)):
        var slot = x % SLOTS
        if slot >= 4 and s + x // SLOTS >= length:
            continue
        if slot >= 5 and _code(current[s + x // SLOTS]) == slot - 5:
            continue
        if slot < 4 and gains[x] < run_gain and _lengthens_run(current, s + x // SLOTS, slot):
            continue
        if gains[x] > best:
            best = gains[x]
            best_x = x
    for b in range(len(block_gain)):
        if block_copy[b] and block_period[b] == 1 and block_gain[b] < run_gain:
            continue
        if block_gain[b] > best:
            best = block_gain[b]
            best_x = len(gains) + b
    var out = _LikelyEdit()
    if best_x < 0:
        return out^
    out.has_edit = True
    if best_x < len(gains):
        var slot = best_x % SLOTS
        var start = s + best_x // SLOTS
        out.start = start
        out.end = start if slot < 4 else start + 1
        if slot < 4:
            out.text.append(LETTERS.as_bytes()[slot])
        elif slot > 4:
            out.text.append(LETTERS.as_bytes()[slot - 5])
    else:
        var b = best_x - len(gains)
        var i = block_at[b]
        var p = block_period[b]
        if block_copy[b]:
            out.start = i
            out.end = i
            for t in range(i - p, i):
                out.text.append(current[t])
        else:
            out.start = i - p
            out.end = i
    return out^


def _likely(
    mut aligner: Aligner,
    var current: List[UInt8],
    sequences: List[List[UInt8]],
    members: List[Int],
    rounds: Int,
    support: Int,
    share: Float64,
    anchors: List[Int],
    mut realigned: _Realigned,
    mut scratch: PolishScratch,
    threads: Int = 1,
) raises -> List[UInt8]:
    """The edits of `polish`, scored by the sequences' summed log-likelihood instead of their edit distance."""
    var count = len(members)
    if count == 0 or rounds <= 0 or len(current) == 0:
        return current^
    ref paths = scratch.paths
    ref bands = scratch.bands
    ref forward = scratch.likely_forward
    ref slotted = scratch.likely_slotted
    ref backward = scratch.likely_backward
    ref forward_scale = scratch.likely_forward_scale
    ref backward_scale = scratch.likely_backward_scale
    ref copied = scratch.likely_copied
    ref ahead = scratch.likely_ahead
    ref gains = scratch.likely_gains
    ref site = scratch.site
    ref periods = scratch.periods
    ref block_at = scratch.block_at
    ref block_period = scratch.block_period
    ref block_copy = scratch.block_copy
    ref block_gain = scratch.likely_block_gain
    var dirty = List[Bool](length=len(current) + 1, fill=True)
    var before = List[UInt8]()
    var run_gain = RUN_GAIN
    for _ in range(rounds):
        paths.fill(aligner, realigned, current, sequences, members, anchors, threads=threads)
        if paths.inserted >= MIN_INSERTED and Float64(paths.copying) > MAX_COPYING * Float64(paths.inserted):
            break
        var length = len(current)
        var letters = Float64(max(1, paths.spanned))
        var pd = min(0.3, max(1e-4, Float64(paths.deleted) / letters))
        var ps = min(0.3, max(1e-4, Float64(paths.substituted) / Float64(max(1, paths.matched + paths.substituted))))
        var per_slot = Float64(paths.inserted) / Float64(max(1, paths.spanned + count))
        var pi = min(0.3, max(1e-4, per_slot / (1.0 + per_slot)))
        var rates = _HmmRates(pd, pi / 4.0, 1.0 - pi, 1.0 - pd, 1.0 - ps, ps / 3.0)
        _sites(paths, support, share, site)
        var mask = RepeatMask(Span(current), REPEAT_K, REPEAT_DISTANCE)
        var bounds = _segments(current, paths.discord, mask)

        # Collect active segments (dirty, with a site, and with no read starting/ending inside).
        var active = List[Int]()
        for g in range(len(bounds) - 1):
            var s = bounds[g]
            var e = bounds[g + 1]
            var top = e + 1 if e == length else e
            var any_site = False
            var any_dirty = False
            for i in range(s, top):
                any_site = any_site or site[i]
                any_dirty = any_dirty or dirty[i]
            if not any_site or not any_dirty:
                continue
            var partial = False
            for k in range(count):
                var overlaps = paths.span_end[k] > s and paths.span_start[k] < e
                partial = partial or (overlaps and (paths.span_start[k] > s or paths.span_end[k] < e))
            if not partial:
                active.append(g)

        # Deduplicate identical members (same letters and anchor via `realigned.same_as`) so pair-HMM DP runs once
        # per distinct read and accumulates in exact ascending member order (100% Float64 bit-exact).
        var uniq_k = List[Int](capacity=count)
        var member_to_uniq = List[Int](length=count, fill=0)
        var seen_first = List[Int](length=len(sequences), fill=-1)
        for k in range(count):
            var first = realigned.same_as[members[k]]
            var u = seen_first[first]
            if u < 0:
                u = len(uniq_k)
                seen_first[first] = u
                uniq_k.append(k)
            member_to_uniq[k] = u
        var u_count = len(uniq_k)

        var edit_start = List[Int]()
        var edit_end = List[Int]()
        var edit_text = List[List[UInt8]]()
        var n_active = len(active)
        if threads > 1 and n_active >= 2 * threads:
            # Plenty of active segments: segment-level dispatch with per-segment unique-read deduplication.
            var seg_workers = thread_count(threads, n_active)
            while len(scratch.likely_workers) < seg_workers:
                scratch.likely_workers.append(_LikelyWorker())
            var seg_order = List[UInt64](capacity=n_active)
            for idx in range(n_active):
                var g = active[idx]
                seg_order.append((UInt64(bounds[g + 1] - bounds[g]) << 32) | UInt64(idx))
            sort(seg_order)
            seg_order.reverse()
            var seg_edits = List[_LikelyEdit](length=n_active, fill=_LikelyEdit())
            var seg_edits_p = seg_edits.unsafe_ptr()
            var pool = Span(scratch.likely_workers)
            var taken = Atomic[Int64](0)

            def score_segments(
                slot: Int,
            ) {
                mut taken,
                imm seg_order,
                imm active,
                imm bounds,
                imm length,
                imm site,
                imm current,
                imm paths,
                imm sequences,
                imm members,
                imm count,
                imm uniq_k,
                imm member_to_uniq,
                imm u_count,
                imm rates,
                imm run_gain,
                imm pool,
                imm seg_edits_p,
            }:
                ref w = pool.unsafe_ptr()[unsafe_offset=slot]
                while True:
                    var pos = Int(taken.fetch_add(1))
                    if pos >= len(seg_order):
                        return
                    var idx = Int(seg_order[pos] & 0xFFFFFFFF)
                    var g = active[idx]
                    var s = bounds[g]
                    var e = bounds[g + 1]
                    var top = e + 1 if e == length else e
                    var n_gains = (top - s) * SLOTS
                    w.gains.resize(n_gains, 0.0)
                    for x in range(n_gains):
                        w.gains[x] = 0.0
                    w.block_at.clear()
                    w.block_period.clear()
                    w.block_copy.clear()
                    for i in range(s, top):
                        if site[i]:
                            _periods(current, i, s, w.periods)
                            for p in w.periods:
                                for copy in range(2):
                                    w.block_at.append(i)
                                    w.block_period.append(p)
                                    w.block_copy.append(copy == 1)
                    var n_blocks = len(w.block_at)
                    w.block_gain.resize(n_blocks, 0.0)
                    for b in range(n_blocks):
                        w.block_gain[b] = 0.0
                    if u_count < count:
                        var active_cols = List[Int]()
                        for c in range(e - s + 1):
                            if s + c < top and site[s + c]:
                                active_cols.append(c)
                        var gstride = len(active_cols) * SLOTS
                        var u_gains = List[Float64](length=u_count * gstride, fill=0.0)
                        var u_blocks = List[Float64](length=u_count * n_blocks, fill=0.0)
                        for u in range(u_count):
                            for a_i in range(len(active_cols)):
                                var row = active_cols[a_i] * SLOTS
                                for si in range(SLOTS):
                                    w.gains[row + si] = 0.0
                            for b in range(n_blocks):
                                w.block_gain[b] = 0.0
                            _score_likely_read(
                                paths,
                                sequences,
                                members,
                                current,
                                s,
                                e,
                                top,
                                length,
                                site,
                                w.block_at,
                                w.block_period,
                                w.block_copy,
                                rates,
                                uniq_k[u],
                                w.bands,
                                w.forward,
                                w.slotted,
                                w.backward,
                                w.forward_scale,
                                w.backward_scale,
                                w.copied,
                                w.ahead,
                                w.gains,
                                w.block_gain,
                            )
                            var gb = u * gstride
                            for a_i in range(len(active_cols)):
                                var row = active_cols[a_i] * SLOTS
                                for si in range(SLOTS):
                                    u_gains[gb + a_i * SLOTS + si] = w.gains[row + si]
                            var bb = u * n_blocks
                            for b in range(n_blocks):
                                u_blocks[bb + b] = w.block_gain[b]
                        for a_i in range(len(active_cols)):
                            var row = active_cols[a_i] * SLOTS
                            for si in range(SLOTS):
                                w.gains[row + si] = 0.0
                        for b in range(n_blocks):
                            w.block_gain[b] = 0.0
                        for k in range(count):
                            var u = member_to_uniq[k]
                            var gb = u * gstride
                            for a_i in range(len(active_cols)):
                                var row = active_cols[a_i] * SLOTS
                                for si in range(SLOTS):
                                    w.gains[row + si] += u_gains[gb + a_i * SLOTS + si]
                            var bb = u * n_blocks
                            for b in range(n_blocks):
                                w.block_gain[b] += u_blocks[bb + b]
                    else:
                        for k in range(count):
                            _score_likely_read(
                                paths,
                                sequences,
                                members,
                                current,
                                s,
                                e,
                                top,
                                length,
                                site,
                                w.block_at,
                                w.block_period,
                                w.block_copy,
                                rates,
                                k,
                                w.bands,
                                w.forward,
                                w.slotted,
                                w.backward,
                                w.forward_scale,
                                w.backward_scale,
                                w.copied,
                                w.ahead,
                                w.gains,
                                w.block_gain,
                            )
                    seg_edits_p[unsafe_offset=idx] = _pick_likely_edit(
                        current, s, length, run_gain, w.gains, w.block_at, w.block_period, w.block_copy, w.block_gain
                    )

            spread(score_segments, seg_workers, seg_workers)
            for idx in range(n_active):
                if seg_edits[idx].has_edit:
                    edit_start.append(seg_edits[idx].start)
                    edit_end.append(seg_edits[idx].end)
                    edit_text.append(seg_edits[idx].text.copy())
        elif threads > 1 and n_active >= 1 and n_active * u_count >= 4:
            # Flatten `(active_seg, unique_read)` across all workers so 1..2*threads-1 active segments (including long
            # repeat stragglers) saturate every core. Each segment reduces in exact `0..count-1` member order (bit-exact).
            var total_tasks = n_active * u_count
            var flat_workers = thread_count(threads, total_tasks)
            while len(scratch.likely_workers) < flat_workers:
                scratch.likely_workers.append(_LikelyWorker())
            var seg_s = List[Int](capacity=n_active)
            var seg_e = List[Int](capacity=n_active)
            var seg_top = List[Int](capacity=n_active)
            var seg_cols = List[List[Int]](capacity=n_active)
            var seg_block_at = List[List[Int]](capacity=n_active)
            var seg_block_period = List[List[Int]](capacity=n_active)
            var seg_block_copy = List[List[Bool]](capacity=n_active)
            var seg_gstride = List[Int](capacity=n_active)
            var seg_bstride = List[Int](capacity=n_active)
            var seg_goff = List[Int](capacity=n_active + 1)
            var seg_boff = List[Int](capacity=n_active + 1)
            seg_goff.append(0)
            seg_boff.append(0)
            for idx in range(n_active):
                var g = active[idx]
                var s = bounds[g]
                var e = bounds[g + 1]
                var top = e + 1 if e == length else e
                seg_s.append(s)
                seg_e.append(e)
                seg_top.append(top)
                var cols = List[Int]()
                for c in range(e - s + 1):
                    if s + c < top and site[s + c]:
                        cols.append(c)
                var bat = List[Int]()
                var bper = List[Int]()
                var bcop = List[Bool]()
                for i in range(s, top):
                    if site[i]:
                        _periods(current, i, s, periods)
                        for p in periods:
                            for copy in range(2):
                                bat.append(i)
                                bper.append(p)
                                bcop.append(copy == 1)
                var gs = len(cols) * SLOTS
                var bs = len(bat)
                seg_gstride.append(gs)
                seg_bstride.append(bs)
                seg_goff.append(seg_goff[idx] + u_count * gs)
                seg_boff.append(seg_boff[idx] + u_count * bs)
                seg_cols.append(cols^)
                seg_block_at.append(bat^)
                seg_block_period.append(bper^)
                seg_block_copy.append(bcop^)
            var seg_order = List[UInt64](capacity=n_active)
            for idx in range(n_active):
                seg_order.append((UInt64(seg_e[idx] - seg_s[idx]) << 32) | UInt64(idx))
            sort(seg_order)
            seg_order.reverse()
            var seg_rank = List[Int](capacity=n_active)
            for pos in range(n_active):
                seg_rank.append(Int(seg_order[pos] & 0xFFFFFFFF))
            var per_read_gains = List[Float64](length=seg_goff[n_active], fill=0.0)
            var per_read_blocks = List[Float64](length=seg_boff[n_active], fill=0.0)
            var pool = Span(scratch.likely_workers)
            var prg_p = per_read_gains.unsafe_ptr()
            var prb_p = per_read_blocks.unsafe_ptr()
            var taken = Atomic[Int64](0)

            def score_likely_flat(
                slot: Int,
            ) {
                mut taken,
                imm total_tasks,
                imm u_count,
                imm uniq_k,
                imm seg_rank,
                imm seg_s,
                imm seg_e,
                imm seg_top,
                imm seg_cols,
                imm seg_block_at,
                imm seg_block_period,
                imm seg_block_copy,
                imm seg_gstride,
                imm seg_bstride,
                imm seg_goff,
                imm seg_boff,
                imm paths,
                imm sequences,
                imm members,
                imm current,
                imm length,
                imm site,
                imm rates,
                imm pool,
                imm prg_p,
                imm prb_p,
            }:
                ref w = pool.unsafe_ptr()[unsafe_offset=slot]
                while True:
                    var task = Int(taken.fetch_add(1))
                    if task >= total_tasks:
                        return
                    var idx = seg_rank[task // u_count]
                    var u = task % u_count
                    var s = seg_s[idx]
                    var e = seg_e[idx]
                    var top = seg_top[idx]
                    var n_gains = (top - s) * SLOTS
                    var bs = seg_bstride[idx]
                    w.gains.resize(n_gains, 0.0)
                    w.block_gain.resize(bs, 0.0)
                    ref cols = seg_cols[idx]
                    for a_i in range(len(cols)):
                        var row = cols[a_i] * SLOTS
                        for si in range(SLOTS):
                            w.gains[row + si] = 0.0
                    for b in range(bs):
                        w.block_gain[b] = 0.0
                    _score_likely_read(
                        paths,
                        sequences,
                        members,
                        current,
                        s,
                        e,
                        top,
                        length,
                        site,
                        seg_block_at[idx],
                        seg_block_period[idx],
                        seg_block_copy[idx],
                        rates,
                        uniq_k[u],
                        w.bands,
                        w.forward,
                        w.slotted,
                        w.backward,
                        w.forward_scale,
                        w.backward_scale,
                        w.copied,
                        w.ahead,
                        w.gains,
                        w.block_gain,
                    )
                    var gbase = seg_goff[idx] + u * seg_gstride[idx]
                    for a_i in range(len(cols)):
                        var row = cols[a_i] * SLOTS
                        for si in range(SLOTS):
                            prg_p[unsafe_offset=gbase + a_i * SLOTS + si] = w.gains[row + si]
                    var bbase = seg_boff[idx] + u * bs
                    for b in range(bs):
                        prb_p[unsafe_offset=bbase + b] = w.block_gain[b]

            spread(score_likely_flat, flat_workers, flat_workers)
            for idx in range(n_active):
                var s = seg_s[idx]
                var top = seg_top[idx]
                var n_gains = (top - s) * SLOTS
                var bs = seg_bstride[idx]
                var gs = seg_gstride[idx]
                var goff = seg_goff[idx]
                var boff = seg_boff[idx]
                ref cols = seg_cols[idx]
                gains.resize(n_gains, 0.0)
                for x in range(n_gains):
                    gains[x] = 0.0
                block_gain.resize(bs, 0.0)
                for b in range(bs):
                    block_gain[b] = 0.0
                for k in range(count):
                    var u = member_to_uniq[k]
                    var gbase = goff + u * gs
                    for a_i in range(len(cols)):
                        var row = cols[a_i] * SLOTS
                        for si in range(SLOTS):
                            gains[row + si] += per_read_gains[gbase + a_i * SLOTS + si]
                    var bbase = boff + u * bs
                    for b in range(bs):
                        block_gain[b] += per_read_blocks[bbase + b]
                var chosen = _pick_likely_edit(
                    current,
                    s,
                    length,
                    run_gain,
                    gains,
                    seg_block_at[idx],
                    seg_block_period[idx],
                    seg_block_copy[idx],
                    block_gain,
                )
                if chosen.has_edit:
                    edit_start.append(chosen.start)
                    edit_end.append(chosen.end)
                    edit_text.append(chosen.text.copy())
        else:
            for g in active:
                var s = bounds[g]
                var e = bounds[g + 1]
                var top = e + 1 if e == length else e
                var n_gains = (top - s) * SLOTS
                gains.resize(n_gains, 0.0)
                for x in range(n_gains):
                    gains[x] = 0.0
                block_at.clear()
                block_period.clear()
                block_copy.clear()
                for i in range(s, top):
                    if site[i]:
                        _periods(current, i, s, periods)
                        for p in periods:
                            for copy in range(2):
                                block_at.append(i)
                                block_period.append(p)
                                block_copy.append(copy == 1)
                var n_blocks = len(block_at)
                block_gain.resize(n_blocks, 0.0)
                for b in range(n_blocks):
                    block_gain[b] = 0.0
                if u_count < count:
                    var active_cols = List[Int]()
                    for c in range(e - s + 1):
                        if s + c < top and site[s + c]:
                            active_cols.append(c)
                    var gstride = len(active_cols) * SLOTS
                    var u_gains = List[Float64](length=u_count * gstride, fill=0.0)
                    var u_blocks = List[Float64](length=u_count * n_blocks, fill=0.0)
                    for u in range(u_count):
                        for a_i in range(len(active_cols)):
                            var row = active_cols[a_i] * SLOTS
                            for si in range(SLOTS):
                                gains[row + si] = 0.0
                        for b in range(n_blocks):
                            block_gain[b] = 0.0
                        _score_likely_read(
                            paths,
                            sequences,
                            members,
                            current,
                            s,
                            e,
                            top,
                            length,
                            site,
                            block_at,
                            block_period,
                            block_copy,
                            rates,
                            uniq_k[u],
                            bands,
                            forward,
                            slotted,
                            backward,
                            forward_scale,
                            backward_scale,
                            copied,
                            ahead,
                            gains,
                            block_gain,
                        )
                        var gb = u * gstride
                        for a_i in range(len(active_cols)):
                            var row = active_cols[a_i] * SLOTS
                            for si in range(SLOTS):
                                u_gains[gb + a_i * SLOTS + si] = gains[row + si]
                        var bb = u * n_blocks
                        for b in range(n_blocks):
                            u_blocks[bb + b] = block_gain[b]
                    for a_i in range(len(active_cols)):
                        var row = active_cols[a_i] * SLOTS
                        for si in range(SLOTS):
                            gains[row + si] = 0.0
                    for b in range(n_blocks):
                        block_gain[b] = 0.0
                    for k in range(count):
                        var u = member_to_uniq[k]
                        var gb = u * gstride
                        for a_i in range(len(active_cols)):
                            var row = active_cols[a_i] * SLOTS
                            for si in range(SLOTS):
                                gains[row + si] += u_gains[gb + a_i * SLOTS + si]
                        var bb = u * n_blocks
                        for b in range(n_blocks):
                            block_gain[b] += u_blocks[bb + b]
                else:
                    for k in range(count):
                        _score_likely_read(
                            paths,
                            sequences,
                            members,
                            current,
                            s,
                            e,
                            top,
                            length,
                            site,
                            block_at,
                            block_period,
                            block_copy,
                            rates,
                            k,
                            bands,
                            forward,
                            slotted,
                            backward,
                            forward_scale,
                            backward_scale,
                            copied,
                            ahead,
                            gains,
                            block_gain,
                        )
                var chosen = _pick_likely_edit(
                    current, s, length, run_gain, gains, block_at, block_period, block_copy, block_gain
                )
                if chosen.has_edit:
                    edit_start.append(chosen.start)
                    edit_end.append(chosen.end)
                    edit_text.append(chosen.text.copy())

        if len(edit_start) == 0:
            break
        var next = List[UInt8](capacity=length + len(edit_start) * MAX_PERIOD)
        var marks = List[Int]()
        var at = 0
        for t in range(len(edit_start)):
            for q in range(at, edit_start[t]):
                next.append(current[q])
            marks.append(len(next))
            next.extend(edit_text[t].copy())
            at = edit_end[t]
        for q in range(at, length):
            next.append(current[q])
        if next == before:
            break
        before = current^
        current = next^
        dirty = List[Bool](length=len(current) + 1, fill=False)
        for m in marks:
            for q in range(max(0, m - 2 * SEGMENT), min(len(current) + 1, m + 2 * SEGMENT)):
                dirty[q] = True
    return current^
