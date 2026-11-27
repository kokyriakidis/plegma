"""Consensus polishing: one-letter and whole-repeat-unit edits to a consensus, each kept when it lowers the sequences'
summed edit distance to it, followed by pair-HMM likelihood polishing (`likelihood.mojo`) and empirical homopolymer
run-length repair (`run_lengths.mojo`)."""

from std.atomic import Atomic

from karan import Aligner
from karan.common import spread, thread_count

from .bands import (
    APART,
    INF,
    MAX_PERIOD,
    OUTLIER,
    REPEAT_DISTANCE,
    REPEAT_K,
    SEGMENT,
    SLOTS,
    PolishScratch,
    RepeatMask,
    _Bands,
    _EditWorker,
    _Paths,
    _Realigned,
    _backward_column,
    _code,
    _forward_column,
    _periods,
    _segments,
    past_the_end,
)
from .homopolymer import HomopolymerModel
from .likelihood import _likely, _sites
from .run_lengths import _groups, _lengthen_ties, _runs
from .trace import POLISH, Trace


@inline(.always)
def _same(seq: List[UInt8], first: Int, last: Int, consensus: List[UInt8], s: Int, e: Int) -> Bool:
    """Whether read letters first ..< last are consensus letters s ..< e exactly."""
    if last - first != e - s:
        return False
    for t in range(e - s):
        if seq[first + t] != consensus[s + t]:
            return False
    return True


struct _SegmentEdits(Copyable, Movable):
    """Edits chosen for one segment in an edit-distance polishing round."""

    var starts: List[Int]
    var ends: List[Int]
    var texts: List[List[UInt8]]

    def __init__(out self):
        self.starts = List[Int]()
        self.ends = List[Int]()
        self.texts = List[List[UInt8]]()


def _score_edit_read(
    paths: _Paths,
    sequences: List[List[UInt8]],
    unique: List[Int],
    weight: List[Int],
    current: List[UInt8],
    s: Int,
    e: Int,
    top: Int,
    length: Int,
    site: List[Bool],
    block_at: List[Int],
    block_period: List[Int],
    block_copy: List[Bool],
    k: Int,
    mut bands: _Bands,
    mut forward: List[Int32],
    mut backward: List[Int32],
    mut copied: List[Int32],
    mut ahead: List[Int32],
    mut gains: List[Int32],
    mut block_gain: List[Int32],
    gain_base: Int = 0,
    block_base: Int = 0,
):
    """Scores one unique read `k` across segment `s..e`, accumulating weighted edit-distance gains into `gains` and
    `block_gain` (at `gain_base` / `block_base`)."""
    if paths.span_end[k] <= s or paths.span_start[k] >= e:
        return
    var columns = e - s + 1
    var c_last = columns - 1
    var free_start = paths.span_start[k] > s
    var free_end = paths.span_end[k] < e
    ref seq = sequences[unique[k]]
    var wk = Int32(weight[k])
    var kb = k * paths.stride
    var first = Int(paths.enter[kb + s])
    var last = Int(paths.enter[kb + e]) if e < length else Int(paths.leave[kb + length])
    if not free_start and not free_end and _same(seq, first, last, current, s, e):
        for c in range(columns):
            var i = s + c
            if i >= top or not site[i]:
                continue
            var row = gain_base + c * SLOTS
            for x in range(4):
                gains[row + x] -= wk
            if i < e:
                gains[row + 4] -= wk
                var own = _code(current[i])
                for x in range(4):
                    gains[row + 5 + x] -= Int32(0) if x == own else wk
        for b in range(len(block_at)):
            block_gain[block_base + b] -= wk * Int32(block_period[b])
        return
    bands.fill(paths, k, s, e, first, last)
    forward.resize(bands.start[columns], INF)
    for j in range(bands.low[0], bands.high[0] + 1):
        forward[j - bands.low[0]] = Int32(j - first)
    for c in range(c_last):
        var letter = current[s + c]
        var lo = bands.low[c + 1]
        var at = bands.start[c + 1]
        _forward_column(forward, bands, c, seq, letter, first)
        if free_start and first >= lo and first <= bands.high[c + 1]:
            forward[at + first - lo] = 0
        if free_end and last >= lo and last <= bands.high[c + 1]:
            forward[at + last - lo] = min(forward[at + last - lo], bands.at(forward, c, last))
    backward.resize(bands.start[columns], INF)
    for j in range(bands.low[c_last], bands.high[c_last] + 1):
        backward[bands.start[c_last] + j - bands.low[c_last]] = Int32(last - j)
    for c in reversed(range(c_last)):
        _backward_column(backward, bands, c, seq, current[s + c], first, last, free_start, free_end)
    var base = backward[bands.start[0]]
    for c in range(columns):
        var i = s + c
        if i >= top or not site[i]:
            continue
        var lo = bands.low[c]
        var hi = bands.high[c]
        var at = bands.start[c]
        var kept = INF
        var inserted = INF
        var inserted_by = SIMD[DType.int32, 4](INF)
        var deleted = INF
        var through = INF
        var through_by = SIMD[DType.int32, 4](INF)
        for j in range(lo, hi + 1):
            var f = forward[at + j - lo]
            if f >= INF:
                continue
            kept = min(kept, f + backward[at + j - lo])
            if i < e:
                deleted = min(deleted, f + bands.at(backward, c + 1, j))
            if j < last:
                var code = _code(seq[j])
                var u = f + backward[at + j + 1 - lo] if j + 1 <= hi else INF
                inserted = min(inserted, u)
                if code < 4:
                    inserted_by[code] = min(inserted_by[code], u)
                if i < e:
                    var t = f + bands.at(backward, c + 1, j + 1)
                    through = min(through, t)
                    if code < 4:
                        through_by[code] = min(through_by[code], t)
        var unread = INF
        var unread_next = INF
        if free_end and last >= lo and last <= hi:
            var f = forward[at + last - lo]
            unread = f + backward[at + last - lo]
            if i < e:
                unread_next = f + bands.at(backward, c + 1, last)
        if free_start and first >= lo and first <= hi:
            var f = forward[at + first - lo]
            unread = min(unread, f + backward[at + first - lo])
            if i < e:
                unread_next = min(unread_next, f + bands.at(backward, c + 1, first))
        var row = gain_base + c * SLOTS
        for x in range(4):
            gains[row + x] += wk * (base - min(min(kept + 1, unread), min(inserted + 1, inserted_by[x])))
        if i < e:
            gains[row + 4] += wk * (base - deleted)
            for x in range(4):
                gains[row + 5 + x] += wk * (base - min(min(deleted + 1, unread_next), min(through + 1, through_by[x])))
    for b in range(len(block_at)):
        var c = block_at[b] - s
        var p = block_period[b]
        var lo = bands.low[c]
        var hi = bands.high[c]
        var best = INF
        if not block_copy[b]:
            for j in range(lo, hi + 1):
                best = min(best, bands.at(forward, c - p, j) + backward[bands.start[c] + j - lo])
        else:
            var wide = min(last, hi + p)
            copied.resize(wide - lo + 1, INF)
            ahead.resize(wide - lo + 1, INF)
            for j in range(lo, wide + 1):
                copied[j - lo] = bands.at(forward, c, j)
            for t in range(p):
                var letter = current[s + c - p + t]
                for j in range(lo, wide + 1):
                    var v = copied[j - lo] + 1
                    if free_start and j == first:
                        v = 0
                    elif free_end and j == last:
                        v = copied[j - lo]
                    if j > lo:
                        v = min(v, copied[j - 1 - lo] + Int32(1 if seq[j - 1] != letter else 0))
                        v = min(v, ahead[j - 1 - lo] + 1)
                    ahead[j - lo] = v
                swap(copied, ahead)
            for j in range(lo, hi + 1):
                best = min(best, copied[j - lo] + backward[bands.start[c] + j - lo])
        block_gain[block_base + b] += wk * (base - min(best, base + Int32(p)))


def _pick_segment_edits(
    current: List[UInt8],
    s: Int,
    alone: Bool,
    gains: List[Int32],
    block_at: List[Int],
    block_period: List[Int],
    block_copy: List[Bool],
    block_gain: List[Int32],
    gain_base: Int = 0,
    gain_len: Int = -1,
    block_base: Int = 0,
) -> _SegmentEdits:
    """Selects the best non-overlapping (`APART` apart, or single if `alone`) edits for segment starting at `s`."""
    comptime LETTERS = "ACGT"
    var n_gains = len(gains) if gain_len < 0 else gain_len
    var cur_len = len(current)
    var keys = List[UInt64]()
    for x in range(n_gains):
        var gx = gains[gain_base + x]
        var slot = x % SLOTS
        var pos = s + x // SLOTS
        if gx <= 0 or (slot >= 4 and pos >= cur_len) or (slot >= 5 and _code(current[pos]) == slot - 5):
            continue
        keys.append((UInt64((1 << 24) - min(Int(gx), (1 << 24) - 1)) << 32) | UInt64(x))
    for b in range(len(block_at)):
        var bg = block_gain[block_base + b]
        if bg > 0:
            var x = n_gains + b
            keys.append((UInt64((1 << 24) - min(Int(bg), (1 << 24) - 1)) << 32) | UInt64(x))
    sort(keys)
    var out = _SegmentEdits()
    var chosen = List[Int]()
    for key in keys:
        var x = Int(key & 0xFFFFFFFF)
        var start: Int
        var end: Int
        var text = List[UInt8]()
        if x < n_gains:
            var slot = x % SLOTS
            start = s + x // SLOTS
            end = start if slot < 4 else start + 1
            if slot < 4:
                text.append(LETTERS.as_bytes()[slot])
            elif slot > 4:
                text.append(LETTERS.as_bytes()[slot - 5])
        else:
            var b = x - n_gains
            var i = block_at[b]
            var p = block_period[b]
            if block_copy[b]:
                start = i
                end = i
                for t in range(i - p, i):
                    text.append(current[t])
            else:
                start = i - p
                end = i
        var far = True
        for q in chosen:
            far = far and abs(q - start) >= APART
        if not far:
            continue
        chosen.append(start)
        out.starts.append(start)
        out.ends.append(end)
        out.texts.append(text^)
        if alone:
            break
    return out^


def polish(
    consensus: List[UInt8],
    sequences: List[List[UInt8]],
    members: List[Int],
    rounds: Int,
    support: Int = 2,
    share: Float64 = 0.5,
    likely_rounds: Int = 20,
    anchors: List[Int] = List[Int](),
    homopolymers: HomopolymerModel = HomopolymerModel(),
    strands: List[Int] = List[Int](),
    threads: Int = 1,
) raises -> List[UInt8]:
    """`consensus` polished against the `members` of `sequences`, for at most `rounds` edit-distance rounds, then at
    most `likely_rounds` likelihood rounds, and (when `homopolymers` is loaded) empirical homopolymer run-length
    repair."""
    var aligner = Aligner(warm=False)
    var scratch = PolishScratch()
    return polish_with_scratch(
        aligner,
        scratch,
        consensus,
        sequences,
        members,
        rounds,
        support,
        share,
        likely_rounds,
        anchors,
        homopolymers,
        strands,
        threads,
    )


def polish_with_scratch(
    mut aligner: Aligner,
    mut scratch: PolishScratch,
    consensus: List[UInt8],
    sequences: List[List[UInt8]],
    members: List[Int],
    rounds: Int,
    support: Int = 2,
    share: Float64 = 0.5,
    likely_rounds: Int = 20,
    anchors: List[Int] = List[Int](),
    homopolymers: HomopolymerModel = HomopolymerModel(),
    strands: List[Int] = List[Int](),
    threads: Int = 1,
) raises -> List[UInt8]:
    """As `polish`, reusing `aligner` and `scratch` across calls."""
    var current = consensus.copy()
    if len(members) == 0:
        return current^
    var unique = List[Int]()
    var weight = List[Int]()
    var same_as = List[Int]()
    _groups(sequences, members, anchors, unique, weight, same_as)
    var count = len(unique)
    var realigned = _Realigned(same_as^)
    ref paths = scratch.paths
    paths.fill(aligner, realigned, current, sequences, unique, anchors, weight, threads=threads)
    var voters = members.copy()
    var kept = List[Int]()
    var kept_weight = List[Int]()
    var outlier = List[Bool](length=len(sequences), fill=False)
    for k in range(count):
        ref found = realigned.found[unique[k]].value()
        var extent = max(found.reference_end - found.reference_start, found.query_end - found.query_start)
        if Float64(found.cost) > OUTLIER * Float64(max(1, extent)):
            outlier[unique[k]] = True
        else:
            kept.append(unique[k])
            kept_weight.append(weight[k])
    var trace = Trace()
    if trace.on(POLISH):
        trace.line(
            POLISH,
            String(
                len(current),
                " letters, ",
                len(members),
                " reads (",
                count,
                " distinct), ",
                count - len(kept),
                " distinct outliers past ",
                OUTLIER,
                " edits a letter",
                "" if len(kept) > 0 else " (all: none left out)",
            ),
        )
    if len(kept) > 0 and len(kept) < count:
        unique = kept^
        weight = kept_weight^
        count = len(unique)
        voters.clear()
        for r in members:
            if not outlier[realigned.same_as[r]]:
                voters.append(r)
    ref bands = scratch.bands
    ref backward = scratch.backward
    ref forward = scratch.forward
    ref copied = scratch.copied
    ref ahead = scratch.ahead
    ref gains = scratch.gains
    ref site = scratch.site
    ref periods = scratch.periods
    ref block_at = scratch.block_at
    ref block_period = scratch.block_period
    ref block_copy = scratch.block_copy
    ref block_gain = scratch.block_gain
    var dirty = List[Bool](length=len(current) + 1, fill=True)
    var previous = List[UInt8]()
    var previous_dirty = List[Bool]()
    var previous_total = Int.MAX
    var one_each = False
    var unchecked = False
    for round in range(rounds):
        paths.fill(aligner, realigned, current, sequences, unique, anchors, weight, threads=threads)
        unchecked = False
        if paths.total >= previous_total:
            if trace.on(POLISH):
                trace.line(
                    POLISH,
                    String(
                        "round ",
                        round,
                        ": edit distance ",
                        paths.total,
                        ", not under ",
                        previous_total,
                        "; the round reverted, ",
                        "done" if one_each else "next one edit a segment",
                    ),
                )
            current = previous^
            previous = List[UInt8]()
            if one_each:
                break
            dirty = previous_dirty^
            previous_dirty = List[Bool]()
            one_each = True
            previous_total = Int.MAX
            continue
        previous_total = paths.total
        var length = len(current)
        _sites(paths, support, share, site)
        var mask = RepeatMask(Span(current), REPEAT_K, REPEAT_DISTANCE)
        var bounds = _segments(current, paths.discord, mask)

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
            if any_site and any_dirty:
                active.append(g)

        var edit_start = List[Int]()
        var edit_end = List[Int]()
        var edit_text = List[List[UInt8]]()
        var n_active = len(active)
        if threads > 1 and n_active >= 2 * threads:
            # Plenty of independent segments to saturate all cores without cross-worker reduction.
            var seg_workers = thread_count(threads, n_active)
            while len(scratch.edit_workers) < seg_workers:
                scratch.edit_workers.append(_EditWorker())
            var seg_order = List[UInt64](capacity=n_active)
            for idx in range(n_active):
                var g = active[idx]
                seg_order.append((UInt64(bounds[g + 1] - bounds[g]) << 32) | UInt64(idx))
            sort(seg_order)
            seg_order.reverse()
            var seg_edits = List[_SegmentEdits](length=n_active, fill=_SegmentEdits())
            var seg_edits_p = seg_edits.unsafe_ptr()
            var pool_p = scratch.edit_workers.unsafe_ptr()
            var taken = Atomic[Int64](0)

            def score_edit_segments(
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
                imm unique,
                imm weight,
                imm count,
                imm one_each,
                imm mask,
                imm pool_p,
                imm seg_edits_p,
            }:
                ref w = pool_p[unsafe_offset=slot]
                while True:
                    var pos = Int(taken.fetch_add(1))
                    if pos >= len(seg_order):
                        return
                    var idx = Int(seg_order[pos] & 0xFFFFFFFF)
                    var g = active[idx]
                    var s = bounds[g]
                    var e = bounds[g + 1]
                    var top = e + 1 if e == length else e
                    w.gains.resize((top - s) * SLOTS, 0)
                    for x in range(len(w.gains)):
                        w.gains[x] = 0
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
                    w.block_gain.resize(len(w.block_at), 0)
                    for b in range(len(w.block_gain)):
                        w.block_gain[b] = 0
                    for k in range(count):
                        _score_edit_read(
                            paths,
                            sequences,
                            unique,
                            weight,
                            current,
                            s,
                            e,
                            top,
                            length,
                            site,
                            w.block_at,
                            w.block_period,
                            w.block_copy,
                            k,
                            w.bands,
                            w.forward,
                            w.backward,
                            w.copied,
                            w.ahead,
                            w.gains,
                            w.block_gain,
                        )
                    var alone = one_each or mask.covers(s, e)
                    seg_edits_p[unsafe_offset=idx] = _pick_segment_edits(
                        current, s, alone, w.gains, w.block_at, w.block_period, w.block_copy, w.block_gain
                    )

            spread(score_edit_segments, seg_workers, seg_workers)
            for idx in range(n_active):
                ref se = seg_edits[idx]
                for t in range(len(se.starts)):
                    edit_start.append(se.starts[t])
                    edit_end.append(se.ends[t])
                    edit_text.append(se.texts[t].copy())
        elif threads > 1 and n_active >= 1 and n_active * count >= 4:
            # Flatten `(active_seg, read_k)` across all workers so 1..2*threads-1 active segments (including long
            # tandem-repeat stragglers) saturate every core. Exact Int32 reduction is 100% bit-exact.
            var total_tasks = n_active * count
            var flat_workers = thread_count(threads, total_tasks)
            while len(scratch.edit_workers) < flat_workers:
                scratch.edit_workers.append(_EditWorker())
            var seg_s = List[Int](capacity=n_active)
            var seg_e = List[Int](capacity=n_active)
            var seg_top = List[Int](capacity=n_active)
            var gain_off = List[Int](capacity=n_active + 1)
            var block_off = List[Int](capacity=n_active + 1)
            var seg_block_at = List[List[Int]](capacity=n_active)
            var seg_block_period = List[List[Int]](capacity=n_active)
            var seg_block_copy = List[List[Bool]](capacity=n_active)
            gain_off.append(0)
            block_off.append(0)
            for idx in range(n_active):
                var g = active[idx]
                var s = bounds[g]
                var e = bounds[g + 1]
                var top = e + 1 if e == length else e
                seg_s.append(s)
                seg_e.append(e)
                seg_top.append(top)
                gain_off.append(gain_off[idx] + (top - s) * SLOTS)
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
                block_off.append(block_off[idx] + len(bat))
                seg_block_at.append(bat^)
                seg_block_period.append(bper^)
                seg_block_copy.append(bcop^)
            var total_gains = gain_off[n_active]
            var total_blocks = block_off[n_active]
            # Order segments longest-first so repeat stragglers launch first across all workers.
            var seg_order = List[UInt64](capacity=n_active)
            for idx in range(n_active):
                seg_order.append((UInt64(seg_e[idx] - seg_s[idx]) << 32) | UInt64(idx))
            sort(seg_order)
            seg_order.reverse()
            var seg_rank = List[Int](capacity=n_active)
            for pos in range(n_active):
                seg_rank.append(Int(seg_order[pos] & 0xFFFFFFFF))
            var pool_p = scratch.edit_workers.unsafe_ptr()
            var taken = Atomic[Int64](0)

            def score_edit_flat(
                slot: Int,
            ) {
                mut taken,
                imm total_tasks,
                imm count,
                imm seg_rank,
                imm seg_s,
                imm seg_e,
                imm seg_top,
                imm gain_off,
                imm block_off,
                imm seg_block_at,
                imm seg_block_period,
                imm seg_block_copy,
                imm total_gains,
                imm total_blocks,
                imm paths,
                imm sequences,
                imm unique,
                imm weight,
                imm current,
                imm length,
                imm site,
                imm pool_p,
            }:
                ref w = pool_p[unsafe_offset=slot]
                w.gains.resize(total_gains, 0)
                for x in range(total_gains):
                    w.gains[x] = 0
                w.block_gain.resize(total_blocks, 0)
                for b in range(total_blocks):
                    w.block_gain[b] = 0
                while True:
                    var task = Int(taken.fetch_add(1))
                    if task >= total_tasks:
                        return
                    var idx = seg_rank[task // count]
                    var k = task % count
                    _score_edit_read(
                        paths,
                        sequences,
                        unique,
                        weight,
                        current,
                        seg_s[idx],
                        seg_e[idx],
                        seg_top[idx],
                        length,
                        site,
                        seg_block_at[idx],
                        seg_block_period[idx],
                        seg_block_copy[idx],
                        k,
                        w.bands,
                        w.forward,
                        w.backward,
                        w.copied,
                        w.ahead,
                        w.gains,
                        w.block_gain,
                        gain_off[idx],
                        block_off[idx],
                    )

            spread(score_edit_flat, flat_workers, flat_workers)
            gains.resize(total_gains, 0)
            for x in range(total_gains):
                gains[x] = 0
            block_gain.resize(total_blocks, 0)
            for b in range(total_blocks):
                block_gain[b] = 0
            for slot in range(flat_workers):
                ref w = scratch.edit_workers[slot]
                for x in range(total_gains):
                    gains[x] += w.gains[x]
                for b in range(total_blocks):
                    block_gain[b] += w.block_gain[b]
            for idx in range(n_active):
                var s = seg_s[idx]
                var e = seg_e[idx]
                var alone = one_each or mask.covers(s, e)
                var se = _pick_segment_edits(
                    current,
                    s,
                    alone,
                    gains,
                    seg_block_at[idx],
                    seg_block_period[idx],
                    seg_block_copy[idx],
                    block_gain,
                    gain_off[idx],
                    gain_off[idx + 1] - gain_off[idx],
                    block_off[idx],
                )
                for t in range(len(se.starts)):
                    edit_start.append(se.starts[t])
                    edit_end.append(se.ends[t])
                    edit_text.append(se.texts[t].copy())
        else:
            for g in active:
                var s = bounds[g]
                var e = bounds[g + 1]
                var top = e + 1 if e == length else e
                gains.resize((top - s) * SLOTS, 0)
                for x in range(len(gains)):
                    gains[x] = 0
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
                block_gain.resize(len(block_at), 0)
                for b in range(len(block_gain)):
                    block_gain[b] = 0
                for k in range(count):
                    _score_edit_read(
                        paths,
                        sequences,
                        unique,
                        weight,
                        current,
                        s,
                        e,
                        top,
                        length,
                        site,
                        block_at,
                        block_period,
                        block_copy,
                        k,
                        bands,
                        forward,
                        backward,
                        copied,
                        ahead,
                        gains,
                        block_gain,
                    )
                var alone = one_each or mask.covers(s, e)
                var se = _pick_segment_edits(current, s, alone, gains, block_at, block_period, block_copy, block_gain)
                for t in range(len(se.starts)):
                    edit_start.append(se.starts[t])
                    edit_end.append(se.ends[t])
                    edit_text.append(se.texts[t].copy())

        if trace.on(POLISH):
            trace.line(
                POLISH,
                String("round ", round, ": edit distance ", paths.total, ", ", len(edit_start), " edits taken"),
            )
        if len(edit_start) == 0:
            break
        var order = List[UInt64]()
        for t in range(len(edit_start)):
            order.append((UInt64(edit_start[t]) << 32) | UInt64(t))
        sort(order)
        var next = List[UInt8](capacity=length + len(order) * MAX_PERIOD)
        var marks = List[Int]()
        var at = 0
        for key in order:
            var t = Int(key & 0xFFFFFFFF)
            for q in range(at, edit_start[t]):
                next.append(current[q])
            marks.append(len(next))
            next.extend(edit_text[t].copy())
            at = edit_end[t]
        for q in range(at, length):
            next.append(current[q])
        previous = current^
        previous_dirty = dirty^
        current = next^
        unchecked = True
        dirty = List[Bool](length=len(current) + 1, fill=False)
        for m in marks:
            for q in range(max(0, m - 2 * SEGMENT), min(len(current) + 1, m + 2 * SEGMENT)):
                dirty[q] = True
    if unchecked:
        paths.fill(aligner, realigned, current, sequences, unique, anchors, weight, threads=threads)
        if paths.total >= previous_total:
            current = previous^
    if rounds > 0:
        current = _lengthen_ties(
            aligner, current^, sequences, unique, anchors, weight, realigned, scratch, threads=threads
        )
    if trace.on(POLISH):
        trace.line(POLISH, String("edit distance stage: ", len(current), " letters"))
    var out = _likely(
        aligner,
        current^,
        sequences,
        voters,
        likely_rounds,
        support,
        share,
        anchors,
        realigned,
        scratch,
        threads=threads,
    )
    if trace.on(POLISH):
        trace.line(POLISH, String("likelihood stage: ", len(out), " letters"))
    if homopolymers.loaded():
        out = _runs(
            aligner, out^, sequences, voters, anchors, homopolymers, strands, realigned, scratch, threads=threads
        )
        if trace.on(POLISH):
            trace.line(POLISH, String("homopolymer stage: ", len(out), " letters"))
    return out^
