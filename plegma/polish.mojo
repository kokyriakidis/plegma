"""Consensus polishing: one-letter edits to a consensus, each kept when it lowers the sequences' summed edit
distance to it.

A POA's consensus goes wrong where the sequences' own errors leave no column a majority: in a tandem repeat or a
short run, one missing letter is an insertion in one slot for some sequences, in the next slot or a substitution
for others, each placement a minority. Their edit distances to the consensus agree all the same. Every round aligns
each sequence to the consensus (karan, unit costs) and cuts the consensus into segments at positions outside tandem
repeats. In a segment holding a disagreement, each sequence's edit distance matrix is filled forward and backward in
a band around its alignment, which gives the sequence's distance after any one-letter substitution, deletion or
insertion exactly: the least sum of a forward cell, the edit and a backward cell. A sequence whose extra copy of a
repeat unit sits anywhere in the band counts for the edit that restores it. Where the consensus repeats with a
period of up to 16 letters, a unit deleted or copied whole is tried too: one letter of a missing unit alone may
cost the sequences more than it saves. Edits in different segments are taken together, since their gains add up
(the old alignment's path through each boundary stays open); within a segment, the best ones far apart, or the
best alone inside a tandem repeat, where one edit gains alike at every copy. Rounds go on, revisiting only the
segments the last round changed, until none gains.
"""

from .bands import (
    APART,
    _Bands,
    INF,
    MAX_PERIOD,
    OUTLIER,
    _Paths,
    REPEAT_DISTANCE,
    REPEAT_K,
    _Realigned,
    RepeatMask,
    SEGMENT,
    SLOTS,
    _backward_column,
    _code,
    _forward_column,
    _periods,
    _segments,
)
from .homopolymer import HomopolymerModel
from .likelihood import _likely, _sites
from .run_lengths import _groups, _lengthen_ties, _runs

from .trace import POLISH, Trace

# Re-exported for callers that import it from here.
from .bands import past_the_end


@always_inline
def _same(seq: List[UInt8], first: Int, last: Int, consensus: List[UInt8], s: Int, e: Int) -> Bool:
    """Whether read letters first ..< last are consensus letters s ..< e exactly."""
    if last - first != e - s:
        return False
    for t in range(e - s):
        if seq[first + t] != consensus[s + t]:
            return False
    return True


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
) raises -> List[UInt8]:
    """`consensus` polished against the `members` of `sequences` (each aligned whole, or, with `anchors` per
    sequence, a read anchored at one end against a prefix or a suffix of it), for at most `rounds` rounds of edits
    scored by
    edit distance, then at most `likely_rounds` scored by likelihood (see `_likely`): edits are tried at positions
    where, within two letters, `support` sequences disagree with it and as many as `share` of them (the reads' own
    errors alone rarely reach half: on real VNTR read sets, 0 to 0.7 mended the same errors). With a `homopolymers`
    model loaded, every run of two letters or more then takes its most likely length (`_runs`), each sequence on its
    strand in `strands` (0, 1 or `STRAND_UNKNOWN`; all unknown when empty)."""
    var current = consensus.copy()
    if len(members) == 0:
        return current^
    # The edit rounds count reads with the same letters and anchor once, weighed by how many they are (a local
    # assembler's reads repeat often); their sums are integers, the same either way. The likelihood rounds take
    # every read, their sums in floating point.
    var unique = List[Int]()
    var weight = List[Int]()
    var same_as = List[Int]()
    _groups(sequences, members, anchors, unique, weight, same_as)
    var count = len(unique)
    var realigned = _Realigned(same_as^)
    var paths = _Paths()
    # A sequence whose alignment to the consensus is more than `OUTLIER` edits is no noisy copy of it (an empty read
    # between the anchors, another allele's): every round would move toward it, the likelihood model's rates with
    # it, so polishing leaves it out (all are kept when every one is such).
    paths.fill(realigned, current, sequences, unique, anchors, weight)
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
    var bands = _Bands()
    var backward = List[Int32]()
    var forward = List[Int32]()
    var copied = List[Int32]()
    var ahead = List[Int32]()
    var gains = List[Int32]()
    var site = List[Bool]()
    var periods = List[Int]()
    # A segment's whole-unit edits: at boundary `block_at`, `block_period` letters deleted or copied.
    var block_at = List[Int]()
    var block_period = List[Int]()
    var block_copy = List[Bool]()
    var block_gain = List[Int32]()
    var dirty = List[Bool](length=len(current) + 1, fill=True)
    var previous = List[UInt8]()
    var previous_dirty = List[Bool]()
    var previous_total = Int.MAX
    var one_each = False
    comptime LETTERS = "ACGT"
    var unchecked = False
    for round in range(rounds):
        paths.fill(realigned, current, sequences, unique, anchors, weight)
        unchecked = False
        if paths.total >= previous_total:
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
            # The last round's edits together cost more than they saved: back to before them. Once, then taken one
            # a segment; when those cost more too, done (taking them again would only go back and forth).
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
        # Positions where enough of the sequences there disagree, within two letters.
        _sites(paths, support, share, site)
        var mask = RepeatMask(Span(current), REPEAT_K, REPEAT_DISTANCE)
        var bounds = _segments(current, paths.discord, mask)
        # The edits taken: letters start ..< end replaced by `edit_text`.
        var edit_start = List[Int]()
        var edit_end = List[Int]()
        var edit_text = List[List[UInt8]]()
        for g in range(len(bounds) - 1):
            var s = bounds[g]
            var e = bounds[g + 1]
            # Insertions before letter e belong to the next segment, except after the last letter.
            var top = e + 1 if e == length else e
            var any_site = False
            var any_dirty = False
            for i in range(s, top):
                any_site = any_site or site[i]
                any_dirty = any_dirty or dirty[i]
            if not any_site or not any_dirty:
                continue
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
            var columns = e - s + 1
            var c_last = columns - 1
            for k in range(count):
                # A read anchored at one end votes on the letters of the segment its alignment spans: the consensus
                # before its start (after its end) costs it nothing, so it neither votes those letters deleted nor
                # leaves the segment, a tandem array its reads end inside, without votes.
                if paths.span_end[k] <= s or paths.span_start[k] >= e:
                    continue
                var free_start = paths.span_start[k] > s
                var free_end = paths.span_end[k] < e
                ref seq = sequences[unique[k]]
                var wk = Int32(weight[k])
                var kb = k * paths.stride
                var first = Int(paths.enter[kb + s])
                var last = Int(paths.enter[kb + e]) if e < length else Int(paths.leave[kb + length])
                if not free_start and not free_end and _same(seq, first, last, current, s, e):
                    # The read spells the segment exactly: every one-letter edit costs it one (a letter substituted
                    # by itself nothing) and a unit of p letters deleted or copied p (a string and itself with one
                    # such edit are exactly that far apart), so its votes need no dynamic programming.
                    for c in range(columns):
                        var i = s + c
                        if i >= top or not site[i]:
                            continue
                        var row = c * SLOTS
                        for x in range(4):
                            gains[row + x] -= wk
                        if i < e:
                            gains[row + 4] -= wk
                            var own = _code(current[i])
                            for x in range(4):
                                gains[row + 5 + x] -= Int32(0) if x == own else wk
                    for b in range(len(block_at)):
                        block_gain[b] -= wk * Int32(block_period[b])
                    continue
                bands.fill(paths, k, s, e, first, last)
                # Forward: column c at read offset j, the distance of consensus s ..< s + c to read first ..< j.
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
                    # The read used up: the consensus after its end is free, so its cost stays as it was.
                    if free_end and last >= lo and last <= bands.high[c + 1]:
                        forward[at + last - lo] = min(forward[at + last - lo], bands.at(forward, c, last))
                # Backward: column c at read offset j, the distance of consensus s + c ..< e to read j ..< last.
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
                    # Delete letter i: forward and backward of the next column at one offset. Substitute x: x
                    # against read letter j (free where they agree) between forward j and backward of the next
                    # column at j + 1, or x alone. Insert x before i likewise, against this column.
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
                    # A letter inserted or substituted where the read has ended (or not yet begun) costs it nothing:
                    # the path through this column with the read used up (or not yet read), the consensus there free.
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
                    var row = c * SLOTS
                    for x in range(4):
                        gains[row + x] += wk * (base - min(min(kept + 1, unread), min(inserted + 1, inserted_by[x])))
                    if i < e:
                        gains[row + 4] += wk * (base - deleted)
                        for x in range(4):
                            gains[row + 5 + x] += wk * (
                                base - min(min(deleted + 1, unread_next), min(through + 1, through_by[x]))
                            )
                # Whole units: deleted, forward at i - p against backward at i; copied, forward at i carried
                # through the p letters again, then against backward at i.
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
                                # Copied letters past the read's end cost it nothing.
                                if free_end and j == last:
                                    v = copied[j - lo]
                                if j > lo:
                                    v = min(v, copied[j - 1 - lo] + Int32(1 if seq[j - 1] != letter else 0))
                                    v = min(v, ahead[j - 1 - lo] + 1)
                                ahead[j - lo] = v
                            swap(copied, ahead)
                        for j in range(lo, hi + 1):
                            best = min(best, copied[j - lo] + backward[bands.start[c] + j - lo])
                    # A unit of p letters deleted or copied changes the read's distance by p at most: its band
                    # may not reach both ends of the unit (best then INF, and nine such reads wrapped the gain to
                    # the top of 32 bits).
                    block_gain[b] += wk * (base - min(best, base + Int32(p)))
            # This segment's best edits: far apart, or the best alone inside a tandem repeat (or after a round
            # whose edits clashed).
            var alone = one_each or mask.covers(s, e)
            var keys = List[UInt64]()
            for x in range(len(gains)):
                var slot = x % SLOTS
                if gains[x] <= 0 or (slot >= 5 and _code(current[s + x // SLOTS]) == slot - 5):
                    continue
                keys.append((UInt64((1 << 24) - min(Int(gains[x]), (1 << 24) - 1)) << 32) | UInt64(x))
            for b in range(len(block_at)):
                if block_gain[b] > 0:
                    var x = len(gains) + b
                    keys.append((UInt64((1 << 24) - min(Int(block_gain[b]), (1 << 24) - 1)) << 32) | UInt64(x))
            sort(keys)
            var chosen = List[Int]()
            for key in keys:
                var x = Int(key & 0xFFFFFFFF)
                var start: Int
                var end: Int
                var text = List[UInt8]()
                if x < len(gains):
                    var slot = x % SLOTS
                    start = s + x // SLOTS
                    end = start if slot < 4 else start + 1
                    if slot < 4:
                        text.append(LETTERS.as_bytes()[slot])
                    elif slot > 4:
                        text.append(LETTERS.as_bytes()[slot - 5])
                else:
                    var b = x - len(gains)
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
                edit_start.append(start)
                edit_end.append(end)
                edit_text.append(text^)
                if alone:
                    break
        trace.line(
            POLISH,
            String("round ", round, ": edit distance ", paths.total, ", ", len(edit_start), " edits taken"),
        )
        if len(edit_start) == 0:
            break
        # The new consensus, the edits in place left to right; the segments near them are revisited next round.
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
        # The last round's edits, past the rounds allowed, kept only if they gained.
        paths.fill(realigned, current, sequences, unique, anchors, weight)
        if paths.total >= previous_total:
            current = previous^
    if rounds > 0:
        current = _lengthen_ties(current^, sequences, unique, anchors, weight, realigned, paths)
    trace.line(POLISH, String("edit distance stage: ", len(current), " letters"))
    var out = _likely(current^, sequences, voters, likely_rounds, support, share, anchors, realigned)
    trace.line(POLISH, String("likelihood stage: ", len(out), " letters"))
    if homopolymers.loaded():
        out = _runs(out^, sequences, voters, anchors, homopolymers, strands, realigned)
        trace.line(POLISH, String("homopolymer stage: ", len(out), " letters"))
    return out^
