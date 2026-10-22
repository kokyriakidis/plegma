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
period of up to 24 letters, a unit deleted or copied whole is tried too: one letter of a missing unit alone may
cost the sequences more than it saves. Edits in different segments are taken together, since their gains add up
(the old alignment's path through each boundary stays open); within a segment, the best ones far apart, or the
best alone inside a tandem repeat, where one edit gains alike at every copy. Rounds go on, revisiting only the
segments the last round changed, until none gains.
"""

from karan import Aligner, Costs, Mode, Ties

from .backbone import RepeatMask
from .graph import DELETION, INSERTION, MATCH, cigar_runs


comptime BAND = 32
"""Read letters either side of a sequence's alignment path that its forward and backward matrices span."""
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
comptime MAX_PERIOD = 24
"""Repeat units up to this long are deleted or copied whole (a sequence's band spans the shift)."""


def _code(letter: UInt8) -> Int:
    """A, C, G and T as 0..3 (upper or lower case); 4 for anything else."""
    var lower = letter | 0x20
    if lower == UInt8(ord("a")):
        return 0
    if lower == UInt8(ord("c")):
        return 1
    if lower == UInt8(ord("g")):
        return 2
    if lower == UInt8(ord("t")):
        return 3
    return 4


def _text(bytes: ImmSpan[UInt8, _]) -> String:
    return String(unsafe_from_utf8=List[UInt8](bytes))


struct _Paths(Movable):
    """Every sequence's alignment to the consensus: at consensus boundary i (before letter i), sequence k's path
    enters at read offset `enter[k * stride + i]` and leaves at `leave[k * stride + i]`, the letters between
    inserted before letter i."""

    var enter: List[Int32]
    var leave: List[Int32]
    var discord: List[Int32]
    """Per consensus boundary i: how many sequences insert before letter i, or delete or substitute it."""
    var total: Int
    var stride: Int

    def __init__(out self):
        self.enter = List[Int32]()
        self.leave = List[Int32]()
        self.discord = List[Int32]()
        self.total = 0
        self.stride = 0

    def fill(
        mut self, mut aligner: Aligner, consensus: List[UInt8], sequences: List[List[UInt8]], members: List[Int]
    ) raises:
        var length = len(consensus)
        self.stride = length + 1
        var count = len(members)
        self.enter.resize(count * self.stride, 0)
        self.leave.resize(count * self.stride, 0)
        self.discord.resize(self.stride, 0)
        for i in range(self.stride):
            self.discord[i] = 0
        self.total = 0
        var reference = _text(Span(consensus))
        for k in range(count):
            ref seq = sequences[members[k]]
            var found = aligner.align(reference, _text(Span(seq)), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
            self.total += found.cost
            var runs = cigar_runs(found.cigar)
            var base = k * self.stride
            var i = 0
            var j = 0
            self.enter[base] = 0
            for r in range(len(runs.ops)):
                var op = runs.ops[r]
                for _ in range(runs.counts[r]):
                    if op == INSERTION:
                        self.discord[i] += 1
                        j += 1
                        continue
                    self.leave[base + i] = Int32(j)
                    if op != MATCH:
                        self.discord[i] += 1
                    if op != DELETION:
                        j += 1
                    i += 1
                    self.enter[base + i] = Int32(j)
            self.leave[base + length] = Int32(j)


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

    @always_inline
    def at(self, matrix: List[Int32], c: Int, j: Int) -> Int32:
        """Column c's cell at read offset j, INF outside the band."""
        if j < self.low[c] or j > self.high[c]:
            return INF
        return matrix[self.start[c] + j - self.low[c]]


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


def polish(
    consensus: List[UInt8],
    sequences: List[List[UInt8]],
    members: List[Int],
    rounds: Int,
    support: Int = 2,
    share: Float64 = 0.5,
) raises -> List[UInt8]:
    """`consensus` polished against the `members` of `sequences`, for at most `rounds` rounds: edits are tried at
    positions where, within two letters, `support` sequences disagree with it and as many as `share` of them (the
    reads' own errors alone rarely reach half: on real VNTR read sets, 0 to 0.7 mended the same errors)."""
    var current = consensus.copy()
    var count = len(members)
    if count == 0 or rounds <= 0:
        return current^
    var aligner = Aligner()
    var paths = _Paths()
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
    for _ in range(rounds):
        paths.fill(aligner, current, sequences, members)
        if paths.total >= previous_total:
            # The last round's edits together cost more than they saved: back to before them, one a segment.
            current = previous^
            previous = List[UInt8]()
            dirty = previous_dirty^
            previous_dirty = List[Bool]()
            one_each = True
            previous_total = Int.MAX
            continue
        previous_total = paths.total
        var length = len(current)
        var stride = length + 1
        # Positions where enough sequences disagree, within two letters.
        site.resize(stride, False)
        var needed = max(support, Int(Float64(count) * share + 0.999))
        var near = 0
        for i in range(min(2, stride)):
            near += Int(paths.discord[i])
        for i in range(stride):
            if i + 2 < stride:
                near += Int(paths.discord[i + 2])
            if i - 3 >= 0:
                near -= Int(paths.discord[i - 3])
            site[i] = near >= needed
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
                ref seq = sequences[members[k]]
                var kb = k * paths.stride
                var first = Int(paths.enter[kb + s])
                var last = Int(paths.enter[kb + e]) if e < length else Int(paths.leave[kb + length])
                bands.fill(paths, k, s, e, first, last)
                # Forward: column c at read offset j, the distance of consensus s ..< s + c to read first ..< j.
                forward.resize(bands.start[columns], INF)
                for j in range(bands.low[0], bands.high[0] + 1):
                    forward[j - bands.low[0]] = Int32(j - first)
                for c in range(c_last):
                    var letter = current[s + c]
                    var lo = bands.low[c + 1]
                    var at = bands.start[c + 1]
                    for j in range(lo, bands.high[c + 1] + 1):
                        var v = bands.at(forward, c, j) + 1
                        if j > first:
                            v = min(v, bands.at(forward, c, j - 1) + Int32(1 if seq[j - 1] != letter else 0))
                            if j - 1 >= lo:
                                v = min(v, forward[at + j - 1 - lo] + 1)
                        forward[at + j - lo] = v
                # Backward: column c at read offset j, the distance of consensus s + c ..< e to read j ..< last.
                backward.resize(bands.start[columns], INF)
                for j in range(bands.low[c_last], bands.high[c_last] + 1):
                    backward[bands.start[c_last] + j - bands.low[c_last]] = Int32(last - j)
                for c in reversed(range(c_last)):
                    var letter = current[s + c]
                    var lo = bands.low[c]
                    var hi = bands.high[c]
                    var at = bands.start[c]
                    for j in reversed(range(lo, hi + 1)):
                        var v = bands.at(backward, c + 1, j) + 1
                        if j < last:
                            if j + 1 <= hi:
                                v = min(v, backward[at + j + 1 - lo] + 1)
                            v = min(v, bands.at(backward, c + 1, j + 1) + Int32(1 if seq[j] != letter else 0))
                        backward[at + j - lo] = v
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
                    var row = c * SLOTS
                    for x in range(4):
                        gains[row + x] += base - min(kept + 1, min(inserted + 1, inserted_by[x]))
                    if i < e:
                        gains[row + 4] += base - deleted
                        for x in range(4):
                            gains[row + 5 + x] += base - min(deleted + 1, min(through + 1, through_by[x]))
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
                                if j > lo:
                                    v = min(v, copied[j - 1 - lo] + Int32(1 if seq[j - 1] != letter else 0))
                                    v = min(v, ahead[j - 1 - lo] + 1)
                                ahead[j - lo] = v
                            swap(copied, ahead)
                        for j in range(lo, hi + 1):
                            best = min(best, copied[j - lo] + backward[bands.start[c] + j - lo])
                    block_gain[b] += base - best
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
        dirty = List[Bool](length=len(current) + 1, fill=False)
        for m in marks:
            for q in range(max(0, m - 2 * SEGMENT), min(len(current) + 1, m + 2 * SEGMENT)):
                dirty[q] = True
    return current^
