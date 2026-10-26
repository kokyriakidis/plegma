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

from std.math import iota, log

from karan import Aligner, Alignment, Costs, Mode, Ties

from .alphabet import encode_base
from .graph import DELETION, INSERTION, MATCH, cigar_runs
from .homopolymer import FLANK_UNKNOWN, STRAND_UNKNOWN, HomopolymerModel
from .kernel import LEFT, RIGHT


comptime BAND = 16
"""Read letters either side of a sequence's alignment path that its forward and backward matrices span. 16 and 32
gave shasta2's chr21 assembly byte for byte the same (8 one error more), the matrices half as wide: polishing is most
of a short local assembly's time."""
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
comptime MAX_PERIOD = 24
"""Repeat units up to this long are deleted or copied whole (the two bands a unit apart overlap: 2 `BAND` > it)."""


def _code(letter: UInt8) -> Int:
    """A, C, G and T as 0..3 (upper or lower case); 4 for anything else."""
    # Byte values, not `ord` of a literal: that decodes the literal's UTF-8 each call, in the votes' inner loop.
    var lower = letter | 0x20
    if lower == 0x61:
        return 0
    if lower == 0x63:
        return 1
    if lower == 0x67:
        return 2
    if lower == 0x74:
        return 3
    return 4


struct _Realigned(Movable):
    """The sequences' alignments to one consensus, as `_Paths.fill` found them, through one aligner. A fill against
    the same consensus takes them as they are (an alignment depends on its two sequences alone), and a sequence with
    the letters and anchor of one before it takes that one's (`same_as`): the likelihood rounds mostly start on the
    consensus the edit rounds end on, and the runs on theirs, two thirds of the letters polishing aligned."""

    var aligner: Aligner
    var consensus: List[UInt8]
    var found: List[Optional[Alignment]]
    """Per sequence, its alignment to `consensus`, once found."""
    var same_as: List[Int]
    """Per sequence, the first member with its letters and anchor (itself if none before it)."""

    def __init__(out self, var same_as: List[Int]):
        self.aligner = Aligner(warm=False)
        self.consensus = List[UInt8]()
        self.found = List[Optional[Alignment]]()
        self.same_as = same_as^

    def on(mut self, consensus: List[UInt8]):
        """The alignments kept to `consensus`: those to another forgotten."""
        if len(self.found) == len(self.same_as) and self.consensus == consensus:
            return
        self.consensus = consensus.copy()
        self.found = List[Optional[Alignment]](length=len(self.same_as), fill=None)


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

    def fill(
        mut self,
        mut realigned: _Realigned,
        consensus: List[UInt8],
        sequences: List[List[UInt8]],
        members: List[Int],
        anchors: List[Int],
        weights: List[Int] = List[Int](),
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
        for k in range(count):
            # A member standing for `weights[k]` reads alike counts that many times.
            var w = weights[k] if len(weights) > 0 else 1
            ref seq = sequences[members[k]]
            var anchor = anchors[members[k]] if len(anchors) > 0 else 0
            var mode = Mode.GLOBAL
            # A read anchored at one end ends inside the consensus (the whole read against a prefix of it, or for a
            # right read a suffix) or past its end (the whole consensus against a prefix of the read, the rest of
            # the read free): both tried, the cheaper kept. Its length alone does not tell, as a read with
            # deletions may be the shorter and still reach past the consensus's end.
            if anchor == LEFT:
                mode = Mode.PREFIX
            elif anchor == RIGHT:
                mode = Mode.SUFFIX
            var first = realigned.same_as[members[k]]
            if not realigned.found[first]:
                var aligned = realigned.aligner.align(Span(consensus), Span(seq), Costs.edit(), mode, ties=Ties.LEFT)
                # Every consensus letter the read lacks costs the other way at least one: tried only where that
                # could still win, and searched only as far as the first one's cost.
                if (anchor == LEFT or anchor == RIGHT) and length - len(seq) < aligned.cost:
                    var other = past_the_end(realigned.aligner, consensus, seq, anchor == LEFT)
                    if other.cost < aligned.cost:
                        aligned = other^
                realigned.found[first] = aligned^
            ref found = realigned.found[first].value()
            self.total += w * found.cost
            var runs = cigar_runs(found.cigar)
            var base = k * self.stride
            var i = found.reference_start
            var j = found.query_start
            self.span_start[k] = i
            self.span_end[k] = found.reference_end
            self.spanned += w * (found.reference_end - found.reference_start)
            for b in range(i + 1):
                self.enter[base + b] = Int32(j)
                self.leave[base + b] = Int32(j)
            for r in range(len(runs.ops)):
                var op = runs.ops[r]
                if op == INSERTION:
                    self.inserted += w * runs.counts[r]
                elif op == DELETION:
                    self.deleted += w * runs.counts[r]
                elif op == MATCH:
                    self.matched += w * runs.counts[r]
                else:
                    self.substituted += w * runs.counts[r]
                for _ in range(runs.counts[r]):
                    if op == INSERTION:
                        if (i > 0 and seq[j] == consensus[i - 1]) or (i < length and seq[j] == consensus[i]):
                            self.copying += w
                        self.discord[i] += Int32(w)
                        j += 1
                        continue
                    self.leave[base + i] = Int32(j)
                    if op != MATCH:
                        self.discord[i] += Int32(w)
                    if op != DELETION:
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

    @always_inline
    def at(self, matrix: List[Int32], c: Int, j: Int) -> Int32:
        """Column c's cell at read offset j, INF outside the band."""
        if j < self.low[c] or j > self.high[c]:
            return INF
        return matrix[self.start[c] + j - self.low[c]]

    @always_inline
    def at_float(self, matrix: List[Float64], c: Int, j: Int) -> Float64:
        """Column c's probability at read offset j, 0 outside the band."""
        if j < self.low[c] or j > self.high[c]:
            return 0.0
        return matrix[self.start[c] + j - self.low[c]]


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
    for _ in range(rounds):
        paths.fill(realigned, current, sequences, unique, anchors, weight)
        unchecked = False
        if paths.total >= previous_total:
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
                    block_gain[b] += wk * (base - best)
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
    var out = _likely(current^, sequences, voters, likely_rounds, support, share, anchors, realigned)
    if homopolymers.loaded():
        out = _runs(out^, sequences, voters, anchors, homopolymers, strands, realigned)
    return out^


def _runs(
    var consensus: List[UInt8],
    sequences: List[List[UInt8]],
    members: List[Int],
    anchors: List[Int],
    model: HomopolymerModel,
    strands: List[Int],
    mut realigned: _Realigned,
) raises -> List[UInt8]:
    """Every run of two letters or more in `consensus` at its most likely length under `model`, given the length of
    the run each member reports: the read's letters aligned to the run and inserted at its two ends, with the letters
    like it next to them; a read with another letter inside the run, or not reaching past both its ends (anchored at
    one end), does not vote (as shasta2's msa1 repair votes: the reads with the run's symbol in its column)."""
    var length = len(consensus)
    var paths = _Paths()
    paths.fill(realigned, consensus, sequences, members, anchors)
    var out = List[UInt8](capacity=length + length // 8)
    var lengths = List[Int]()
    var on = List[Int]()
    var s = 0
    while s < length:
        var e = s + 1
        while e < length and consensus[e] == consensus[s]:
            e += 1
        var base = _code(consensus[s])
        if e - s < 2 or base > 3:
            for i in range(s, e):
                out.append(consensus[i])
            s = e
            continue
        var left = _code(consensus[s - 1]) if s > 0 else FLANK_UNKNOWN
        var right = _code(consensus[e]) if e < length else FLANK_UNKNOWN
        # A read votes when its alignment spans the run and a letter either side of it (where there is one).
        var need_start = s - 1 if s > 0 else 0
        var need_end = e + 1 if e < length else length
        lengths.clear()
        on.clear()
        for k in range(len(members)):
            if paths.span_start[k] > need_start or paths.span_end[k] < need_end:
                continue
            ref seq = sequences[members[k]]
            var kb = k * paths.stride
            var a = Int(paths.enter[kb + s])
            var b = Int(paths.leave[kb + e])
            # Letters unlike the run at its ends are the flanks' (inserted there, or the run's ends substituted).
            while a < b and _code(seq[a]) != base:
                a += 1
            while b > a and _code(seq[b - 1]) != base:
                b -= 1
            var pure = True
            for j in range(a, b):
                pure = pure and _code(seq[j]) == base
            if not pure:
                continue
            if a < b:
                while a > 0 and _code(seq[a - 1]) == base:
                    a -= 1
                while b < len(seq) and _code(seq[b]) == base:
                    b += 1
            lengths.append(b - a)
            var r = members[k]
            on.append(strands[r] if r < len(strands) else STRAND_UNKNOWN)
        var n = e - s
        if len(lengths) > 0:
            var best = model.most_likely(base, left, right, lengths, on)
            if best >= 2:
                n = best
        for _ in range(n):
            out.append(consensus[s])
        s = e
    return out^


def _groups(
    sequences: List[List[UInt8]],
    members: List[Int],
    anchors: List[Int],
    mut unique: List[Int],
    mut weight: List[Int],
    mut same_as: List[Int],
) raises:
    """`members` grouped by letters and anchor: the first of each group in `unique` (in the order they come), how
    many it stands for in `weight`, and per sequence the first of its group in `same_as` (itself if not a member)."""
    var at = Dict[String, Int]()
    same_as = List[Int](capacity=len(sequences))
    for r in range(len(sequences)):
        same_as.append(r)
    for r in members:
        var anchor = anchors[r] if len(anchors) > 0 else 0
        var key = String(anchor) + ":" + String(unsafe_from_utf8=sequences[r])
        var k = at.get(key, -1)
        if k < 0:
            at[key] = len(unique)
            unique.append(r)
            weight.append(1)
        else:
            weight[k] += 1
            same_as[r] = unique[k]


def past_the_end(mut aligner: Aligner, consensus: List[UInt8], seq: List[UInt8], left: Bool) raises -> Alignment:
    """`seq`, anchored at the consensus's start (`left`) or end, aligned past the consensus's other end: the whole
    consensus against a prefix (or a suffix) of `seq`, the rest of `seq` free. Aligned the other way round, `seq` as
    the reference, where karan's bit-parallel sweep serves it (free ends in the reference alone, no cap; the
    wavefront it replaces took most of polishing's time on reads past the end), then turned back."""
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
    """Forward column c + 1 from column c, as one cell at a time takes it: each cell the least of the cell beside it
    in column c plus one (the consensus letter deleted), the one diagonally before plus the two letters' mismatch
    (from the read's second letter on, `first` + 1), and the one above plus one (a read letter inserted). The first
    two from column c alone, a vector at a time where both are inside its band; the third a running minimum, the
    least over t <= j of a cell's own plus j - t, a vector at a time."""
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
    var above = FAR  # the cell above the block's first, plus one
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
    """Backward column c from column c + 1, as `_forward_column` builds forward ones, mirrored: each cell the least
    of the cell beside it in column c + 1 plus one, the one diagonally after plus the mismatch (before the read's
    end), and the one below plus one, a running minimum from the bottom. The read used up (`free_end`), its last cell
    is 0, which the cells above it see; not yet started (`free_start`), its first takes the cell beside it as it is
    (the consensus before the read's start free)."""
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
    var below = FAR  # the cell below the block's last, plus one
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


@always_inline
def _same(seq: List[UInt8], first: Int, last: Int, consensus: List[UInt8], s: Int, e: Int) -> Bool:
    """Whether read letters first ..< last are consensus letters s ..< e exactly."""
    if last - first != e - s:
        return False
    for t in range(e - s):
        if seq[first + t] != consensus[s + t]:
            return False
    return True


comptime MAX_COPYING = 0.7
"""Past this share of inserted letters copying a neighbour, the sequences' indels lengthen and shorten runs (as
nanopore reads' do), which the pair HMM's uniform insertions misjudge: on such read sets its edits made the
consensus worse, and `_likely` leaves it as it is. Uniform errors copy a neighbour some 45% of the time."""
comptime RUN_GAIN = 1.0
"""A likelihood edit lengthening a run (a letter inserted beside its like, or a one-letter unit copied) gains at
least this much: on simulated local assemblies 36 of the pair HMM's 38 errors the edit distance would have undone
were runs it lengthened, and with this floor their errors fell 88 -> 67 (0 and 8 lets every one or none through,
88 and 118)."""
comptime MIN_GAIN = 0.01
"""A likelihood edit gains at least this much (natural log, summed over the sequences)."""


@always_inline
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


def _likely(
    var current: List[UInt8],
    sequences: List[List[UInt8]],
    members: List[Int],
    rounds: Int,
    support: Int,
    share: Float64,
    anchors: List[Int],
    mut realigned: _Realigned,
) raises -> List[UInt8]:
    """The edits of `polish`, scored by the sequences' summed log-likelihood instead of their edit distance.

    The model, a pair HMM: before each template letter and after the last, k letters are inserted with probability
    pi^k (1 - pi), each uniform; a template letter is deleted with probability pd, else read as itself with
    probability 1 - ps, as each other letter ps / 3. Every alignment counts, so a read one letter short in a run of
    L counts L ways, where its edit distance counts one; on read sets whose edit distance prefers the consensus to
    the truth, this likelihood prefers the truth. The rates come from the sequences' alignments to the consensus,
    every round. Each sequence's forward and backward matrices, in probability, each column rescaled, span the
    band around its alignment that `polish`'s do, in the same segments; a segment takes its best edit alone a round
    (gains add up only nearly within one)."""
    var count = len(members)
    if count == 0 or rounds <= 0 or len(current) == 0:
        return current^
    var paths = _Paths()
    var bands = _Bands()
    var forward = List[Float64]()
    var slotted = List[Float64]()
    var backward = List[Float64]()
    var forward_scale = List[Float64]()
    var backward_scale = List[Float64]()
    var copied = List[Float64]()
    var ahead = List[Float64]()
    var gains = List[Float64]()
    var site = List[Bool]()
    var periods = List[Int]()
    var block_at = List[Int]()
    var block_period = List[Int]()
    var block_copy = List[Bool]()
    var block_gain = List[Float64]()
    var dirty = List[Bool](length=len(current) + 1, fill=True)
    comptime LETTERS = "ACGT"
    var before = List[UInt8]()
    var run_gain = RUN_GAIN
    for _ in range(rounds):
        paths.fill(realigned, current, sequences, members, anchors)
        if Float64(paths.copying) > MAX_COPYING * Float64(max(1, paths.inserted)):
            break
        var length = len(current)
        # The model's rates, each kept within [1e-4, 0.3].
        var letters = Float64(max(1, paths.spanned))
        var pd = min(0.3, max(1e-4, Float64(paths.deleted) / letters))
        var ps = min(0.3, max(1e-4, Float64(paths.substituted) / Float64(max(1, paths.matched + paths.substituted))))
        var per_slot = Float64(paths.inserted) / Float64(max(1, paths.spanned + count))
        var pi = min(0.3, max(1e-4, per_slot / (1.0 + per_slot)))
        var stay = pi / 4.0
        var leave = 1.0 - pi
        var kept = 1.0 - pd
        var right = 1.0 - ps
        var wrong = ps / 3.0
        _sites(paths, support, share, site)
        var mask = RepeatMask(Span(current), REPEAT_K, REPEAT_DISTANCE)
        var bounds = _segments(current, paths.discord, mask)
        var edit_start = List[Int]()
        var edit_end = List[Int]()
        var edit_text = List[List[UInt8]]()
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
            gains.resize((top - s) * SLOTS, 0.0)
            for x in range(len(gains)):
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
            block_gain.resize(len(block_at), 0.0)
            for b in range(len(block_gain)):
                block_gain[b] = 0.0
            # A read the segment holds only part of (anchored at one end, ending or starting inside) has no
            # likelihood here: the segment is left to the edit distance, which frees its end, rather than to the
            # reads spanning it alone, whose own errors would win.
            var partial = False
            for k in range(count):
                var overlaps = paths.span_end[k] > s and paths.span_start[k] < e
                partial = partial or (overlaps and (paths.span_start[k] > s or paths.span_end[k] < e))
            if partial:
                continue
            var columns = e - s + 1
            var c_last = columns - 1
            forward_scale.resize(columns, 0.0)
            backward_scale.resize(columns, 0.0)
            for k in range(count):
                if paths.span_start[k] > s or paths.span_end[k] < e:
                    continue
                ref seq = sequences[members[k]]
                var kb = k * paths.stride
                var first = Int(paths.enter[kb + s])
                var last = Int(paths.enter[kb + e]) if e < length else Int(paths.leave[kb + length])
                bands.fill(paths, k, s, e, first, last)
                var cells = bands.start[columns]
                forward.resize(cells, 0.0)
                slotted.resize(cells, 0.0)
                backward.resize(cells, 0.0)
                # Forward: column c holds the read's first ..< j generated by consensus s ..< s + c, before the
                # slot at s + c; `slotted` the same after it.
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
                    var most = 0.0
                    for j in range(nlo, nhi + 1):
                        var v = slotted[at + j - lo] * pd if j >= lo and j <= hi else 0.0
                        if j - 1 >= lo and j - 1 <= hi:
                            v += slotted[at + j - 1 - lo] * kept * (right if seq[j - 1] == letter else wrong)
                        forward[nat + j - nlo] = v
                        most = max(most, v)
                    if most <= 0:
                        most = 1e-300
                    for j in range(nlo, nhi + 1):
                        forward[nat + j - nlo] /= most
                    forward_scale[c + 1] = forward_scale[c] + log(most)
                # Backward: column c holds the read's j ..< last generated from the slot at s + c on.
                var lo_last = bands.low[c_last]
                var at_last = bands.start[c_last]
                for j in reversed(range(lo_last, bands.high[c_last] + 1)):
                    if e == length:
                        backward[at_last + j - lo_last] = (
                            leave if j == last else stay * backward[at_last + j + 1 - lo_last]
                        )
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
                    for j in range(lo, hi + 1):
                        backward[at + j - lo] /= most
                    backward_scale[c] = backward_scale[c + 1] + log(most)
                var base = _log(backward[bands.start[0] + first - bands.low[0]]) + backward_scale[0]
                for c in range(columns):
                    var i = s + c
                    if i >= top or not site[i]:
                        continue
                    var lo = bands.low[c]
                    var hi = bands.high[c]
                    var at = bands.start[c]
                    # Delete letter i: forward before slot i against backward from slot i + 1 (the two slots one).
                    # Substitute x: after slot i, x read as read letter j or deleted, then backward from slot i + 1.
                    # Insert x before i: after slot i, x, then backward from slot i (a new slot after x).
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
                # Whole units: deleted, forward before slot i - p against backward from slot i; copied, the units'
                # letters (each with its slot after) between slot i and backward from slot i again.
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
                    block_gain[b] += _log(total) + scale - base
            # The segment's best edit alone.
            var best = MIN_GAIN
            var best_x = -1
            for x in range(len(gains)):
                var slot = x % SLOTS
                # After the last letter only insertions are scored.
                if slot >= 4 and s + x // SLOTS >= length:
                    continue
                if slot >= 5 and _code(current[s + x // SLOTS]) == slot - 5:
                    continue
                # A letter lengthening a run only where it gains `RUN_GAIN` or more: the pair HMM's alignments multiply
                # with a run's length, so it favours runs longer than the reads hold.
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
            if best_x < 0:
                continue
            var text = List[UInt8]()
            if best_x < len(gains):
                var slot = best_x % SLOTS
                var start = s + best_x // SLOTS
                edit_start.append(start)
                edit_end.append(start if slot < 4 else start + 1)
                if slot < 4:
                    text.append(LETTERS.as_bytes()[slot])
                elif slot > 4:
                    text.append(LETTERS.as_bytes()[slot - 5])
            else:
                var b = best_x - len(gains)
                var i = block_at[b]
                var p = block_period[b]
                if block_copy[b]:
                    edit_start.append(i)
                    edit_end.append(i)
                    for t in range(i - p, i):
                        text.append(current[t])
                else:
                    edit_start.append(i - p)
                    edit_end.append(i)
            edit_text.append(text^)
        if len(edit_start) == 0:
            break
        # Segments are disjoint and listed left to right, and so are their edits.
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
            # Back to the consensus of two rounds ago: the edits only go back and forth.
            break
        before = current^
        current = next^
        dirty = List[Bool](length=len(current) + 1, fill=False)
        for m in marks:
            for q in range(max(0, m - 2 * SEGMENT), min(len(current) + 1, m + 2 * SEGMENT)):
                dirty[q] = True
    return current^
