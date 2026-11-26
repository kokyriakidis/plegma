"""Joining the left and right one-sided consensus sequences of an unspanned local assembly (`join_sides`).

Separated from `poa.mojo` so period detection (`repeat_period`), phase matching, overlap scoring, tandem-array copy
hypothesis scoring, and multi-threaded candidate evaluation can be tested independently.
"""

from std.atomic import Atomic

from karan import Aligner, Alignment, Costs, Mode
from karan.common import spread, thread_count

from .alphabet import encode_base
from .cigar import DELETION, INSERTION, MATCH
from .pinned import unique_pins
from .pins import pins
from .trace import POA, Trace

comptime JOIN_ENDS = 2000
"""Without pins between the two sides' consensus sequences, their last and first this many letters are aligned for an
overlap."""
comptime JOIN_LEAST = 200
"""The fewest letters an overlap joining the two sides spans."""
comptime JOIN_IDENTITY = 0.75
"""The least share of an overlap's letters a karan overlap matches."""
comptime JOIN_EDIT = 19
"""What an edit costs in that overlap, a match earning 1."""
comptime JOIN_EDIT_STRICT = 100
"""The same, for the overlap of the ends found again with edits dearer still."""
comptime JOIN_MAX_PERIOD = 500
comptime JOIN_PERIOD_K = 8
comptime JOIN_PERIOD_SHARE = 0.3
"""An overlap whose `JOIN_PERIOD_K`-mers recur one period on (at most `JOIN_MAX_PERIOD`) this often is a tandem
array's (`repeat_period`)."""
comptime JOIN_SHIFTS = 64
"""The most copies either way a placement inside a tandem array is shifted by, to be scored too."""
comptime JOIN_DRIFT = 16
"""Pins whose diagonals are within this many letters of each other place right alike."""
comptime JOIN_MARGIN = 1.5
comptime JOIN_MARGIN_FLOOR = 0.002
"""A placement on another diagonal agreeing within `JOIN_MARGIN` times the best's edits a letter, plus
`JOIN_MARGIN_FLOOR`, leaves the two sides unjoined."""


def _support(reach: List[Int], at: Int) -> Int:
    """How many reads reaching `reach` letters from their anchor get past `at` letters from it."""
    var n = 0
    for length in reach:
        if length > at:
            n += 1
    return n


@fieldwise_init
struct _Placing(Copyable, Movable):
    """A way to join the two sides: right's letters from `cut_right` on after left's first `cut_left` (right
    starting at `diagonal` in left's letters there); how well the two agree over the overlap it implies (`edits` in
    `letters`), and where it puts each of left's letters in right (`path`)."""

    var diagonal: Int
    var cut_left: Int
    var cut_right: Int
    var edits: Int
    var letters: Int
    var how: String
    var path: List[Int]
    var path_from: Int
    var edits_after: List[Int]

    def edits_from(self, at_left: Int) -> Int:
        """This placement's edits at left's letters from `at_left` on."""
        var i = at_left - self.path_from
        if i <= 0:
            return self.edits_after[0] if len(self.edits_after) > 0 else self.edits
        return self.edits_after[min(i, len(self.edits_after) - 1)]

    def lands(self, at_left: Int) -> Int:
        """The right letter left's letter `at_left` lands on (outside the overlap: as the diagonal puts it)."""
        var k = at_left - self.path_from
        if k >= 0 and k < len(self.path):
            return self.path[k]
        return at_left - self.diagonal

    def rate(self) -> Float64:
        return Float64(self.edits) / Float64(max(1, self.letters))


@fieldwise_init
struct _Candidate(Copyable, Movable):
    """A candidate join meeting point `(cut_left, cut_right)` and its trace label."""

    var cut_left: Int
    var cut_right: Int
    var how: String


def _walk_and_tally(found: Alignment, right_from: Int, first: Int, mut path: List[Int], mut at: List[Int]):
    """Walks `found.cigar` directly in one pass without allocation, appending right coordinates into `path` and
    recording edit counts per left position into `at`."""
    var j = right_from + found.reference_start
    var i = first
    var at_last = len(at) - 1
    var run_n = 0
    for ch in found.cigar.as_bytes():
        if ch >= 48 and ch <= 57:
            run_n = run_n * 10 + Int(ch - 48)
            continue
        var n = run_n if run_n > 0 else 1
        run_n = 0
        if ch == DELETION:
            j += n
            at[min(i, at_last)] += n
        elif ch == INSERTION:
            for _ in range(n):
                path.append(j)
                at[i] += 1
                i += 1
        elif ch == MATCH:
            for _ in range(n):
                path.append(j)
                j += 1
                i += 1
        else:
            for _ in range(n):
                path.append(j)
                at[i] += 1
                j += 1
                i += 1


def _placing(
    mut aligner: Aligner, left: List[UInt8], right: List[UInt8], cut_left: Int, cut_right: Int, var how: String
) raises -> _Placing:
    """Right placed so that left's letter `cut_left` meets right's `cut_right`."""
    var lo = max(0, cut_left - cut_right)
    var out = _Placing(
        cut_left - cut_right,
        cut_left,
        cut_right,
        0,
        len(left) - lo,
        how^,
        List[Int](capacity=len(left) - lo),
        lo,
        List[Int](),
    )
    var at = List[Int](length=len(left) - lo + 1, fill=0)
    if cut_left > lo:
        var back = aligner.align(Span(right)[0:cut_right], Span(left)[lo:cut_left], Costs.edit(), Mode.SUFFIX)
        out.edits += back.cost
        _walk_and_tally(back, 0, 0, out.path, at)
    if cut_left < len(left):
        var on = aligner.align(Span(right)[cut_right:], Span(left)[cut_left:], Costs.edit(), Mode.PREFIX)
        out.edits += on.cost
        _walk_and_tally(on, cut_right, cut_left - lo, out.path, at)
    for k in reversed(range(len(at) - 1)):
        at[k] += at[k + 1]
    out.edits_after = at^
    return out^


def _eval_candidates(
    mut aligner: Aligner,
    left: List[UInt8],
    right: List[UInt8],
    candidates: List[_Candidate],
    threads: Int,
) -> List[_Placing]:
    """Evaluates all join candidates, in parallel when `threads > 1`, preserving candidate order and dropping any
    candidate whose alignment raises."""
    var n = len(candidates)
    var placings = List[_Placing]()
    if n == 0:
        return placings^
    var workers = thread_count(threads, n)
    if workers > 1 and n >= 2:
        var results = List[Optional[_Placing]](length=n, fill=None)
        var results_p = results.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def eval_worker(slot: Int) {mut taken, imm n, imm candidates, imm left, imm right, imm results_p}:
            var local_aligner = Aligner(warm=False)
            while True:
                var idx = Int(taken.fetch_add(1))
                if idx >= n:
                    return
                try:
                    results_p[unsafe_offset=idx] = _placing(
                        local_aligner,
                        left,
                        right,
                        candidates[idx].cut_left,
                        candidates[idx].cut_right,
                        candidates[idx].how.copy(),
                    )
                except:
                    pass

        spread(eval_worker, workers, workers)
        for idx in range(n):
            if results[idx]:
                placings.append(results[idx].take())
        return placings^

    for idx in range(n):
        try:
            placings.append(
                _placing(
                    aligner,
                    left,
                    right,
                    candidates[idx].cut_left,
                    candidates[idx].cut_right,
                    candidates[idx].how.copy(),
                )
            )
        except:
            pass
    return placings^


def _copy_join(
    left: List[UInt8], right: List[UInt8], hypotheses: List[_Placing], period: Int, trace: Trace
) -> List[UInt8]:
    """The two sides joined at the copy number whose placement has the fewest edits a letter, when every other copy
    number has more than `JOIN_MARGIN` times as many plus `JOIN_MARGIN_FLOOR`, on the left letters both align."""
    var best = _best(hypotheses)
    if best < 0:
        trace.line(POA, String("not joined: no copy number ", JOIN_IDENTITY, " alike over ", JOIN_LEAST, " letters"))
        return List[UInt8]()
    ref chosen = hypotheses[best]
    for c in range(len(hypotheses)):
        ref h = hypotheses[c]
        if c == best or h.letters < JOIN_LEAST or abs(h.diagonal - chosen.diagonal) < period // 2:
            continue
        var start = max(h.path_from, chosen.path_from)
        var e = _versus(h, chosen)
        var letters = Float64(max(1, len(left) - start))
        var mine = Float64(e[1]) / letters
        var theirs = Float64(e[0]) / letters
        if theirs <= mine * JOIN_MARGIN + JOIN_MARGIN_FLOOR:
            trace.line(
                POA,
                String(
                    "not joined: ",
                    chosen.how,
                    " has ",
                    e[1],
                    " edits, ",
                    h.how,
                    " ",
                    e[0],
                    " in the ",
                    len(left) - start,
                    " left letters both align",
                ),
            )
            return List[UInt8]()
    trace.line(POA, String("join at ", chosen.how, ": ", chosen.edits, " edits in ", chosen.letters, " letters"))
    var out = List[UInt8](capacity=chosen.cut_left + len(right) - chosen.cut_right)
    for i in range(chosen.cut_left):
        out.append(left[i])
    for j in range(chosen.cut_right, len(right)):
        out.append(right[j])
    return out^


def _best(placings: List[_Placing]) -> Int:
    """The placement with the fewest edits a letter (of two as good, the longer overlap), among those over `JOIN_LEAST`
    letters and `JOIN_IDENTITY` alike; -1 for none."""
    var best = -1
    for c in range(len(placings)):
        ref p = placings[c]
        if p.letters < JOIN_LEAST or p.rate() > 1.0 - JOIN_IDENTITY:
            continue
        if (
            best < 0
            or p.rate() < placings[best].rate()
            or (p.rate() == placings[best].rate() and p.letters > placings[best].letters)
        ):
            best = c
    return best


def repeat_period(seq: List[UInt8], start: Int, end: Int) -> Int:
    """The period `seq[start:end]` repeats with, or 0 for none."""
    var n = end - start
    if n < 2 * JOIN_PERIOD_K:
        return 0
    var at_distance = List[Int](length=JOIN_MAX_PERIOD + 1, fill=0)
    var last = Dict[UInt64, Int]()
    var kmers = 0
    var key = UInt64(0)
    var valid = 0
    var mask = (UInt64(1) << UInt64(2 * JOIN_PERIOD_K)) - 1
    for x in range(start, end):
        var code = encode_base(seq[x])
        valid = valid + 1 if code < 4 else 0
        key = ((key << 2) | UInt64(code & 3)) & mask
        if valid < JOIN_PERIOD_K:
            continue
        kmers += 1
        var seen = last.get(key)
        if seen and x - seen.value() <= JOIN_MAX_PERIOD:
            at_distance[x - seen.value()] += 1
        last[key] = x
    var best = 0
    var best_count = 0
    for p in range(2, min(JOIN_MAX_PERIOD, n // 2) + 1):
        var width = max(1, p // 20)
        var count = 0
        for d in range(max(2, p - width), min(JOIN_MAX_PERIOD, p + width) + 1):
            count += at_distance[d]
        if count > best_count:
            best = p
            best_count = count
    if Float64(best_count) < JOIN_PERIOD_SHARE * Float64(kmers):
        return 0
    return best


def _same_phase(right: List[UInt8], at: Int, near: Int, half: Int) -> Int:
    """The letter of right within `half` of `near` whose surroundings (`JOIN_PERIOD_K` either side) are most like
    those of right's letter `at`."""
    var best = near
    var best_same = -1
    for j in range(max(0, near - half), min(len(right), near + half + 1)):
        var same = 0
        for t in range(-JOIN_PERIOD_K, JOIN_PERIOD_K):
            var a = at + t
            var b = j + t
            if a >= 0 and b >= 0 and a < len(right) and b < len(right) and right[a] == right[b]:
                same += 1
        if same > best_same or (same == best_same and abs(j - near) < abs(best - near)):
            best = j
            best_same = same
    return best


def _versus(a: _Placing, b: _Placing) -> Tuple[Int, Int]:
    """Two copy hypotheses' edits on the left letters both align (from the later start of the two)."""
    var start = max(a.path_from, b.path_from)
    return (a.edits_from(start), b.edits_from(start))


def join_sides(
    mut aligner: Aligner,
    left: List[UInt8],
    right: List[UInt8],
    left_reach: List[Int],
    right_reach: List[Int],
    threads: Int = 1,
) raises -> List[UInt8]:
    """`left` and `right` (the two sides' consensus sequences) joined where `right`'s start overlaps `left`'s end, by
    `JOIN_LEAST` letters or more; empty when they do not overlap, or when where they do is not clear."""
    var trace = Trace()
    var found = pins(left, right)
    var unique = unique_pins(left, right, found)
    var keys = List[UInt64]()
    for m in range(len(found.r)):
        var start = found.r[m] - found.q[m]
        var end = start + len(right)
        if unique[m] and end >= len(left) and min(len(left), end) - max(0, start) >= JOIN_LEAST:
            keys.append((UInt64(start + (1 << 30)) << 32) | UInt64(m))
    sort(keys)
    var candidates = List[_Candidate]()
    var k = 0
    while k < len(keys):
        var group = List[Int]()
        var last = Int(keys[k] >> 32) - (1 << 30)
        while k < len(keys) and Int(keys[k] >> 32) - (1 << 30) - last <= JOIN_DRIFT:
            last = Int(keys[k] >> 32) - (1 << 30)
            group.append(Int(keys[k] & 0xFFFFFFFF))
            k += 1
        var best = -1
        var best_support = -1
        var best_off = 0
        for g in range(len(group)):
            var m = group[g]
            var support = min(_support(left_reach, found.r[m]), _support(right_reach, len(right) - found.q[m]))
            var off = abs(g - len(group) // 2)
            if support > best_support or (support == best_support and off < best_off):
                best = m
                best_support = support
                best_off = off
        candidates.append(
            _Candidate(
                found.r[best],
                found.q[best],
                String(
                    "the pin left ",
                    found.r[best],
                    " / right ",
                    found.q[best],
                    " (",
                    len(group),
                    " pins; ",
                    best_support,
                    " reads either side)",
                ),
            )
        )
    # The overlaps of the ends.
    var a_start = max(0, len(left) - JOIN_ENDS)
    var a = Span(left)[a_start:]
    var b = Span(right)[0 : min(len(right), JOIN_ENDS)]
    var mode = Mode.ends_free(reference_start=len(a), query_end=len(b)).with_match_score(1)
    for edit in [JOIN_EDIT, JOIN_EDIT_STRICT]:
        var got = aligner.align(a, b, Costs.linear(edit, edit), mode)
        if got.reference_end != len(a) or got.query_start != 0 or got.query_end - got.query_start < JOIN_LEAST:
            continue
        candidates.append(
            _Candidate(
                len(left),
                got.query_end,
                String("an overlap of the ends, ", got.query_end, " letters (edits ", edit, ")"),
            )
        )
    var placings = _eval_candidates(aligner, left, right, candidates, threads)
    var best = _best(placings)
    if best >= 0:
        var at = placings[best].diagonal
        var period = repeat_period(left, max(0, at), len(left))
        if period > 0:
            var cut_left = placings[best].cut_left
            var cut_right = placings[best].cut_right
            var from_how = placings[best].how.copy()
            var shift_candidates = List[_Candidate]()
            for shift in range(-JOIN_SHIFTS, JOIN_SHIFTS + 1):
                var meets = cut_right if shift == 0 else _same_phase(
                    right, cut_right, cut_right - shift * period, period // 2
                )
                var diagonal = cut_left - meets
                if meets < 0 or meets > len(right) or diagonal + len(right) < len(left):
                    continue
                if len(left) - max(0, diagonal) < JOIN_LEAST:
                    continue
                shift_candidates.append(
                    _Candidate(
                        cut_left,
                        meets,
                        String("copy shift ", shift, " (period ", period, ") of ", from_how),
                    )
                )
            var hypotheses = _eval_candidates(aligner, left, right, shift_candidates, threads)
            if len(hypotheses) > 0:
                return _copy_join(left, right, hypotheses, period, trace)
            best = _best(placings)
    if best < 0:
        trace.line(
            POA,
            String(
                "not joined: ", len(placings), " placings, none ", JOIN_IDENTITY, " alike over ", JOIN_LEAST, " letters"
            ),
        )
        return List[UInt8]()
    ref chosen = placings[best]
    for c in range(len(placings)):
        ref p = placings[c]
        var at_end = p.cut_left >= len(left)
        var meets = chosen.lands(p.cut_left - 1 if at_end else p.cut_left)
        var same = abs(meets - (p.cut_right - 1 if at_end else p.cut_right)) <= JOIN_DRIFT
        if not same and p.letters >= JOIN_LEAST:
            if p.rate() <= chosen.rate() * JOIN_MARGIN + JOIN_MARGIN_FLOOR:
                trace.line(
                    POA,
                    String(
                        "not joined: ",
                        chosen.how,
                        " at ",
                        chosen.rate(),
                        " edits a letter, ",
                        p.how,
                        " ",
                        p.diagonal - chosen.diagonal,
                        " letters off at ",
                        p.rate(),
                    ),
                )
                return List[UInt8]()
    trace.line(POA, String("join at ", chosen.how, ": ", chosen.edits, " edits in ", chosen.letters, " letters"))
    var out = List[UInt8](capacity=chosen.cut_left + len(right) - chosen.cut_right)
    for i in range(chosen.cut_left):
        out.append(left[i])
    for j in range(chosen.cut_right, len(right)):
        out.append(right[j])
    return out^
