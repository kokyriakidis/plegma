"""minimap2's pinned alignment of a sequence to the reference, so `on_reference` writes each sequence's differences
against the reference as minimap2 (2.31, -x asm5 as dipcall runs it: -z200000,10000 -r2k), and so dipcall and the
benchmarks dipcall makes, write them.

Sketching and reference indexing live in `seeds.mojo`; RMQ chaining in `chain.mojo`; this module applies minimap2's
bad-end and bad-seed filters, chooses pins, aligns the stretches between them exactly (karan), and left-shifts gaps
as `mm_fix_cigar` does.
"""

from std.atomic import Atomic

from karan import Aligner, Costs, Mode, Ties
from karan.common import spread, thread_count

from .chain import BW, MAX_GAP, MIN_CHAIN_SCORE, primary_chain
from .cigar import append_pairing_from_cigar
from .seeds import IGNORE, K, LONG_JOIN, POSITION_BITS, ReferenceIndex, Seed, Sketch, TANDEM, sketch
from .trace import PINS, Trace

comptime MIN_KSW_LEN = 200
comptime PIN_MEMORY = 2 << 30


struct Pins(Movable):
    """The points (reference letter, query letter) an alignment passes through, in order."""

    var r: List[Int]
    var q: List[Int]

    def __init__(out self):
        self.r = List[Int]()
        self.q = List[Int]()


def _fix_bad_ends(a: List[Seed], bw: Int, min_match: Int) -> Tuple[Int, Int]:
    """mm_fix_bad_ends: the first seed kept and the seeds' count."""
    var cnt = len(a)
    if cnt < 3:
        return (0, cnt)
    var mlen = a[0].span
    for i in range(1, cnt):
        var tl = a[i].x - a[i - 1].x
        var ql = a[i].q - a[i - 1].q
        mlen += a[i].span if tl > a[i].span and ql > a[i].span else min(tl, ql)
    var first = 0
    var m = a[0].span
    var l = m
    for i in range(1, cnt - 1):
        if a[i].flags & LONG_JOIN:
            break
        var lr = a[i].x - a[i - 1].x
        var lq = a[i].q - a[i - 1].q
        var mn = min(lr, lq)
        if max(lr, lq) - mn > l >> 1:
            first = i
        l += mn
        m += min(mn, a[i].span)
        if l >= bw << 1 or (m >= min_match and m >= bw) or m >= mlen >> 1:
            break
    var count = cnt - first
    m = a[cnt - 1].span
    l = m
    var i = cnt - 2
    while i > first:
        if a[i + 1].flags & LONG_JOIN:
            break
        var lr = a[i + 1].x - a[i].x
        var lq = a[i + 1].q - a[i].q
        var mn = min(lr, lq)
        if max(lr, lq) - mn > l >> 1:
            count = i + 1 - first
        l += mn
        m += min(mn, a[i + 1].span)
        if l >= bw << 1 or (m >= min_match and m >= bw) or m >= mlen >> 1:
            break
        i -= 1
    return (first, count)


def _gap(a: List[Seed], i: Int) -> Int:
    return (a[i].q - a[i - 1].q) - (a[i].x - a[i - 1].x)


def _long_gaps(a: List[Seed], min_gap: Int) -> List[Int]:
    var out = List[Int]()
    for i in range(1, len(a)):
        var g = _gap(a, i)
        if g > min_gap or g < -min_gap:
            out.append(i)
    if len(out) <= 1:
        return List[Int]()
    return out^


def _filter_bad_seeds(mut a: List[Seed]):
    """mm_filter_bad_seeds(min_gap 10, diff_thres 40, max_ext_len MAX_GAP / 2, max_ext_cnt 10): the seeds of a
    stretch whose long insertions and deletions mostly cancel are not pinned."""
    var ks = _long_gaps(a, 10)
    var n = len(ks)
    if n == 0:
        return
    var best = 0
    var max_st = -1
    var max_en = -1
    var k = 0
    while True:
        if k == n or k >= max_en:
            if max_en > 0:
                for i in range(ks[max_st], ks[max_en]):
                    a[i].flags |= IGNORE
            best = 0
            max_st = -1
            max_en = -1
            if k == n:
                break
        var i = ks[k]
        var g = _gap(a, i)
        var n_ins = max(g, 0)
        var n_del = max(-g, 0)
        var qs = a[i - 1].q
        var rs = a[i - 1].x
        var max_diff = 0
        var max_diff_l = -1
        for l in range(k + 1, min(n, k + 11)):
            var j = ks[l]
            if a[j].q - qs > MAX_GAP >> 1 or a[j].x - rs > MAX_GAP >> 1:
                break
            g = _gap(a, j)
            if g > 0:
                n_ins += g
            else:
                n_del -= g
            var diff = n_ins + n_del - abs(n_ins - n_del)
            if max_diff < diff:
                max_diff = diff
                max_diff_l = l
        if max_diff > 40 and max_diff > best:
            best = max_diff
            max_st = k
            max_en = max_diff_l
        k += 1


def _filter_bad_seeds_alt(mut a: List[Seed]):
    """mm_filter_bad_seeds_alt(min_gap 30, max_ext MAX_GAP / 2): long gaps close together joined into one."""
    var ks = _long_gaps(a, 30)
    var n = len(ks)
    var k = 0
    while k < n:
        var i = ks[k]
        var gap1 = abs(_gap(a, i))
        var re1 = a[i].x
        var qe1 = a[i].q
        var l = k + 1
        while l < n:
            var j = ks[l]
            if a[j].q - qe1 > MAX_GAP >> 1 or a[j].x - re1 > MAX_GAP >> 1:
                break
            var gap2 = abs(_gap(a, j))
            var pre = a[j - 1].span
            if min(a[j - 1].x + pre - re1, a[j - 1].q + pre - qe1) > gap1 + gap2:
                break
            re1 = a[j].x
            qe1 = a[j].q
            gap1 = gap2
            l += 1
        if l > k + 1:
            var end = ks[l - 1]
            for j in range(ks[k], end):
                a[j].flags |= IGNORE
            a[end].flags |= LONG_JOIN
        k = l


def _pins_from_chain(chain: List[Seed], ref_len: Int, query_len: Int) -> Pins:
    """Selects pins from a primary chain after bad-end and bad-seed filtering."""
    var out = Pins()
    var trace = Trace()
    if trace.on(PINS):
        trace.line(PINS, String("chain of ", len(chain), " seeds, ", ref_len, " x ", query_len, " letters"))
    if len(chain) == 0:
        return out^
    var kept = _fix_bad_ends(chain, BW, MIN_CHAIN_SCORE * 2)
    var a = List[Seed](capacity=kept[1])
    for i in range(kept[0], kept[0] + kept[1]):
        a.append(chain[i])
    _filter_bad_seeds(a)
    _filter_bad_seeds_alt(a)
    var half = K >> 1
    var rs = a[0].x - half
    var qs = a[0].q - half
    var points_r = List[Int]()
    var points_q = List[Int]()
    points_r.append(rs)
    points_q.append(qs)
    for i in range(1, len(a)):
        var last = i == len(a) - 1
        if a[i].flags & (IGNORE | TANDEM) and not last:
            continue
        var re = a[i].x - half
        var qe = a[i].q - half
        if last or a[i].flags & LONG_JOIN or (qe - qs >= MIN_KSW_LEN and re - rs >= MIN_KSW_LEN):
            points_r.append(re)
            points_q.append(qe)
            rs = re
            qs = qe
    for m in range(len(points_r)):
        if 0 < points_r[m] < ref_len and 0 < points_q[m] < query_len:
            out.r.append(points_r[m])
            out.q.append(points_q[m])
    return out^


def pins_with_index(index: ReferenceIndex, query: List[UInt8]) -> Pins:
    """`pins` against a pre-built `ReferenceIndex`."""
    var trace = Trace()
    if index.length >= 1 << POSITION_BITS or len(query) >= 1 << POSITION_BITS:
        trace.line(PINS, String("no pins: ", index.length, " x ", len(query), " letters, past 2^", POSITION_BITS))
        return Pins()
    return _pins_from_chain(primary_chain(index, query), index.length, len(query))


def pins(reference: List[UInt8], query: List[UInt8]) -> Pins:
    """mm_align1's pins for `query` against `reference`: the middles of the primary chain's kept seeds where both
    sequences advanced MIN_KSW_LEN letters or more since the last pin (or a long join), the last seed's always."""
    var trace = Trace()
    if len(reference) >= 1 << POSITION_BITS or len(query) >= 1 << POSITION_BITS:
        trace.line(PINS, String("no pins: ", len(reference), " x ", len(query), " letters, past 2^", POSITION_BITS))
        return Pins()
    return _pins_from_chain(primary_chain(reference, query), len(reference), len(query))


def _concat3(a: List[UInt8], b: List[UInt8], c: List[UInt8]) -> List[UInt8]:
    """Concatenates `a + b + c` into one buffer with a single allocation."""
    var out = List[UInt8](capacity=len(a) + len(b) + len(c))
    for x in a:
        out.append(x)
    for x in b:
        out.append(x)
    for x in c:
        out.append(x)
    return out^


def _crop_context_pins(found: Pins, left_len: Int, ref_len: Int, query_len: Int) -> Pins:
    var out = Pins()
    for m in range(len(found.r)):
        var r = found.r[m] - left_len
        var q = found.q[m] - left_len
        if 0 < r < ref_len and 0 < q < query_len:
            if len(out.r) == 0 or (r >= out.r[len(out.r) - 1] and q >= out.q[len(out.q) - 1]):
                out.r.append(r)
                out.q.append(q)
    return out^


def pins_in_context(reference: List[UInt8], query: List[UInt8], left: List[UInt8], right: List[UInt8]) -> Pins:
    """The pins of `left + reference + right` against `left + query + right` (the flanks context for the seeds and the
    chain only) that fall inside `reference` and `query`, in their own coordinates, each past the last on both."""
    var whole_ref = _concat3(left, reference, right)
    var whole_query = _concat3(left, query, right)
    return _crop_context_pins(pins(whole_ref, whole_query), len(left), len(reference), len(query))


def pins_in_context_with_index(
    ref_index: ReferenceIndex, ref_len: Int, query: List[UInt8], left: List[UInt8], right: List[UInt8]
) -> Pins:
    """`pins_in_context` using a pre-built `ReferenceIndex` of `left + reference + right`."""
    var whole_query = _concat3(left, query, right)
    return _crop_context_pins(pins_with_index(ref_index, whole_query), len(left), ref_len, len(query))


def _align_pinned_windows(
    mut aligner: Aligner,
    reference: List[UInt8],
    query: List[UInt8],
    found: Pins,
    right_flank: Bool,
    costs: Costs,
    threads: Int = 1,
) raises -> List[Int]:
    var trace = Trace()
    if trace.on(PINS):
        var text = String(len(found.r), " pins (reference,query):")
        for m in range(len(found.r)):
            text.write(" ", found.r[m], ",", found.q[m])
        trace.line(PINS, text)
    var points_r = List[Int](capacity=len(found.r) + 2)
    var points_q = List[Int](capacity=len(found.q) + 2)
    points_r.append(0)
    points_q.append(0)
    for x in found.r:
        points_r.append(x)
    for x in found.q:
        points_q.append(x)
    points_r.append(len(reference))
    points_q.append(len(query))
    var num_windows = len(points_r) - 1
    var workers = thread_count(threads, num_windows)
    var out = List[Int](capacity=len(query))
    if workers > 1 and num_windows >= 2 and len(query) >= 1024:
        var win_order = List[UInt64](capacity=num_windows)
        for m in range(num_windows):
            var span_m = max(points_r[m + 1] - points_r[m], points_q[m + 1] - points_q[m])
            win_order.append((UInt64(span_m) << 32) | UInt64(m))
        sort(win_order)
        win_order.reverse()
        var parts = List[List[Int]](length=num_windows, fill=List[Int]())
        var failed = List[Bool](length=num_windows, fill=False)
        var parts_p = parts.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def run_windows(
            slot: Int,
        ) {
            mut taken,
            imm win_order,
            imm points_r,
            imm points_q,
            imm reference,
            imm query,
            imm costs,
            imm parts_p,
            imm failed_p,
        }:
            var local_aligner = Aligner(warm=False)
            while True:
                var pos = Int(taken.fetch_add(1))
                if pos >= len(win_order):
                    return
                var m = Int(win_order[pos] & 0xFFFFFFFF)
                var r0 = points_r[m]
                var r1 = points_r[m + 1]
                var q0 = points_q[m]
                var q1 = points_q[m + 1]
                if q1 == q0:
                    continue
                var local_out = List[Int](capacity=q1 - q0)
                if r1 == r0:
                    for _ in range(q0, q1):
                        local_out.append(-1)
                    parts_p[unsafe_offset=m] = local_out^
                    continue
                try:
                    var aligned = local_aligner.align(
                        Span(reference)[r0:r1],
                        Span(query)[q0:q1],
                        costs,
                        Mode.GLOBAL,
                        ties=Ties.LEFT,
                        max_memory=PIN_MEMORY,
                    )
                    append_pairing_from_cigar(aligned.cigar, r0, local_out)
                    parts_p[unsafe_offset=m] = local_out^
                except:
                    failed_p[unsafe_offset=m] = True

        spread(run_windows, workers, workers)
        for m in range(num_windows):
            if failed[m]:
                var r0 = points_r[m]
                var r1 = points_r[m + 1]
                var q0 = points_q[m]
                var q1 = points_q[m + 1]
                var aligned = aligner.align(
                    Span(reference)[r0:r1],
                    Span(query)[q0:q1],
                    costs,
                    Mode.GLOBAL,
                    ties=Ties.LEFT,
                    max_memory=PIN_MEMORY,
                )
                append_pairing_from_cigar(aligned.cigar, r0, out)
            else:
                for p in parts[m]:
                    out.append(p)
        return _left_aligned(out, reference, query, right_flank)

    for m in range(num_windows):
        var r0 = points_r[m]
        var r1 = points_r[m + 1]
        var q0 = points_q[m]
        var q1 = points_q[m + 1]
        if q1 == q0:
            continue
        if r1 == r0:
            for _ in range(q0, q1):
                out.append(-1)
            continue
        var aligned = aligner.align(
            Span(reference)[r0:r1], Span(query)[q0:q1], costs, Mode.GLOBAL, ties=Ties.LEFT, max_memory=PIN_MEMORY
        )
        append_pairing_from_cigar(aligned.cigar, r0, out)
    return _left_aligned(out, reference, query, right_flank)


def pinned_pairing(
    mut aligner: Aligner,
    reference: List[UInt8],
    query: List[UInt8],
    left: List[UInt8],
    right: List[UInt8],
    costs: Costs,
    threads: Int = 1,
) raises -> List[Int]:
    """For each letter of `query`, the letter of `reference` its pinned alignment pairs it with (-1: inserted): the
    stretches between the pins (`pins_in_context`) each aligned exactly at `costs`, gaps left. Raises past 2^25
    letters (flanks and sequence), where the seeds' positions no longer fit: no pins would leave one exact alignment
    of the whole, far past any memory."""
    var longest = len(left) + max(len(reference), len(query)) + len(right)
    if longest >= 1 << POSITION_BITS:
        raise Error(t"pinned alignment: {longest} letters with the flanks, past the 2^{POSITION_BITS} pins reach")
    var found = pins_in_context(reference, query, left, right)
    return _align_pinned_windows(aligner, reference, query, found, len(right) > 0, costs, threads)


def pinned_pairing_with_index(
    mut aligner: Aligner,
    ref_index: ReferenceIndex,
    reference: List[UInt8],
    query: List[UInt8],
    left: List[UInt8],
    right: List[UInt8],
    costs: Costs,
    threads: Int = 1,
) raises -> List[Int]:
    """`pinned_pairing` reusing a pre-built `ReferenceIndex` of `left + reference + right`."""
    var longest = len(left) + max(len(reference), len(query)) + len(right)
    if longest >= 1 << POSITION_BITS:
        raise Error(t"pinned alignment: {longest} letters with the flanks, past the 2^{POSITION_BITS} pins reach")
    var found = pins_in_context_with_index(ref_index, len(reference), query, left, right)
    return _align_pinned_windows(aligner, reference, query, found, len(right) > 0, costs, threads)


comptime _M = 0
comptime _I = 1
comptime _D = 2


def _left_aligned(pairing: List[Int], reference: List[UInt8], query: List[UInt8], right_flank: Bool) -> List[Int]:
    """`pairing` with every insertion and deletion shifted left through the matched run before it while the letters
    allow, as minimap2's mm_fix_cigar does to an alignment last (its first loop, one pass, each indel between two
    match runs; the windows' own alignments leave a gap left only within the window, so a gap minimap2 writes left
    of a pin was written right of it). The right flank (`right_flank`) stands for a matched run after the sequence,
    so an indel at its end may shift too; none shifts into the left flank, which no row holds (a gap minimap2 shifts
    into it is the same letters gapped at the sequence's start)."""
    var ops = List[Int]()
    var counts = List[Int]()

    def push(op: Int, n: Int) {mut ops, mut counts}:
        if n <= 0:
            return
        if len(ops) > 0 and ops[len(ops) - 1] == op:
            counts[len(counts) - 1] += n
        else:
            ops.append(op)
            counts.append(n)

    var next_ref = 0
    for j in range(len(pairing)):
        var p = pairing[j]
        if p < 0:
            push(_I, 1)
            continue
        push(_D, p - next_ref)
        push(_M, 1)
        next_ref = p + 1
    push(_D, len(reference) - next_ref)
    var toff = 0
    var qoff = 0
    for k in range(len(ops)):
        var op = ops[k]
        var n = counts[k]
        if op == _M:
            toff += n
            qoff += n
            continue
        var after = (k + 1 < len(ops) and ops[k + 1] == _M) or (k + 1 == len(ops) and right_flank)
        if k > 0 and ops[k - 1] == _M and after:
            var before = counts[k - 1]
            var l = 0
            while l < before:
                if op == _I and query[qoff - 1 - l] != query[qoff + n - 1 - l]:
                    break
                if op == _D and reference[toff - 1 - l] != reference[toff + n - 1 - l]:
                    break
                l += 1
            if l > 0:
                counts[k - 1] -= l
                if k + 1 < len(ops):
                    counts[k + 1] += l
                else:
                    ops.append(_M)
                    counts.append(l)
                qoff -= l
                toff -= l
        if op == _I:
            qoff += n
        else:
            toff += n
    var out = List[Int](capacity=len(pairing))
    var i = 0
    for k in range(len(ops)):
        for _ in range(counts[k]):
            if ops[k] == _M:
                out.append(i)
                i += 1
            elif ops[k] == _I:
                out.append(-1)
            else:
                i += 1
    return out^
