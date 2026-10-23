"""Tandem-repeat arrays aligned copy by copy: a start for an MSA whose copies stack in phase.

A region's sequences that each hold one tandem array between shared flanks are cut into the array's copies, all
phased by one consensus unit, and the copies, each distinct one a symbol, are aligned across the sequences on their
sum of pairs: a copy present in some sequences and not others is gapped whole, and copy i of one sequence faces
copy i of another wherever that is cheapest. Letter-level alignment cannot see this: in a run of near-identical
copies, a stretch placed one copy early or late costs almost the same locally, and refinement one row at a time
settles in such places.

- The period: for every period from 2 up, the best stretch of letters matching the letter one period on; the
  period whose stretch scores best, or the smallest near-divisor of it scoring nearly as well (a period's multiples
  match as well, so the smallest explains the array).
- The copies: each sequence's array aligned to the unit repeated (a wrap-around DP at unit costs, the first copy
  starting and the last ending anywhere in the unit); a copy ends where the unit wraps, so every sequence's copies
  start at the same unit position.
- The copies' MSA: rows of symbols aligned progressively, longest first, each exactly to the profile of those
  before it (substituting copy a for copy b costs their edit distance, gapping a copy its length: the letters'
  sum of pairs at unit costs), then each realigned to the others until no row gains.
- Letters: each column of copies expanded through the unit, a column per unit position and the copies' insertions
  between two positions left-aligned.
`tandem_layout` gives the arrays' rows and where each array starts and ends; the flanks are the caller's.
"""

from karan import Aligner

comptime GAP = UInt8(ord("-"))
comptime MIN_PERIOD = 2
comptime MAX_PERIOD = 1000
comptime LOSS = 2
"""A letter differing from the one a period on costs this many matching letters in the period scan: arrays down to
two thirds identity score."""
comptime TOLERANT_LOSS = 4
"""The same with a letter matching one either side of a period on too: random letters match 58% of the time so,
and must score well below nothing."""
comptime MAX_SYMBOLS = 2000
"""Past this many distinct copies, no layout (the copies' distances cost the square)."""
comptime BIG = Int32(1 << 30)


@fieldwise_init
struct Repeat(Copyable, Movable):
    """A tandem array: letters [start, end) of a sequence repeat with `period`; `identity` of them match the letter
    one period on. Period 0: none."""

    var period: Int
    var start: Int
    var end: Int
    var identity: Float64


struct Layout(Movable):
    """The arrays of a region's sequences aligned copy by copy: `rows[k]` aligns letters [`starts[k]`, `ends[k]`) of
    sequence k; `period` is the unit's length. Empty `rows`: no layout (see `tandem_layout`)."""

    var period: Int
    var starts: List[Int]
    var ends: List[Int]
    var rows: List[List[UInt8]]

    def __init__(out self):
        self.period = 0
        self.starts = List[Int]()
        self.ends = List[Int]()
        self.rows = List[List[UInt8]]()


def _segment(seq: List[UInt8], p: Int, tolerant: Bool = False) -> Tuple[Int, Int, Int, Int]:
    """The best-scoring stretch of letters matching the letter `p` on (a match 1, a mismatch -LOSS; `tolerant`: the
    letter `p` - 1 or `p` + 1 on matches too, so an indel does not break the stretch for a whole copy, and a mismatch
    costs TOLERANT_LOSS): its score,
    its first letter, its end (past its last letter's partner) and its matches."""
    var n = len(seq)
    var loss = TOLERANT_LOSS if tolerant else LOSS
    var best = 0
    var cur = 0
    var cur_start = 0
    var cur_matches = 0
    var best_start = 0
    var best_end = 0
    var best_matches = 0
    for i in range(n - p):
        if (
            seq[i] == seq[i + p]
            or (tolerant and seq[i] == seq[i + p - 1])
            or (tolerant and i + p + 1 < n and seq[i] == seq[i + p + 1])
        ):
            cur += 1
            cur_matches += 1
        else:
            cur -= loss
        if cur <= 0:
            cur = 0
            cur_start = i + 1
            cur_matches = 0
        elif cur > best:
            best = cur
            best_start = cur_start
            best_end = i + 1
            best_matches = cur_matches
    if best == 0:
        return (0, 0, 0, 0)
    return (best, best_start, best_end + p, best_matches)


def _repeat_at(seq: List[UInt8], p: Int) -> Repeat:
    """The array of period `p`: letters matching one either side of a period on too (`_segment`'s `tolerant`)."""
    var seg = _segment(seq, p, True)
    if seg[0] == 0:
        return Repeat(0, 0, 0, 0.0)
    return Repeat(p, seg[1], seg[2], Float64(seg[3]) / Float64(max(1, seg[2] - seg[1] - p)))


def tandem_repeat(seq: List[UInt8], max_period: Int = MAX_PERIOD) -> Repeat:
    """The sequence's tandem array: for every period from 2 up to `max_period` (and a third of the sequence), the best
    stretch of letters matching a letter one period on, give or take one (`_segment`'s `tolerant`: an array whose
    copies carry indels still scores by the letters it covers); the period whose stretch scores best, or the smallest
    near-divisor of it (within a letter per copy) scoring at least 85% as well, since a period's multiples score
    nearly as well; then of it and its two neighbours, which tolerant matching scores alike, the one matching best
    exactly one period on."""
    var top = min(max_period, len(seq) // 3)
    if top < MIN_PERIOD:
        return Repeat(0, 0, 0, 0.0)
    var scores = List[Int](length=top + 1, fill=0)
    var best_p = 0
    for p in range(MIN_PERIOD, top + 1):
        scores[p] = _segment(seq, p, True)[0]
        if scores[p] > scores[best_p]:
            best_p = p
    if best_p == 0:
        return Repeat(0, 0, 0, 0.0)
    var chosen = best_p
    for p in range(MIN_PERIOD, best_p):
        var k = (best_p + p // 2) // p
        if k >= 2 and abs(best_p - k * p) <= k and scores[p] * 100 >= scores[best_p] * 85:
            chosen = p
            break
    var exact = chosen
    var exact_score = _segment(seq, chosen)[0]
    for p in [chosen - 1, chosen + 1]:
        if p >= MIN_PERIOD and p <= top:
            var score = _segment(seq, p)[0]
            if score > exact_score:
                exact = p
                exact_score = score
    return _repeat_at(seq, exact)


def _decompose(seq: List[UInt8], a: Int, b: Int, unit: List[UInt8]) -> Tuple[List[Int], List[Int32]]:
    """Letters [a, b) of `seq` aligned to `unit` repeated, at unit costs, the first copy starting and the last ending
    anywhere in the unit: each copy's first letter (offsets from a, rising) and each letter's place, its unit
    position or -(k + 1) for a letter inserted before position k (k = len(unit): after the last position, letters
    ending the last copy past its last unit letter). A copy ends where the unit wraps."""
    var m = b - a
    var p = len(unit)
    var cost = List[Int32](length=(m + 1) * p, fill=0)
    var move = List[UInt8](length=(m + 1) * p, fill=0)
    # State u: the next unit letter to face is unit[u]. Moves: 0 the letter faces unit[u - 1], 1 it is inserted,
    # 2 unit[u - 1] is deleted.
    for i in range(1, m + 1):
        var x = seq[a + i - 1]
        var row = i * p
        var prev = (i - 1) * p
        for u in range(p):
            var v = u - 1 if u > 0 else p - 1
            var diag = cost[prev + v] + (Int32(0) if unit[v] == x else Int32(1))
            var ins = cost[prev + u] + 1
            if diag <= ins:
                cost[row + u] = diag
                move[row + u] = 0
            else:
                cost[row + u] = ins
                move[row + u] = 1
        # Deletions: one sweep around the unit from its cheapest state, which no deletion improves.
        var start = 0
        for u in range(1, p):
            if cost[row + u] < cost[row + start]:
                start = u
        for k in range(1, p):
            var u = (start + k) % p
            var v = u - 1 if u > 0 else p - 1
            if cost[row + v] + 1 < cost[row + u]:
                cost[row + u] = cost[row + v] + 1
                move[row + u] = 2
    var u = 0
    for w in range(1, p):
        if cost[m * p + w] < cost[m * p + u]:
            u = w
    var places = List[Int32](length=m, fill=0)
    var cycle_of = List[Int](length=m, fill=0)
    var cycle = 0
    var last = p
    var i = m
    while i > 0:
        var mv = move[i * p + u]
        var v = u - 1 if u > 0 else p - 1
        if mv == 1:
            # Trailing letters past a last copy ending at the unit's end follow its last position.
            places[i - 1] = Int32(-(p + 1)) if last == p and u == 0 else Int32(-(u + 1))
            cycle_of[i - 1] = cycle
            i -= 1
            continue
        # A unit position above the last one taken (walking back) is the previous copy's.
        if v > last:
            cycle -= 1
        last = v
        if mv == 0:
            places[i - 1] = Int32(v)
            cycle_of[i - 1] = cycle
            i -= 1
        u = v
    var starts = List[Int]()
    for k in range(m):
        if k == 0 or cycle_of[k] != cycle_of[k - 1]:
            starts.append(k)
    return (starts^, places^)


def _consensus_unit(seq: List[UInt8], a: Int, b: Int, p: Int) -> List[UInt8]:
    """The array's unit: its first `p` letters, then twice the commonest letter facing each unit position."""
    var unit = List[UInt8](capacity=p)
    for k in range(p):
        unit.append(seq[a + k])
    for _ in range(2):
        var got = _decompose(seq, a, b, unit)
        var counts = List[Int](length=p * 256, fill=0)
        for k in range(b - a):
            var place = Int(got[1][k])
            if place >= 0:
                counts[place * 256 + Int(seq[a + k])] += 1
        for u in range(p):
            var best = Int(unit[u])
            for x in range(256):
                if counts[u * 256 + x] > counts[u * 256 + best]:
                    best = x
            unit[u] = UInt8(best)
    return unit^


def _column_costs(
    rows: List[List[Int]],
    skip: Int,
    dist: List[Int32],
    length: List[Int32],
    n_symbols: Int,
    mine: List[Int],
    weights: List[Int],
) -> Tuple[List[Int32], List[Int32], Int32]:
    """For a row `mine` (symbols) against the others' columns (row `skip` left out, columns no other row holds kept,
    harmlessly): each column's cost of each distinct symbol of `mine` facing it (the others' symbols' distances to
    it, their gaps its length), laid out column by column over `mine`'s distinct symbols in first-seen order; each
    column's cost of a gap facing it (the others' symbols' lengths); and what a symbol alone in a column costs per
    letter (the other rows). Every other row counts its weight."""
    var columns = len(rows[0]) if len(rows) > 0 else 0
    var others = 0
    for r in range(len(rows)):
        if r != skip:
            others += weights[r]
    var local = List[Int](length=n_symbols, fill=-1)
    var distinct = List[Int]()
    for x in mine:
        if local[x] < 0:
            local[x] = len(distinct)
            distinct.append(x)
    var nx = len(distinct)
    var facing = List[Int32](length=columns * nx, fill=0)
    var gap = List[Int32](length=columns, fill=0)
    var held = List[Int]()
    var counts = List[Int]()
    for c in range(columns):
        held.clear()
        counts.clear()
        var gaps = 0
        for r in range(len(rows)):
            if r == skip:
                continue
            var y = rows[r][c]
            if y < 0:
                gaps += weights[r]
                continue
            var found = False
            for t in range(len(held)):
                if held[t] == y:
                    counts[t] += weights[r]
                    found = True
                    break
            if not found:
                held.append(y)
                counts.append(weights[r])
        var g = Int32(0)
        for t in range(len(held)):
            g += Int32(counts[t]) * length[held[t]]
        gap[c] = g
        for xi in range(nx):
            var x = distinct[xi]
            var total = Int32(gaps) * length[x]
            for t in range(len(held)):
                total += Int32(counts[t]) * dist[x * n_symbols + held[t]]
            facing[c * nx + xi] = total
    return (facing^, gap^, Int32(others))


def _realign_symbols(
    mut rows: List[List[Int]],
    s: Int,
    mine: List[Int],
    dist: List[Int32],
    length: List[Int32],
    n_symbols: Int,
    weights: List[Int],
) -> Int32:
    """Row `s` (its symbols `mine`) realigned exactly to the others' columns: the rows rewritten (columns no row holds
    dropped). Returns the new row's cost against the others."""
    var costs = _column_costs(rows, s, dist, length, n_symbols, mine, weights)
    var local = List[Int](length=n_symbols, fill=-1)
    var nx = 0
    for x in mine:
        if local[x] < 0:
            local[x] = nx
            nx += 1
    var columns = len(rows[0]) if len(rows) > 0 else 0
    var n = len(mine)
    var width = n + 1
    var cost = List[Int32](length=(columns + 1) * width, fill=BIG)
    var move = List[UInt8](length=(columns + 1) * width, fill=0)
    cost[0] = 0
    for c in range(columns + 1):
        for j in range(n + 1):
            var here = c * width + j
            if c == 0 and j == 0:
                continue
            var best = BIG
            var how = UInt8(0)
            if c > 0 and j > 0:
                var v = cost[here - width - 1]
                if v < BIG:
                    v += costs[0][(c - 1) * nx + local[mine[j - 1]]]
                    if v < best:
                        best = v
                        how = 0
            if c > 0:
                var v = cost[here - width]
                if v < BIG:
                    v += costs[1][c - 1]
                    if v < best:
                        best = v
                        how = 1
            if j > 0:
                var v = cost[here - 1]
                if v < BIG:
                    v += costs[2] * length[mine[j - 1]]
                    if v < best:
                        best = v
                        how = 2
            cost[here] = best
            move[here] = how
    var steps = List[Int]()
    var letters = List[Int]()
    var c = columns
    var j = n
    while c > 0 or j > 0:
        var how = move[c * width + j]
        if how == 0:
            steps.append(c - 1)
            letters.append(mine[j - 1])
            c -= 1
            j -= 1
        elif how == 1:
            steps.append(c - 1)
            letters.append(-1)
            c -= 1
        else:
            steps.append(-1)
            letters.append(mine[j - 1])
            j -= 1
    var fresh = List[List[Int]]()
    for r in range(len(rows)):
        var row = List[Int](capacity=len(steps))
        for t in reversed(range(len(steps))):
            if r == s:
                row.append(letters[t])
            elif steps[t] >= 0:
                row.append(rows[r][steps[t]])
            else:
                row.append(-1)
        fresh.append(row^)
    var kept = List[Int]()
    for t in range(len(steps)):
        for r in range(len(fresh)):
            if fresh[r][t] >= 0:
                kept.append(t)
                break
    if len(kept) < len(steps):
        for r in range(len(fresh)):
            var row = List[Int](capacity=len(kept))
            for t in kept:
                row.append(fresh[r][t])
            fresh[r] = row^
    rows = fresh^
    return cost[columns * width + n]


def _row_cost(
    rows: List[List[Int]], s: Int, dist: List[Int32], length: List[Int32], n_symbols: Int, weights: List[Int]
) -> Int:
    """Row `s`'s cost against the others, as `_realign_symbols` counts it."""
    var total = 0
    for r in range(len(rows)):
        if r == s:
            continue
        var pair = 0
        for c in range(len(rows[s])):
            var x = rows[s][c]
            var y = rows[r][c]
            if x < 0 and y < 0:
                continue
            if x < 0:
                pair += Int(length[y])
            elif y < 0:
                pair += Int(length[x])
            else:
                pair += Int(dist[x * n_symbols + y])
        total += pair * weights[r]
    return total


def _symbol_msa(
    seqs: List[List[Int]], dist: List[Int32], length: List[Int32], n_symbols: Int, passes: Int, weights: List[Int]
) -> List[List[Int]]:
    """The symbol strings aligned on their sum of pairs, each pair counting the product of its `weights`:
    progressively, the longest first, then each row realigned to the others for at most `passes` passes, until none
    gains. Rows in `seqs`' order, -1 for gaps."""
    var n = len(seqs)
    var keys = List[UInt64]()
    for k in range(n):
        keys.append((UInt64((1 << 32) - len(seqs[k])) << 32) | UInt64(k))
    sort(keys)
    var order = List[Int]()
    for key in keys:
        order.append(Int(key & 0xFFFFFFFF))
    # rows[t] is sequence order[t]'s row, weighing w[t].
    var w = List[Int]()
    for t in range(n):
        w.append(weights[order[t]])
    var rows = List[List[Int]]()
    for t in range(n):
        var k = order[t]
        if t == 0:
            rows.append(seqs[k].copy())
            continue
        # The new row starts as all its symbols in columns of their own, then is placed exactly.
        var columns = len(rows[0])
        for r in range(len(rows)):
            for _ in range(len(seqs[k])):
                rows[r].append(-1)
        var mine = List[Int](length=columns, fill=-1)
        mine.extend(seqs[k].copy())
        rows.append(mine^)
        _ = _realign_symbols(rows, len(rows) - 1, seqs[k], dist, length, n_symbols, w)
    for _ in range(passes):
        var gained = 0
        for t in range(n):
            var before = _row_cost(rows, t, dist, length, n_symbols, w)
            var after = _realign_symbols(rows, t, seqs[order[t]], dist, length, n_symbols, w)
            gained += before - Int(after)
        if gained <= 0:
            break
    var out = List[List[Int]]()
    for _ in range(n):
        out.append(List[Int]())
    for t in range(n):
        out[order[t]] = rows[t].copy()
    return out^


def tandem_layout(seqs: List[List[UInt8]], passes: Int = 3, weights: List[Int] = List[Int]()) raises -> Layout:
    """The sequences' arrays aligned copy by copy (module docstring), or an empty layout where the sequences do not
    each hold one clear array of one period: the longest sequence's array (`tandem_repeat`) must match itself one
    period on at 75% or more, hold four copies and cover 40% of the sequence, and every sequence's array at that
    period must hold two copies at 70%. `passes` bounds the copies' MSA's refinement; a pair of sequences counts the
    product of their `weights` (1 each where empty)."""
    var out = Layout()
    var n = len(seqs)
    if n < 2:
        return out^
    var weight = weights.copy()
    if len(weight) != n:
        weight = List[Int](length=n, fill=1)
    var longest = 0
    for k in range(1, n):
        if len(seqs[k]) > len(seqs[longest]):
            longest = k
    var top = tandem_repeat(seqs[longest])
    var p = top.period
    if (
        p == 0
        or top.identity < 0.75
        or top.end - top.start < 4 * p
        or (top.end - top.start) * 10 < len(seqs[longest]) * 4
    ):
        return out^
    var arrays = List[Repeat]()
    for k in range(n):
        var r = _repeat_at(seqs[k], p)
        if r.period == 0 or r.end - r.start < 2 * p or r.identity < 0.7:
            return out^
        arrays.append(r^)
    var unit = _consensus_unit(seqs[longest], top.start, top.end, p)
    # Every sequence's copies, each distinct one a symbol; a symbol's letters' places from its first sighting.
    var symbol_of = Dict[String, Int]()
    var texts = List[List[UInt8]]()
    var places = List[List[Int32]]()
    var strings = List[List[Int]]()
    for k in range(n):
        var a = arrays[k].start
        var b = arrays[k].end
        var got = _decompose(seqs[k], a, b, unit)
        var mine = List[Int]()
        for c in range(len(got[0])):
            var lo = got[0][c]
            var hi = got[0][c + 1] if c + 1 < len(got[0]) else b - a
            var letters = List[UInt8](capacity=hi - lo)
            var placed_at = List[Int32](capacity=hi - lo)
            for t in range(lo, hi):
                letters.append(seqs[k][a + t])
                placed_at.append(got[1][t])
            var key = String(unsafe_from_utf8=letters.copy())
            if key not in symbol_of:
                symbol_of[key] = len(texts)
                texts.append(letters^)
                places.append(placed_at^)
            mine.append(symbol_of[key])
        strings.append(mine^)
    var n_symbols = len(texts)
    if n_symbols > MAX_SYMBOLS:
        return out^
    var aligner = Aligner()
    var length = List[Int32](capacity=n_symbols)
    for x in range(n_symbols):
        length.append(Int32(len(texts[x])))
    var dist = List[Int32](length=n_symbols * n_symbols, fill=0)
    for x in range(n_symbols):
        for y in range(x + 1, n_symbols):
            var d = Int32(aligner.distance(Span(texts[x]), Span(texts[y])))
            dist[x * n_symbols + y] = d
            dist[y * n_symbols + x] = d
    var msa = _symbol_msa(strings, dist, length, n_symbols, passes, weight)
    # Letters: each column of copies through the unit, insertions before each position left-aligned.
    var rows = List[List[UInt8]]()
    for _ in range(n):
        rows.append(List[UInt8]())
    # Slots 0 to p: insertions before each unit position, and after the last.
    var width = List[Int](length=p + 1, fill=0)
    var used = List[Bool](length=p, fill=False)
    var count = List[Int](length=p + 1, fill=0)
    for c in range(len(msa[0]) if n > 0 else 0):
        for u in range(p + 1):
            width[u] = 0
            if u < p:
                used[u] = False
        for k in range(n):
            var x = msa[k][c]
            if x < 0:
                continue
            for u in range(p + 1):
                count[u] = 0
            for place in places[x]:
                if place >= 0:
                    used[Int(place)] = True
                else:
                    count[Int(-place) - 1] += 1
            for u in range(p + 1):
                width[u] = max(width[u], count[u])
        # Block columns: for each position u, its insertions' columns, then its own (where any copy holds it); then
        # the insertions after the last position.
        var offset = List[Int](length=p + 1, fill=0)
        var block = 0
        for u in range(p + 1):
            offset[u] = block
            block += width[u] + (1 if u < p and used[u] else 0)
        for k in range(n):
            var cells = List[UInt8](length=block, fill=GAP)
            var x = msa[k][c]
            if x >= 0:
                for u in range(p + 1):
                    count[u] = 0
                for t in range(len(places[x])):
                    var place = Int(places[x][t])
                    if place >= 0:
                        cells[offset[place] + width[place]] = texts[x][t]
                    else:
                        var u = -place - 1
                        cells[offset[u] + count[u]] = texts[x][t]
                        count[u] += 1
            rows[k].extend(cells^)
    out.period = p
    for k in range(n):
        out.starts.append(arrays[k].start)
        out.ends.append(arrays[k].end)
    out.rows = rows^
    return out^
