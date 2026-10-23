"""An MSA refined on its sum of pairs: each distinct sequence realigned in turn to the profile of the others.

The sum of pairs at unit costs, over every pair of distinct sequences, counts the columns where the two rows
differ, a column where both are gaps costing nothing: what vntr-eval's all_cost_over_opt divides by the pairs'
optimal costs. With the other rows fixed, a sequence's best row is exact by a DP over the columns: its letter
against column c costs the others that do not hold that letter there (their gaps included), its gap there costs
the others that hold a letter, and its letter in a column of its own costs every other row. The DP keeps within
`band` of the sequence's present place, which is one of its choices, so no pass raises the sum. Passes go on until
one gains nothing. On a hotspot VNTR whose POA (abPOA's, or plegma's longest first) cost 1.29 times the optima,
three passes took it to 1.13.
"""

from karan import Aligner, Costs, Mode, Ties

from .align import anchored_alignment
from .backbone import Level, MUMS
from .chain import ChainCosts
from .graph import DELETION, INSERTION, cigar_runs


comptime GAP = UInt8(ord("-"))
comptime SYMBOLS = 6
"""A, C, G, T, any other letter, and the gap."""
comptime BIG = Int32(1 << 29)
comptime MATCHED = UInt8(0)
comptime GAPPED = UInt8(1)
comptime OWN = UInt8(2)
comptime EDGE = 4


@always_inline
def _symbol(letter: UInt8) -> Int:
    var lower = letter | 0x20
    if letter == GAP:
        return 5
    if lower == UInt8(ord("a")):
        return 0
    if lower == UInt8(ord("c")):
        return 1
    if lower == UInt8(ord("g")):
        return 2
    if lower == UInt8(ord("t")):
        return 3
    return 4


def gap_opens(rows: List[List[UInt8]]) -> Int:
    """The gaps the rows' pairs open, summed over every pair: in each pair's alignment (columns where both are gaps
    left out), the runs of columns where just one of the two has a gap."""
    var total = 0
    var n = len(rows)
    for a in range(n):
        for b in range(a + 1, n):
            var state = 0
            for c in range(len(rows[a])):
                var x = rows[a][c] == GAP
                var y = rows[b][c] == GAP
                if x and y:
                    continue
                var now = 1 if x else (2 if y else 0)
                if now != 0 and now != state:
                    total += 1
                state = now
    return total


def sum_of_pairs(rows: List[List[UInt8]]) -> Int:
    """The rows' sum of pairs at unit costs (see the module docstring). Letters other than ACGT count as one."""
    if len(rows) == 0:
        return 0
    var n = len(rows)
    var total = 0
    var counts = List[Int](length=SYMBOLS, fill=0)
    for c in range(len(rows[0])):
        for x in range(SYMBOLS):
            counts[x] = 0
        for r in range(n):
            counts[_symbol(rows[r][c])] += 1
        var same = 0
        for x in range(SYMBOLS):
            same += counts[x] * (counts[x] - 1) // 2
        total += n * (n - 1) // 2 - same
    return total


comptime PARTNERS = 24
"""Past this many other distinct sequences, each sequence is aligned to its nearest this many alone (by unit edit
distance), as are those it is nearest to: the anchors' votes need a few dozen, and all pairs cost the square."""


def _pair_runs(
    mut aligner: Aligner, partner: ImmSpan[UInt8, _], query: ImmSpan[UInt8, _], costs: Costs, seeded: Bool
) raises -> Tuple[List[UInt8], List[Int]]:
    """The runs of `query`'s alignment to `partner`: optimal under `costs`, or with `seeded` through anchors."""
    if seeded:
        var levels = List[Level]()
        for k in [51, 31, 21]:
            levels.append(Level(MUMS, k))
        var through = anchored_alignment(
            aligner,
            partner,
            query,
            levels,
            ChainCosts(substitution=0.1),
            costs,
            0,
            True,
            repeat_k=12,
            repeat_distance=500,
            mask_depth=1,
        )
        return (through.ops.copy(), through.counts.copy())
    var found = aligner.align(partner, query, costs, Mode.GLOBAL, ties=Ties.LEFT)
    var runs = cigar_runs(found.cigar)
    return (runs.ops.copy(), runs.counts.copy())


struct Partners(Movable):
    """Each distinct sequence's alignment to its partners (every other, or past PARTNERS others its nearest), as
    the letter of the partner each of its letters faces: `facing[s][r][j]` (-1 where letter j is inserted; the list
    empty where r is no partner of s). They do not depend on the MSA, so they are made once."""

    var facing: List[List[List[Int32]]]

    def __init__(out self, spelled: List[List[UInt8]], costs: Costs = Costs.edit(), seeded: Bool = False) raises:
        """Optimal alignments under `costs` (karan, indels placed left); with `seeded`, through anchors instead
        (`anchored_alignment`): MUMs of 51, 31 and 21 letters outside tandem repeats pin the flanks, then inside each
        gap between them, a repeat array, MUMs unique within the gap's two pieces, which sit on copies that carry
        something of their own and so pair homologous copies; only the stretches between anchors are aligned for
        least cost."""
        var n = len(spelled)
        self.facing = List[List[List[Int32]]]()
        for _ in range(n):
            var row = List[List[Int32]]()
            for _ in range(n):
                row.append(List[Int32]())
            self.facing.append(row^)
        var aligner = Aligner()
        # Which pairs to align: all, or past PARTNERS each sequence's nearest.
        var wanted = List[Bool](length=n * n, fill=n - 1 <= PARTNERS)
        if n - 1 > PARTNERS:
            var distance = List[Int](length=n * n, fill=0)
            for a in range(n):
                for b in range(a + 1, n):
                    var d = aligner.distance(Span(spelled[a]), Span(spelled[b]))
                    distance[a * n + b] = d
                    distance[b * n + a] = d
            var keys = List[UInt64]()
            for a in range(n):
                keys.clear()
                for b in range(n):
                    if b != a:
                        keys.append((UInt64(distance[a * n + b]) << 32) | UInt64(b))
                sort(keys)
                for t in range(PARTNERS):
                    var b = Int(keys[t] & 0xFFFFFFFF)
                    wanted[a * n + b] = True
                    wanted[b * n + a] = True
        for a in range(n):
            for b in range(a + 1, n):
                if not wanted[a * n + b] or len(spelled[a]) == 0 or len(spelled[b]) == 0:
                    continue
                var ab = List[Int32](length=len(spelled[a]), fill=-1)
                var ba = List[Int32](length=len(spelled[b]), fill=-1)
                var runs = _pair_runs(aligner, Span(spelled[b]), Span(spelled[a]), costs, seeded)
                var i = 0
                var j = 0
                for r in range(len(runs[0])):
                    var op = runs[0][r]
                    for _ in range(runs[1][r]):
                        if op == DELETION:
                            i += 1
                        elif op == INSERTION:
                            j += 1
                        else:
                            ab[j] = Int32(i)
                            ba[i] = Int32(j)
                            i += 1
                            j += 1
                self.facing[a][b] = ab^
                self.facing[b][a] = ba^


def _anchored_centre(rows: List[List[UInt8]], s: Int, partners: Partners, mut centre: List[Int]):
    """Where each column's band is centred for row `s`: the sequence's letters before the column as its confident
    anchors put them. Each other row votes, through its optimal pairwise alignment to the sequence, a column for
    every letter of the sequence; a letter is anchored where more than half of the rows vote one column, and the
    anchors kept are the most that rise in both letter and column (a longest increasing chain). Between them the
    centre is interpolated."""
    var columns = len(rows[0])
    var n_rows = len(rows)
    # Each row's letters' columns.
    var column_of = List[List[Int32]]()
    var n = 0
    for r in range(n_rows):
        var cols = List[Int32]()
        for c in range(columns):
            if rows[r][c] != GAP:
                cols.append(Int32(c))
        if r == s:
            n = len(cols)
        column_of.append(cols^)
    var anchor_j = List[Int]()
    var anchor_c = List[Int]()
    var votes = List[Int32]()
    var voters = 0
    for r in range(n_rows):
        if r != s and len(partners.facing[s][r]) > 0:
            voters += 1
    for j in range(n):
        votes.clear()
        for r in range(n_rows):
            if r == s or len(partners.facing[s][r]) == 0:
                continue
            var k = Int(partners.facing[s][r][j])
            if k >= 0:
                votes.append(column_of[r][k])
        if len(votes) * 2 <= voters:
            continue
        sort(votes)
        var best = -1
        var best_count = 0
        var run = 0
        for t in range(len(votes)):
            run = run + 1 if t > 0 and votes[t] == votes[t - 1] else 1
            if run > best_count:
                best_count = run
                best = Int(votes[t])
        if best_count * 2 > voters:
            anchor_j.append(j)
            anchor_c.append(best)
    # The longest chain rising in column (letters already rise): patience sorting.
    var tails = List[Int]()
    var back = List[Int](length=len(anchor_c), fill=-1)
    for t in range(len(anchor_c)):
        var lo = 0
        var hi = len(tails)
        while lo < hi:
            var mid = (lo + hi) // 2
            if anchor_c[tails[mid]] < anchor_c[t]:
                lo = mid + 1
            else:
                hi = mid
        if lo > 0:
            back[t] = tails[lo - 1]
        if lo == len(tails):
            tails.append(t)
        else:
            tails[lo] = t
    var chain = List[Int]()
    var t = tails[len(tails) - 1] if len(tails) > 0 else -1
    while t >= 0:
        chain.append(t)
        t = back[t]
    # Points (column, letters before it): the start, each anchor's column and the one after it, the end.
    var point_c = List[Int]()
    var point_j = List[Int]()
    point_c.append(0)
    point_j.append(0)
    for q in reversed(range(len(chain))):
        var a = chain[q]
        if anchor_c[a] > point_c[len(point_c) - 1] and anchor_j[a] >= point_j[len(point_j) - 1]:
            point_c.append(anchor_c[a])
            point_j.append(anchor_j[a])
            point_c.append(anchor_c[a] + 1)
            point_j.append(anchor_j[a] + 1)
    if point_c[len(point_c) - 1] < columns:
        point_c.append(columns)
        point_j.append(n)
    centre.resize(columns + 1, 0)
    var p = 0
    for c in range(columns + 1):
        while p + 1 < len(point_c) - 1 and point_c[p + 1] <= c:
            p += 1
        var c0 = point_c[p]
        var c1 = point_c[min(p + 1, len(point_c) - 1)]
        var j0 = point_j[p]
        var j1 = point_j[min(p + 1, len(point_j) - 1)]
        if c1 <= c0:
            centre[c] = min(n, max(0, j0))
        else:
            centre[c] = min(n, max(0, j0 + (c - c0) * (j1 - j0) // (c1 - c0)))


def _realign(
    mut rows: List[List[UInt8]],
    s: Int,
    band: Int,
    max_band: Int,
    centre: List[Int],
    gap_open: Int32,
    mut cost: List[Int32],
    mut moves: List[UInt8],
):
    """Row `s` realigned to the others' profile within `band` letters of its place (or of `centre`, when given),
    widened while the best path touches the band's edge, up to `max_band`; the rows rewritten (all-gap columns
    dropped), unless the new row would cost more than the present one.

    Three states per cell: the sequence's letter in a column (M), its gap there (G), its letter in a column of its
    own (O). Besides the unit costs, every gap a pair newly opens costs `gap_open`, counted from the columns' profile:
    the sequence's gap opening against the rows holding a letter there, a row's own gap starting against the
    sequence's letter, the sequence's own column against every row holding a letter beside it."""
    var columns = len(rows[0])
    var others = Int32(len(rows) - 1)
    # The others' symbol counts per column, and per column how many of them hold a letter, start a gap (a letter
    # in the column before), or end one (a gap before).
    var counts = List[Int32](length=columns * SYMBOLS, fill=0)
    var held = List[Int32](length=columns, fill=0)
    var starts = List[Int32](length=columns, fill=0)
    var ends = List[Int32](length=columns, fill=0)
    for r in range(len(rows)):
        if r == s:
            continue
        var was_letter = True
        for c in range(columns):
            var x = rows[r][c]
            counts[c * SYMBOLS + _symbol(x)] += 1
            if x != GAP:
                held[c] += 1
                if not was_letter:
                    ends[c] += 1
            elif was_letter:
                starts[c] += 1
            was_letter = x != GAP
    var letters = List[UInt8]()
    var before = List[Int](length=columns + 1, fill=0)
    for c in range(columns):
        before[c] = len(letters)
        if rows[s][c] != GAP:
            letters.append(rows[s][c])
    before[columns] = len(letters)
    var n = len(letters)
    var mid = centre.copy() if len(centre) == columns + 1 else before.copy()
    var g = gap_open
    # The present row's cost under the same model, which the new one must not pass.
    var present = Int32(0)
    var in_gap = False
    for c in range(columns):
        if rows[s][c] != GAP:
            present += others - counts[c * SYMBOLS + _symbol(rows[s][c])] + g * starts[c]
            in_gap = False
        else:
            present += held[c] + (g * ends[c] if in_gap else g * held[c])
            in_gap = True
    var b = band
    var steps = List[Int]()
    var mine = List[UInt8]()
    while True:
        var width = 2 * b + 1
        var cells = (columns + 1) * width
        cost.resize(3 * cells, BIG)
        moves.resize(3 * cells, MATCHED)
        # State t's cell (c, j) sits at t * cells + c * width + j - (mid[c] - b).
        for c in range(columns + 1):
            var lo = mid[c] - b
            var base = c * width
            var x = c - 1
            var lo_prev = mid[c - 1] - b if c > 0 else 0
            var prev = (c - 1) * width
            # Opening against the rows holding a letter in this column (before the first: every row).
            var held_here = held[x] if c > 0 else others
            for k in range(width):
                var j = lo + k
                var m_best = BIG
                var m_how = MATCHED
                var g_best = BIG
                var g_how = MATCHED
                var o_best = BIG
                var o_how = MATCHED
                if j >= 0 and j <= n:
                    if c == 0 and j == 0:
                        m_best = 0
                    if c > 0:
                        # M: letter j - 1 in column x.
                        if j > 0:
                            var km = j - 1 - lo_prev
                            if km >= 0 and km < width:
                                var unit = others - counts[x * SYMBOLS + _symbol(letters[j - 1])]
                                var from_m = cost[prev + km]
                                var from_g = cost[cells + prev + km]
                                var from_o = cost[2 * cells + prev + km]
                                if from_m < BIG and from_m + g * starts[x] < m_best:
                                    m_best = from_m + g * starts[x]
                                    m_how = MATCHED
                                if from_g < BIG and from_g + g * starts[x] < m_best:
                                    m_best = from_g + g * starts[x]
                                    m_how = GAPPED
                                if from_o < BIG and from_o < m_best:
                                    m_best = from_o
                                    m_how = OWN
                                if m_best < BIG:
                                    m_best += unit
                        # G: the sequence gapped in column x.
                        var kp = j - lo_prev
                        if kp >= 0 and kp < width:
                            var from_m = cost[prev + kp]
                            var from_g = cost[cells + prev + kp]
                            var from_o = cost[2 * cells + prev + kp]
                            if from_m < BIG and from_m + g * held[x] < g_best:
                                g_best = from_m + g * held[x]
                                g_how = MATCHED
                            if from_g < BIG and from_g + g * ends[x] < g_best:
                                g_best = from_g + g * ends[x]
                                g_how = GAPPED
                            if from_o < BIG and from_o + g * held[x] < g_best:
                                g_best = from_o + g * held[x]
                                g_how = OWN
                            if g_best < BIG:
                                g_best += held[x]
                    # O: letter j - 1 in a column of its own after column x.
                    if j > 0 and k > 0:
                        var from_m = cost[base + k - 1]
                        var from_g = cost[cells + base + k - 1]
                        var from_o = cost[2 * cells + base + k - 1]
                        if from_m < BIG and from_m + g * held_here < o_best:
                            o_best = from_m + g * held_here
                            o_how = MATCHED
                        if from_g < BIG and from_g + g * held_here < o_best:
                            o_best = from_g + g * held_here
                            o_how = GAPPED
                        if from_o < BIG and from_o < o_best:
                            o_best = from_o
                            o_how = OWN
                        if o_best < BIG:
                            o_best += others
                cost[base + k] = m_best
                moves[base + k] = m_how
                cost[cells + base + k] = g_best
                moves[cells + base + k] = g_how
                cost[2 * cells + base + k] = o_best
                moves[2 * cells + base + k] = o_how
        var end = n - (mid[columns] - b)
        var final_state = 0
        var final_cost = BIG
        if end >= 0 and end < width:
            for t in range(3):
                if cost[t * cells + columns * width + end] < final_cost:
                    final_cost = cost[t * cells + columns * width + end]
                    final_state = t
        if final_cost >= BIG:
            if b >= max_band or b >= max(n, columns):
                return
            b = min(2 * b, max_band)
            continue
        steps.clear()
        mine.clear()
        var tight = False
        var c = columns
        var j = n
        var state = final_state
        while c > 0 or j > 0:
            var k = j - (mid[c] - b)
            if (k < EDGE and mid[c] - b > 0) or (width - 1 - k < EDGE and mid[c] + b < n):
                tight = True
            var came = Int(moves[state * cells + c * width + k])
            if state == 0:
                steps.append(c - 1)
                mine.append(letters[j - 1])
                c -= 1
                j -= 1
            elif state == 1:
                steps.append(c - 1)
                mine.append(GAP)
                c -= 1
            else:
                steps.append(-1)
                mine.append(letters[j - 1])
                j -= 1
            state = came
        if tight and b < max_band and b < max(n, columns):
            b = min(2 * b, max_band)
            continue
        if final_cost > present:
            return
        break
    var fresh = List[List[UInt8]]()
    for r in range(len(rows)):
        var row = List[UInt8](capacity=len(steps))
        for t in reversed(range(len(steps))):
            var old = steps[t]
            if r == s:
                row.append(mine[t])
            elif old >= 0:
                row.append(rows[r][old])
            else:
                row.append(GAP)
        fresh.append(row^)
    var kept = List[Int]()
    for t in range(len(steps)):
        for r in range(len(fresh)):
            if fresh[r][t] != GAP:
                kept.append(t)
                break
    if len(kept) < len(steps):
        for r in range(len(fresh)):
            var row = List[UInt8](capacity=len(kept))
            for t in kept:
                row.append(fresh[r][t])
            fresh[r] = row^
    rows = fresh^


def refine(
    rows: List[List[UInt8]],
    passes: Int,
    band: Int = 64,
    max_band: Int = 4096,
    anchored: Bool = True,
    anchor_costs: Costs = Costs.edit(),
    gap_open: Int = 0,
    seeded: Bool = False,
) raises -> List[List[UInt8]]:
    """`rows` (one length, `-` gaps) refined on their sum of pairs for at most `passes` passes over the distinct
    sequences, each within `band` letters of its place, widened up to `max_band` while its best path touches the
    band's edge; with `gap_open`, every gap a pair opens costs that much more (see `_realign`). Identical
    sequences share their row. The rows keep their order and spell what they did."""
    var n = len(rows)
    if n < 2 or passes <= 0:
        return rows.copy()
    # One row per distinct sequence.
    var spelled = List[List[UInt8]]()
    var distinct = List[Int]()
    var row_of = List[Int](length=n, fill=-1)
    for r in range(n):
        var letters = List[UInt8]()
        for x in rows[r]:
            if x != GAP:
                letters.append(x)
        var found = -1
        for d in range(len(distinct)):
            if spelled[distinct[d]] == letters:
                found = d
                break
        spelled.append(letters^)
        if found < 0:
            row_of[r] = len(distinct)
            distinct.append(r)
        else:
            row_of[r] = found
    var work = List[List[UInt8]]()
    for d in distinct:
        work.append(rows[d].copy())
    if len(work) < 2:
        return rows.copy()
    # Columns gapped in every distinct row go first.
    var keep = List[Int]()
    for c in range(len(work[0])):
        for r in range(len(work)):
            if work[r][c] != GAP:
                keep.append(c)
                break
    for r in range(len(work)):
        var row = List[UInt8](capacity=len(keep))
        for c in keep:
            row.append(work[r][c])
        work[r] = row^
    var cost = List[Int32]()
    var moves = List[UInt8]()
    var centre = List[Int]()
    var spelled_distinct = List[List[UInt8]]()
    for d in distinct:
        spelled_distinct.append(spelled[d].copy())
    if not anchored:
        spelled_distinct.clear()
    var partners = Partners(spelled_distinct, anchor_costs, seeded)
    var last = sum_of_pairs(work) + gap_open * gap_opens(work) if gap_open > 0 else sum_of_pairs(work)
    for _ in range(passes):
        for s in range(len(work)):
            centre.clear()
            if anchored:
                _anchored_centre(work, s, partners, centre)
            _realign(work, s, band, max_band, centre, Int32(gap_open), cost, moves)
        var now = sum_of_pairs(work) + gap_open * gap_opens(work) if gap_open > 0 else sum_of_pairs(work)
        # A pass gaining under 0.1% ends the refinement: the later passes of a long run gain little.
        if now >= last or (last - now) * 1000 < last:
            break
        last = now
    var out = List[List[UInt8]]()
    for r in range(n):
        out.append(work[row_of[r]].copy())
    return out^
