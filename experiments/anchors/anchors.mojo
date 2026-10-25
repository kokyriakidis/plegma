"""Seed anchors against the truth: which anchors pair letters no optimal alignment pairs, and what marks them.

usage: anchors <sequences.fa> <period> <max pairs> <out.tsv> <plan> [<plan> ...]
       anchors --compose <sequences.fa> <a> <b>
       anchors --compose2 <sequences.fa> <a> <b> <period>
       anchors --brute <sequences.fa> <a> <b> <period> [bin] [--prices x,o,e]
       anchors --seed-audit <dump> [--prices x,o,e]
       anchors --coopt <pairs.tsv> [--two-piece x,o1,e1,o2,e2]
       anchors --coopt-set <pairs.tsv> [--two-piece x,o1,e1,o2,e2]

For an even sample of at most `max pairs` pairs of distinct sequences (the longer as the partner, as placement
seeds a new sequence against a longer placed one), each plan seeds the pair as plegma's backbone does, traced: every
anchor with the depth of the recursion and the level that found it, its length, its string's copies in each piece,
the window it was found in, whether the partner's repeat mask covers it, how many anchors its window chained, and
its diagonal against the window's two ends. Each anchored letter pair (past the trim) is checked against every
optimal unit-cost alignment of the pair, vntr-eval's metric: on one exactly when its forward cost, its own and its
backward cost sum to the optimum. A letter pair on none is off; its offset is the distance to the nearest partner
an optimal alignment gives the letter, in whole periods (`copies`), within two letters (`near`) or else (`other`).

A plan is `levels[;option ...]`: levels `mum51,mum31,mum21` (MUMs), `smem31` (SMEMs), `ref31` (unique in the
partner), `all31` (every MEM); a level may end in `:m` (its anchors kept out of the partner's tandem repeats at every depth: the flanks) and `:e`
(inside the repeats, the one-indel rule below); with either, a plan's mask and rule apply per level, not per
depth; options `stay` (a level again inside its own gaps until one yields nothing), `M=k:d` (no anchor in the
partner's tandem repeats: k-mers recurring within d letters), `D=n` (that mask at the first n depths only), `T=n`
(letters an anchor gives back at each end, 10), `B=n` (no seeding inside a window of n letters or fewer: left to the
exact aligner), `ends=e` (past the mask's depths, an anchor in the repeat is kept only if its diagonal lies within e
letters of one of its window's ends: the window's length difference all after it or all before it, one indel).
Or `optimal` (the optimal two-piece alignment, 4,4,2,24,1, what placement's last source uses) or `edit` (the
optimal unit-cost one): every matched letter pair, as runs of matches. Or `copies[;min=L]`: the flanks' MUMs (51,
31, 21, masked in the repeats at every depth), then inside the arrays (`units.repeat_at` at the period, cut into
copies against the region's unit, the longest sequence's consensus) the copies paired by `units.pair_copies`, and
inside each pair of copies the runs of at least L (1) matching letters of their optimal unit-cost alignment.

Being off an optimal alignment is a strict test in a repeat: an alignment a fraction of a percent worse may pair
most letters a copy apart. So each plan is also charged what its anchors force: the pair's best alignment through
every anchored letter pair (each piece between consecutive anchors aligned optimally, as vntr-eval scores an MSA's
shared letters) minus its optimum, summed over the pairs, against the summed optima.

Writes one TSV line per anchor (plan, pair, its features, letters, off letters by offset class) and prints per plan
the anchors, the letters and how many are off, the excess the anchors force, and per depth the anchors and letters
off.
"""

from std.os import getenv
from std.sys import argv

from karan import Aligner, Costs, Mode, Ties
from sporos import read_fasta

from plegma.anchors import Anchor, find_anchors, mums, mums_or_smems, smems
from plegma.backbone import ALL, Level, MUMS, REFERENCE_UNIQUE, SMEMS, RepeatMask, split_anchors
from plegma.chain import ChainCosts, chain, without_overlaps
from plegma.graph import DELETION, INSERTION, cigar_runs
from plegma.units import Copies, Repeat, consensus_unit, copies, pair_copies, repeat_at


struct Plan(Copyable, Movable):
    var text: String
    var levels: List[Level]
    var level_mask: List[Bool]
    """Per level (`:m`): no anchor in the partner's tandem repeats, at any depth."""
    var level_ends: List[Bool]
    """Per level (`:e`): inside the repeats, the one-indel rule (`ends`, 2 when not given)."""
    var per_level: Bool
    """Some level carries `:m` or `:e`: the mask and the rule apply per level."""
    var len_d: Float64
    """`split`'s `lenD=f`: a window's match must hold f times its length difference letters or more."""
    var trim_mask: Bool
    """`trimflanks`: a masked level's match reaching into a repeat cut back to its longest unmasked end (the flank
    part, up to the repeat's boundary) instead of dropped."""
    var max_depth: Int
    """`split`'s `maxdepth=k`: no splitting past depth k."""
    var stay: Bool
    var trim: Int
    var mask_k: Int
    var mask_distance: Int
    var mask_depth: Int
    var budget: Int
    var ends: Int
    var min_run: Int
    var kind: String
    """"seeds", "optimal", "edit", "copies" or "msa"."""
    var msa_path: String
    var rows: List[List[UInt8]]
    """With `msa`, each distinct sequence's row of that MSA."""

    def __init__(out self, text: String) raises:
        self.text = text
        self.ends = -1
        self.min_run = 1
        self.kind = "seeds"
        self.msa_path = String()
        self.rows = List[List[UInt8]]()
        self.levels = List[Level]()
        self.level_mask = List[Bool]()
        self.level_ends = List[Bool]()
        self.per_level = False
        self.len_d = 0.0
        self.trim_mask = False
        self.max_depth = Int.MAX
        self.stay = False
        self.trim = 10
        self.mask_k = 0
        self.mask_distance = 0
        self.mask_depth = Int.MAX
        self.budget = 0
        var parts = text.split(";")
        for k in range(1, len(parts)):
            var part = String(parts[k])
            if part == "stay":
                self.stay = True
            elif part.startswith("T="):
                self.trim = Int(String(part[byte=2:]))
            elif part.startswith("D="):
                self.mask_depth = Int(String(part[byte=2:]))
            elif part.startswith("B="):
                self.budget = Int(String(part[byte=2:]))
            elif part.startswith("ends="):
                self.ends = Int(String(part[byte=5:]))
            elif part.startswith("min="):
                self.min_run = Int(String(part[byte=4:]))
            elif (
                part.startswith("stop=")
                or part == "smems"
                or part == "flanks"
                or part.startswith("rule=")
                or part.startswith("lenD=")
                or part == "trimflanks"
                or part == "arrays"
                or part.startswith("maxdepth=")
            ):
                pass  # the greedy plan's, read with it below
            elif part.startswith("M="):
                var halves = String(part[byte=2:]).split(":")
                self.mask_k = Int(String(halves[0]))
                self.mask_distance = Int(String(halves[1]))
            else:
                raise Error(t"unknown option {part} in {text}")
        if String(parts[0]) == "optimal" or String(parts[0]) == "edit":
            self.kind = String(parts[0])
            self.levels.append(Level(MUMS, 1))
            return
        if String(parts[0]).startswith("msa="):
            # An MSA's stacked letters: a pair's anchors are the runs of columns holding one letter in both rows.
            self.kind = "msa"
            self.msa_path = String(String(parts[0])[byte=4:])
            self.trim = 0
            self.levels.append(Level(MUMS, 1))
            return
        if String(parts[0]) == "libsplit":
            # plegma's own `backbone.split_anchors`, its defaults.
            self.kind = "libsplit"
            self.trim = 0
            self.levels.append(Level(MUMS, 21))
            return
        if String(parts[0]) == "split":
            # One anchor at a time: in each window, the longest match unique in both of its pieces (frequencies
            # counted in the window alone), then the same inside the two windows it leaves; `flanks`, the flanks'
            # MUMs (outside the repeats) first; `min=L`, matches of L letters or more (21).
            self.kind = "split"
            self.trim = 0
            self.min_run = 21
            for k in range(1, len(parts)):
                var part = String(parts[k])
                if part.startswith("min="):
                    self.min_run = Int(String(part[byte=4:]))
                elif part == "flanks":
                    self.mask_k = 12
                    self.mask_distance = 500
                    for f in [51, 31, 21]:
                        self.levels.append(Level(MUMS, f, 1, True))
                elif part == "rule=inside":
                    self.ends = 1_000_000
                elif part.startswith("lenD="):
                    self.len_d = Float64(String(part[byte=5:]))
                elif part == "trimflanks":
                    self.trim_mask = True
                elif part == "arrays":
                    self.kind = "split_arrays"
                elif part.startswith("maxdepth="):
                    self.max_depth = Int(String(part[byte=9:]))
                elif part.startswith("rule=ends"):
                    self.ends = 2
            if len(self.levels) == 0:
                self.levels.append(Level(MUMS, 21))
            return
        if String(parts[0]) == "greedy":
            # Every MUM of 21 letters or more (unique in both sequences), the longest first, each kept if it lies
            # before or after every one kept already on both sequences; past stop=k, only those of k periods or
            # more (stop=0: all).
            self.kind = "greedy"
            self.trim = 0
            self.levels.append(Level(MUMS, 21))
            for k in range(1, len(parts)):
                var part = String(parts[k])
                if part.startswith("stop="):
                    self.min_run = Int(String(part[byte=5:]))
                elif part == "smems":
                    # The candidates: SMEMs as well as MUMs, ranked together by length.
                    self.stay = True
            return
        if String(parts[0]) == "copymums":
            # The flanks' MUMs, then whole copies whose string is unique in both arrays (copy-level MUMs), the most
            # that rise in both copy orders.
            self.kind = "copymums"
            self.trim = 0
            self.mask_k = 12
            self.mask_distance = 500
            for k in [51, 31, 21]:
                self.levels.append(Level(MUMS, k))
            self.stay = True
            return
        if String(parts[0]) == "copies":
            # The flanks' MUMs, masked at every depth; the copies' anchors are added to them.
            self.kind = "copies"
            self.trim = 0
            self.mask_k = 12
            self.mask_distance = 500
            for k in [51, 31, 21]:
                self.levels.append(Level(MUMS, k))
            self.stay = True
            return
        for item in String(parts[0]).split(","):
            var token = String(item)
            var mask = False
            var ends = False
            while token.endswith(":m") or token.endswith(":e"):
                if token.endswith(":m"):
                    mask = True
                else:
                    ends = True
                var shorter = String(token[byte = 0 : token.byte_length() - 2])
                token = shorter^
            self.level_mask.append(mask)
            self.level_ends.append(ends)
            if mask or ends:
                self.per_level = True
                if self.mask_k == 0:
                    self.mask_k = 12
                    self.mask_distance = 500
                if self.ends < 0 and ends:
                    self.ends = 2
            if token.startswith("mum"):
                self.levels.append(Level(MUMS, Int(String(token[byte=3:]))))
            elif token.startswith("smem"):
                self.levels.append(Level(SMEMS, Int(String(token[byte=4:]))))
            elif token.startswith("ref"):
                self.levels.append(Level(REFERENCE_UNIQUE, Int(String(token[byte=3:]))))
            elif token.startswith("all"):
                self.levels.append(Level(ALL, Int(String(token[byte=3:]))))
            else:
                raise Error(t"unknown level {token} in {text}")


@fieldwise_init
struct Traced(Copyable, Movable):
    var ref_pos: Int
    var query_pos: Int
    var length: Int
    var depth: Int
    var level: Int
    var ref_copies: Int
    var query_copies: Int
    var win_ref: Int
    var win_query: Int
    var in_repeat: Bool
    var siblings: Int
    var diag_start: Int
    """The anchor's diagonal against its window's start: query letters before it minus partner letters."""
    var diag_end: Int
    """Likewise against the window's end: query letters after it minus partner letters."""


def _trace(
    a: List[UInt8],
    b: List[UInt8],
    r0: Int,
    r1: Int,
    q0: Int,
    q1: Int,
    plan: Plan,
    level: Int,
    depth: Int,
    repeats: RepeatMask,
    mut out: List[Traced],
) raises:
    """`backbone`'s recursion over window [r0, r1) x [q0, q1), every anchor kept traced."""
    if level >= len(plan.levels) or max(r1 - r0, q1 - q0) <= plan.budget:
        return
    var current = plan.levels[level]
    var found = current.select(find_anchors(Span(a)[r0:r1], Span(b)[q0:q1], current.min_len, current.counted_to()))
    var per_level = plan.per_level
    var masked = per_level and plan.level_mask[level]
    var one_indel = per_level and plan.level_ends[level]
    if repeats.k > 0 and (masked or (depth < plan.mask_depth and not per_level)):
        var outside = List[Anchor]()
        for x in found:
            if not repeats.covers(r0 + x.ref_pos, r0 + x.ref_end()):
                outside.append(x)
            elif plan.trim_mask:
                # Cut back to the longest end clear of the repeat: up to its boundary.
                var lo = 0
                var hi = x.length
                while lo < hi:
                    var mid = (lo + hi + 1) // 2
                    if repeats.covers(r0 + x.ref_pos, r0 + x.ref_pos + mid):
                        hi = mid - 1
                    else:
                        lo = mid
                var prefix = lo
                lo = 0
                hi = x.length
                while lo < hi:
                    var mid = (lo + hi + 1) // 2
                    if repeats.covers(r0 + x.ref_end() - mid, r0 + x.ref_end()):
                        hi = mid - 1
                    else:
                        lo = mid
                var suffix = lo
                if max(prefix, suffix) >= current.min_len:
                    if prefix >= suffix:
                        outside.append(Anchor(x.ref_pos, x.query_pos, prefix, x.ref_copies, x.query_copies))
                    else:
                        var cut = x.length - suffix
                        outside.append(Anchor(x.ref_pos + cut, x.query_pos + cut, suffix, x.ref_copies, x.query_copies))
        found = outside^
    elif repeats.k > 0 and plan.ends >= 0 and (one_indel or not per_level):
        # Inside the repeat, only anchors leaving the window's length difference to one side of them.
        var kept = List[Anchor]()
        for x in found:
            var start = x.query_pos - x.ref_pos
            var end = ((q1 - q0) - x.query_end()) - ((r1 - r0) - x.ref_end())
            if not repeats.covers(r0 + x.ref_pos, r0 + x.ref_end()) or abs(start) <= plan.ends or abs(end) <= plan.ends:
                kept.append(x)
        found = kept^
    var picked = List[Anchor]()
    if len(found) <= 50_000:
        picked = without_overlaps(chain(found, r1 - r0, q1 - q0, ChainCosts(substitution=0.1), 0).anchors)
    var next = level if plan.stay else level + 1
    if len(picked) == 0:
        _trace(a, b, r0, r1, q0, q1, plan, level + 1, depth, repeats, out)
        return
    var ref_at = 0
    var query_at = 0
    for x in picked:
        _trace(a, b, r0 + ref_at, r0 + x.ref_pos, q0 + query_at, q0 + x.query_pos, plan, next, depth + 1, repeats, out)
        out.append(
            Traced(
                r0 + x.ref_pos,
                q0 + x.query_pos,
                x.length,
                depth,
                level,
                x.ref_copies,
                x.query_copies,
                r1 - r0,
                q1 - q0,
                repeats.covers(r0 + x.ref_pos, r0 + x.ref_end()),
                len(picked),
                x.query_pos - x.ref_pos,
                ((q1 - q0) - x.query_end()) - ((r1 - r0) - x.ref_end()),
            )
        )
        ref_at = x.ref_end()
        query_at = x.query_end()
    _trace(a, b, r0 + ref_at, r1, q0 + query_at, q1, plan, next, depth + 1, repeats, out)


def _aligned(a: List[UInt8], b: List[UInt8], unit: Bool, mut out: List[Traced]) raises:
    """The matched letter pairs of the pair's optimal alignment (two-piece, or unit costs), as runs of matches."""
    var aligner = Aligner()
    var costs = Costs.edit() if unit else Costs.two_piece(4, 4, 2, 24, 1)
    var found = aligner.align(Span(a), Span(b), costs, Mode.GLOBAL, ties=Ties.LEFT)
    var runs = cigar_runs(found.cigar)
    var i = 0
    var j = 0
    var run_i = -1
    var run_j = -1
    var run = 0

    @parameter
    def close():
        if run > 0:
            out.append(Traced(run_i, run_j, run, 0, 0, 1, 1, len(a), len(b), False, 0, run_j - run_i, 0))

    for x in range(len(runs.ops)):
        for _ in range(runs.counts[x]):
            var op = runs.ops[x]
            if op == DELETION:
                close()
                run = 0
                i += 1
            elif op == INSERTION:
                close()
                run = 0
                j += 1
            else:
                if a[i] == b[j]:
                    if run == 0:
                        run_i = i
                        run_j = j
                    run += 1
                else:
                    close()
                    run = 0
                i += 1
                j += 1
    close()


def _split(
    a: List[UInt8],
    b: List[UInt8],
    r0: Int,
    r1: Int,
    q0: Int,
    q1: Int,
    min_len: Int,
    rule: Int,
    len_d: Float64,
    max_depth: Int,
    depth: Int,
    mut out: List[Traced],
) raises:
    """In window [r0, r1) x [q0, q1): the longest match of `min_len` letters or more unique in both pieces (its
    frequencies counted in the window alone) that `rule` accepts, anchored; then the same in the windows left of it
    and right of it. `rule` >= 0: the match's diagonal must lie within the window's (between its two ends, so the
    window's length difference D is taken by gaps on one side or both, never by opposing ones); below 1,000,000,
    also within `rule` letters of one end (D all to one side: one indel). -1: any."""
    if r1 - r0 < min_len or q1 - q0 < min_len or depth > max_depth:
        return
    var d = (q1 - q0) - (r1 - r0)
    # The larger the window's length difference, the more letters a match must hold (`len_d`).
    var need = max(min_len, Int(len_d * Float64(abs(d))))
    var found = mums(find_anchors(Span(a)[r0:r1], Span(b)[q0:q1], need, 2))
    var best = -1
    for k in range(len(found)):
        var start = found[k].query_pos - found[k].ref_pos
        if rule >= 0:
            if start < min(0, d) - 2 or start > max(0, d) + 2:
                continue
            if rule < 1_000_000 and abs(start) > rule and abs(d - start) > rule:
                continue
        if best < 0 or found[k].length > found[best].length:
            best = k
    if best < 0:
        return
    var x = found[best]
    _split(a, b, r0, r0 + x.ref_pos, q0, q0 + x.query_pos, min_len, rule, len_d, max_depth, depth + 1, out)
    out.append(
        Traced(
            r0 + x.ref_pos,
            q0 + x.query_pos,
            x.length,
            depth,
            0,
            1,
            1,
            r1 - r0,
            q1 - q0,
            False,
            len(found),
            x.query_pos - x.ref_pos,
            ((q1 - q0) - x.query_end()) - ((r1 - r0) - x.ref_end()),
        )
    )
    _split(a, b, r0 + x.ref_end(), r1, q0 + x.query_end(), q1, min_len, rule, len_d, max_depth, depth + 1, out)


def _greedy(a: List[UInt8], b: List[UInt8], min_len: Int, with_smems: Bool, mut out: List[Traced]) raises:
    """The pair's MUMs (21 letters or more, `min_len` or more; with `with_smems`, its SMEMs too, every occurrence),
    the longest first, each kept if it lies before or after every one kept already on both sequences: the longest
    pin the backbone, the shorter fill what is left."""
    var found = mums_or_smems(find_anchors(Span(a), Span(b), max(21, min_len), 2)) if with_smems else mums(
        find_anchors(Span(a), Span(b), max(21, min_len), 2)
    )
    var keys = List[UInt64]()
    for k in range(len(found)):
        keys.append((UInt64((1 << 32) - found[k].length) << 32) | UInt64(k))
    sort(keys)
    var kept = List[Anchor]()
    for key in keys:
        var x = found[Int(key & 0xFFFFFFFF)]
        var fits = True
        for y in kept:
            var before = x.ref_end() <= y.ref_pos and x.query_end() <= y.query_pos
            var after = x.ref_pos >= y.ref_end() and x.query_pos >= y.query_end()
            if not before and not after:
                fits = False
                break
        if fits:
            kept.append(x)
    for x in kept:
        out.append(Traced(x.ref_pos, x.query_pos, x.length, 0, 0, 1, 1, len(a), len(b), False, len(kept), 0, 0))


def _copy_mums(a: List[UInt8], b: List[UInt8], period: Int, unit: List[UInt8], mut out: List[Traced]) raises:
    """Whole copies whose string is once in each array, the most that rise in both copy orders (a longest
    increasing chain), each an anchor of its letters; those crossing an anchor already in `out` dropped."""
    var ra = repeat_at(a, period)
    var rb = repeat_at(b, period)
    if ra.period == 0 or rb.period == 0:
        return
    var ca = copies(a, ra, unit)
    var cb = copies(b, rb, unit)
    var count_a = Dict[String, Int]()
    var count_b = Dict[String, Int]()
    var where_b = Dict[String, Int]()
    for c in range(len(ca.starts)):
        var key = String(unsafe_from_utf8=List[UInt8](Span(a)[ca.starts[c] : ca.ends[c]]))
        count_a[key] = count_a.get(key, 0) + 1
    for d in range(len(cb.starts)):
        var key = String(unsafe_from_utf8=List[UInt8](Span(b)[cb.starts[d] : cb.ends[d]]))
        count_b[key] = count_b.get(key, 0) + 1
        where_b[key] = d
    var pa = List[Int]()
    var pb = List[Int]()
    for c in range(len(ca.starts)):
        var key = String(unsafe_from_utf8=List[UInt8](Span(a)[ca.starts[c] : ca.ends[c]]))
        if count_a.get(key, 0) == 1 and count_b.get(key, 0) == 1:
            pa.append(c)
            pb.append(where_b[key])
    # The longest chain rising in b's copy order (a's already rises): patience sorting.
    var tails = List[Int]()
    var back = List[Int](length=len(pb), fill=-1)
    for t in range(len(pb)):
        var lo = 0
        var hi = len(tails)
        while lo < hi:
            var mid = (lo + hi) // 2
            if pb[tails[mid]] < pb[t]:
                lo = mid + 1
            else:
                hi = mid
        if lo > 0:
            back[t] = tails[lo - 1]
        if lo == len(tails):
            tails.append(t)
        else:
            tails[lo] = t
    var keep = List[Int]()
    var t = tails[len(tails) - 1] if len(tails) > 0 else -1
    while t >= 0:
        keep.append(t)
        t = back[t]
    var flanks = len(out)
    for q in reversed(range(len(keep))):
        var u = keep[q]
        var i0 = ca.starts[pa[u]]
        var j0 = cb.starts[pb[u]]
        var length = ca.ends[pa[u]] - i0
        var crosses = False
        for f in range(flanks):
            ref fa = out[f]
            var before = i0 + length <= fa.ref_pos and j0 + length <= fa.query_pos
            var after = i0 >= fa.ref_pos + fa.length and j0 >= fa.query_pos + fa.length
            if not before and not after:
                crosses = True
                break
        if not crosses:
            out.append(Traced(i0, j0, length, 1, 0, 1, 1, len(a), len(b), True, len(keep), j0 - i0, 0))


def _copy_anchors(
    a: List[UInt8], b: List[UInt8], period: Int, unit: List[UInt8], min_run: Int, mut out: List[Traced]
) raises:
    """Inside the two arrays, the copies paired (`pair_copies`) and the matching runs of each pair's optimal
    unit-cost alignment, at least `min_run` long; those crossing an anchor already in `out` (the flanks') dropped."""
    var ra = repeat_at(a, period)
    var rb = repeat_at(b, period)
    if ra.period == 0 or rb.period == 0:
        return
    var ca = copies(a, ra, unit)
    var cb = copies(b, rb, unit)
    var aligner = Aligner()
    var pairs = pair_copies(a, ca, b, cb, aligner)
    var flanks = len(out)
    for pr in pairs:
        var i0 = ca.starts[pr[0]]
        var j0 = cb.starts[pr[1]]
        var found = aligner.align(
            Span(a)[i0 : ca.ends[pr[0]]], Span(b)[j0 : cb.ends[pr[1]]], Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT
        )
        var runs = cigar_runs(found.cigar)
        var i = i0
        var j = j0
        var run_i = -1
        var run_j = -1
        var run = 0
        for x in range(len(runs.ops) + 1):
            var count = runs.counts[x] if x < len(runs.ops) else 1
            for _ in range(count):
                var matched = (
                    x < len(runs.ops) and runs.ops[x] != DELETION and runs.ops[x] != INSERTION and a[i] == b[j]
                )
                if matched:
                    if run == 0:
                        run_i = i
                        run_j = j
                    run += 1
                elif run >= min_run:
                    # Kept unless it crosses a flank anchor (not before or after it on both sequences).
                    var crosses = False
                    for f in range(flanks):
                        ref fa = out[f]
                        var before = run_i + run <= fa.ref_pos and run_j + run <= fa.query_pos
                        var after = run_i >= fa.ref_pos + fa.length and run_j >= fa.query_pos + fa.length
                        if not before and not after:
                            crosses = True
                            break
                    if not crosses:
                        out.append(
                            Traced(run_i, run_j, run, 1, 0, 1, 1, len(a), len(b), True, len(pairs), run_j - run_i, 0)
                        )
                    run = 0
                else:
                    run = 0
                if x >= len(runs.ops):
                    break
                var op = runs.ops[x]
                if op == DELETION:
                    i += 1
                elif op == INSERTION:
                    j += 1
                else:
                    i += 1
                    j += 1


def _stacked(ra: List[UInt8], rb: List[UInt8], mut out: List[Traced]):
    """The runs of columns where rows `ra` and `rb` hold one letter: what vntr-eval's graph shares between the two."""
    var gap = UInt8(ord("-"))
    var i = 0
    var j = 0
    var run_i = -1
    var run_j = -1
    var run = 0
    for c in range(len(ra) + 1):
        var same = c < len(ra) and ra[c] != gap and ra[c] == rb[c]
        if same:
            if run == 0:
                run_i = i
                run_j = j
            run += 1
        elif run > 0:
            out.append(Traced(run_i, run_j, run, 0, 0, 1, 1, 0, 0, False, 0, run_j - run_i, 0))
            run = 0
        if c < len(ra):
            if ra[c] != gap:
                i += 1
            if rb[c] != gap:
                j += 1


def _compose(args: List[String]) raises:
    """`--compose <sequences.fa> <a> <b>`: for every third distinct sequence x, the optimal unit-cost alignments of a
    and of b to x composed through x's letters (the three-row MSA an MSA aligning both optimally to x implies), and
    the cost vntr-eval would charge the pair (a, b) there (each piece between their stacked letters aligned
    optimally) against their optimum: whether some x lets a and b keep their optimal alignment."""
    var records = read_fasta(args[2])
    var seqs = List[List[UInt8]]()
    for r in range(len(records)):
        var letters = List[UInt8](records.seq(r))
        for k in range(len(letters)):
            if letters[k] >= UInt8(ord("a")) and letters[k] <= UInt8(ord("z")):
                letters[k] -= 32
        var dup = False
        for x in seqs:
            if x == letters:
                dup = True
                break
        if not dup:
            seqs.append(letters^)
    var ia = Int(args[3])
    var ib = Int(args[4])
    ref a = seqs[ia]
    ref b = seqs[ib]
    var aligner = Aligner()
    var best = aligner.distance(Span(a), Span(b))
    print("pair", ia, ib, "optimum", best)
    for x in range(len(seqs)):
        if x == ia or x == ib:
            continue
        ref c = seqs[x]
        # Each letter of x: the letter of a, and of b, its optimal alignment pairs with it (-1: none).
        var to_a = List[Int](length=len(c), fill=-1)
        var to_b = List[Int](length=len(c), fill=-1)
        for side in range(2):
            var found = aligner.align(
                Span(c), Span(a) if side == 0 else Span(b), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT
            )
            var runs = cigar_runs(found.cigar)
            var i = 0
            var j = 0
            for t in range(len(runs.ops)):
                for _ in range(runs.counts[t]):
                    var op = runs.ops[t]
                    if op == DELETION:
                        i += 1
                    elif op == INSERTION:
                        j += 1
                    else:
                        if side == 0:
                            to_a[i] = j
                        else:
                            to_b[i] = j
                        i += 1
                        j += 1
        # a's and b's letters stacked through x: identical letters both paired with one letter of x.
        var cost = 0
        var pa = 0
        var pb = 0
        for i in range(len(c)):
            var ja = to_a[i]
            var jb = to_b[i]
            if ja >= pa and jb >= pb and a[ja] == b[jb]:
                cost += aligner.distance(Span(a)[pa:ja], Span(b)[pb:jb])
                pa = ja + 1
                pb = jb + 1
        cost += aligner.distance(Span(a)[pa : len(a)], Span(b)[pb : len(b)])
        print("  through", x, "cost", cost, "excess", cost - best)


def _homology_marks(
    a: List[UInt8], b: List[UInt8], period: Int, unit: List[UInt8]
) raises -> Tuple[List[Int], List[Int]]:
    """For each letter of a and of b, the index of the shared unique copy it lies in (the same index on both sides
    for a copy of a and the copy of b with its string, once in each array), or -1."""
    var mark_a = List[Int](length=len(a), fill=-1)
    var mark_b = List[Int](length=len(b), fill=-1)
    var ra = repeat_at(a, period)
    var rb = repeat_at(b, period)
    if ra.period == 0 or rb.period == 0:
        return (mark_a^, mark_b^)
    var ca = copies(a, ra, unit)
    var cb = copies(b, rb, unit)
    var count_a = Dict[String, Int]()
    var count_b = Dict[String, Int]()
    var where_b = Dict[String, Int]()
    for c in range(len(ca.starts)):
        var key = String(unsafe_from_utf8=List[UInt8](Span(a)[ca.starts[c] : ca.ends[c]]))
        count_a[key] = count_a.get(key, 0) + 1
    for d in range(len(cb.starts)):
        var key = String(unsafe_from_utf8=List[UInt8](Span(b)[cb.starts[d] : cb.ends[d]]))
        count_b[key] = count_b.get(key, 0) + 1
        where_b[key] = d
    var next = 0
    for c in range(len(ca.starts)):
        var key = String(unsafe_from_utf8=List[UInt8](Span(a)[ca.starts[c] : ca.ends[c]]))
        if count_a.get(key, 0) == 1 and count_b.get(key, 0) == 1:
            var d = where_b[key]
            for k in range(ca.starts[c], ca.ends[c]):
                mark_a[k] = next
            for k in range(cb.starts[d], cb.ends[d]):
                mark_b[k] = next
            next += 1
    return (mark_a^, mark_b^)


def _tie_broken(a: List[UInt8], b: List[UInt8], mark_a: List[Int], mark_b: List[Int]) -> List[Int]:
    """An optimal unit-cost alignment of a and b, of all of them the one matching the most letters of a with the
    same shared unique copy's letters of b (`_homology_marks`): each letter of a's partner in b (-1: deleted)."""
    var n = len(a)
    var m = len(b)
    var w = m + 1
    # Lexicographic: cost first, then the most homologous matches; one key, cost * SCALE - matches.
    comptime SCALE = 1 << 20
    var key = List[Int64](length=(n + 1) * w, fill=0)
    var move = List[UInt8](length=(n + 1) * w, fill=0)
    for j in range(1, w):
        key[j] = Int64(j) * SCALE
        move[j] = 2
    for i in range(1, n + 1):
        key[i * w] = Int64(i) * SCALE
        move[i * w] = 1
        for j in range(1, w):
            var same = a[i - 1] == b[j - 1]
            var bonus = Int64(1) if same and mark_a[i - 1] >= 0 and mark_a[i - 1] == mark_b[j - 1] else Int64(0)
            var d = key[(i - 1) * w + j - 1] + (Int64(0) if same else Int64(SCALE)) - bonus
            var u = key[(i - 1) * w + j] + SCALE
            var l = key[i * w + j - 1] + SCALE
            var best = d
            var how = UInt8(0)
            if u < best:
                best = u
                how = 1
            if l < best:
                best = l
                how = 2
            key[i * w + j] = best
            move[i * w + j] = how
    var partner = List[Int](length=n, fill=-1)
    var i = n
    var j = m
    while i > 0 or j > 0:
        var how = move[i * w + j]
        if how == 0:
            partner[i - 1] = j - 1
            i -= 1
            j -= 1
        elif how == 1:
            i -= 1
        else:
            j -= 1
    return partner^


def _compose2(args: List[String]) raises:
    """`--compose2 <sequences.fa> <a> <b> <period>`: `--compose` with every alignment tie-broken toward pairing
    shared unique copies (`_tie_broken`), and the pair's own optimum beside it."""
    var records = read_fasta(args[2])
    var seqs = List[List[UInt8]]()
    for r in range(len(records)):
        var letters = List[UInt8](records.seq(r))
        for k in range(len(letters)):
            if letters[k] >= UInt8(ord("a")) and letters[k] <= UInt8(ord("z")):
                letters[k] -= 32
        var dup = False
        for x in seqs:
            if x == letters:
                dup = True
                break
        if not dup:
            seqs.append(letters^)
    var ia = Int(args[3])
    var ib = Int(args[4])
    var period = Int(args[5])
    var longest = 0
    for k in range(1, len(seqs)):
        if len(seqs[k]) > len(seqs[longest]):
            longest = k
    var top = repeat_at(seqs[longest], period)
    var unit = consensus_unit(seqs[longest], top.start, top.end, period)
    ref a = seqs[ia]
    ref b = seqs[ib]
    var aligner = Aligner()
    var best = aligner.distance(Span(a), Span(b))
    print("pair", ia, ib, "optimum", best)
    var total_left = 0
    var total_tied = 0
    for x in range(len(seqs)):
        if x == ia or x == ib:
            continue
        ref c = seqs[x]
        var ma = _homology_marks(c, a, period, unit)
        var mb = _homology_marks(c, b, period, unit)
        var to_a = _tie_broken(c, a, ma[0], ma[1])
        var to_b = _tie_broken(c, b, mb[0], mb[1])
        var cost = 0
        var pa = 0
        var pb = 0
        for i in range(len(c)):
            var ja = to_a[i]
            var jb = to_b[i]
            if ja >= pa and jb >= pb and a[ja] == b[jb]:
                cost += aligner.distance(Span(a)[pa:ja], Span(b)[pb:jb])
                pa = ja + 1
                pb = jb + 1
        cost += aligner.distance(Span(a)[pa : len(a)], Span(b)[pb : len(b)])
        print("  through", x, "cost", cost, "excess", cost - best)
        total_tied += cost - best


comptime FAR = Int32(1 << 29)


struct Prices(ImplicitlyCopyable):
    """Gap-affine prices: a mismatch costs `mismatch`, a gap of k letters `open + extend * k` (unit cost: 1, 0, 1)."""

    var mismatch: Int32
    var open: Int32
    var extend: Int32

    def __init__(out self, mismatch: Int32, open: Int32, extend: Int32):
        self.mismatch = mismatch
        self.open = open
        self.extend = extend


struct Truth(Movable):
    """The full forward (D, E, F: ending in a diagonal step, a gap in a, a gap in b) and backward (H') DP of a pair
    under `Prices`, and which cells an optimal alignment crosses."""

    var n: Int
    var m: Int
    var w: Int
    var d: List[Int32]
    var e: List[Int32]
    var f: List[Int32]
    var hb: List[Int32]
    var best: Int32
    var lo_j: List[Int]
    """Per boundary i of a: the least j of a cell (i, j) some optimal alignment crosses."""
    var hi_j: List[Int]
    var count: List[Int]
    """Per boundary i of a: how many cells optimal alignments cross; a pinch has one."""
    var forced: List[Bool]
    """Letters of a every optimal alignment matches to the same letter of b."""
    var forced_letters: Int

    def __init__(out self, a: List[UInt8], b: List[UInt8], p: Prices):
        self.n = len(a)
        self.m = len(b)
        var n = self.n
        var m = self.m
        var w = m + 1
        self.w = w
        var cells = (n + 1) * w
        self.d = List[Int32](length=cells, fill=FAR)
        self.e = List[Int32](length=cells, fill=FAR)
        self.f = List[Int32](length=cells, fill=FAR)
        self.hb = List[Int32](length=cells, fill=FAR)
        var go = p.open + p.extend
        self.d[0] = 0
        for j in range(1, w):
            self.e[j] = p.open + p.extend * Int32(j)
        for i in range(1, n + 1):
            var row = i * w
            var up = row - w
            self.f[row] = p.open + p.extend * Int32(i)
            var x = a[i - 1]
            for j in range(1, w):
                var c = row + j
                var hd = min(self.d[up + j - 1], min(self.e[up + j - 1], self.f[up + j - 1]))
                self.d[c] = hd + (Int32(0) if x == b[j - 1] else p.mismatch)
                self.e[c] = min(min(self.d[c - 1], self.f[c - 1]) + go, self.e[c - 1] + p.extend)
                self.f[c] = min(min(self.d[up + j], self.e[up + j]) + go, self.f[up + j] + p.extend)
        var last = n * w + m
        self.best = min(self.d[last], min(self.e[last], self.f[last]))
        self.lo_j = List[Int](length=n + 1, fill=Int.MAX)
        self.hi_j = List[Int](length=n + 1, fill=Int.MIN)
        self.count = List[Int](length=n + 1, fill=0)
        # Backward: H' in full, E' (in a gap in a, so extending it costs `extend`) along the row and F' (in a gap in
        # b) one row at a time; a cell is crossed by an optimal alignment when min(D + H', E + E', F + F') is the
        # optimum.
        var fb = List[Int32](length=w, fill=FAR)
        for i in reversed(range(n + 1)):
            var row = i * w
            var eb = FAR
            for j in reversed(range(w)):
                var c = row + j
                var h: Int32
                if i == n and j == m:
                    h = 0
                else:
                    h = FAR
                    if i < n and j < m:
                        h = self.hb[c + w + 1] + (Int32(0) if a[i] == b[j] else p.mismatch)
                    if j < m:
                        h = min(h, eb + go)
                    if i < n:
                        h = min(h, fb[j] + go)
                self.hb[c] = h
                eb = min(h, eb + p.extend) if j < m else h
                fb[j] = min(h, fb[j] + p.extend) if i < n else h
                var through = min(self.d[c] + h, min(self.e[c] + eb, self.f[c] + fb[j]))
                if through == self.best:
                    self.lo_j[i] = min(self.lo_j[i], j)
                    self.hi_j[i] = max(self.hi_j[i], j)
                    self.count[i] += 1
        self.forced = List[Bool](length=n, fill=False)
        self.forced_letters = 0
        for i in range(n):
            if (
                self.count[i] == 1
                and self.count[i + 1] == 1
                and self.hi_j[i + 1] == self.lo_j[i] + 1
                and self.lo_j[i] < m
                and a[i] == b[self.lo_j[i]]
            ):
                self.forced[i] = True
                self.forced_letters += 1

    def h(self, i: Int, j: Int) -> Int32:
        var c = i * self.w + j
        return min(self.d[c], min(self.e[c], self.f[c]))

    def excess(self, i: Int, j: Int, length: Int) -> Int:
        """The best alignment through an exact match of `length` letters at (i, j), over the optimum."""
        return Int(self.h(i, j) + self.hb[(i + length) * self.w + j + length] - self.best)

    def matched(self, a: List[UInt8], b: List[UInt8], i: Int) -> Bool:
        """Whether some optimal alignment matches letter i of a to an equal letter of b."""
        for j in range(max(0, self.lo_j[i]), min(self.m, self.hi_j[i] + 1)):
            if a[i] == b[j] and self.h(i, j) + self.hb[(i + 1) * self.w + j + 1] == self.best:
                return True
        return False


def _path_rows(n: Int, m: Int, anchors: List[Anchor], mut jlo: List[Int], mut jhi: List[Int]):
    """Per boundary i of a, the least and greatest letter of b the path from (0, 0) through the anchors to (n, m)
    (straight lines between them) crosses."""
    var pi = List[Int]()
    var pj = List[Int]()
    pi.append(0)
    pj.append(0)
    for x in anchors:
        pi.append(x.ref_pos)
        pj.append(x.query_pos)
        pi.append(x.ref_pos + x.length)
        pj.append(x.query_pos + x.length)
    pi.append(n)
    pj.append(m)
    jlo = List[Int](length=n + 1, fill=Int.MAX)
    jhi = List[Int](length=n + 1, fill=Int.MIN)
    for s in range(len(pi) - 1):
        var i1 = pi[s]
        var j1 = pj[s]
        var i2 = pi[s + 1]
        var j2 = pj[s + 1]
        if i1 == i2:
            jlo[i1] = min(jlo[i1], min(j1, j2))
            jhi[i1] = max(jhi[i1], max(j1, j2))
            continue
        for i in range(i1, i2 + 1):
            jlo[i] = min(jlo[i], j1 + ((j2 - j1) * (i - i1)) // (i2 - i1))
            jhi[i] = max(jhi[i], j1 + ((j2 - j1) * (i - i1) + (i2 - i1) - 1) // (i2 - i1))


def _band_forward(
    a: List[UInt8],
    b: List[UInt8],
    lo: List[Int],
    hi: List[Int],
    p: Prices,
    mut alpha: List[Int32],
    mut beta: List[Int32],
) -> Int32:
    """The best alignment whose cells (i, j) all keep lo[i] <= j <= hi[i] (lo and hi never decreasing, the band
    holding (0, 0) and (n, m)); and per row i, over the band's exits there (the cells a step can leave it from: (i,
    hi[i]) above, (i, lo[i] .. lo[i + 1] - 1) below), the least of H - extend * d and of H + extend * d (d = j - i)
    into alpha and beta."""
    var n = len(a)
    var m = len(b)
    var go = p.open + p.extend
    var pd = List[Int32](length=m + 1, fill=FAR)
    var pe = List[Int32](length=m + 1, fill=FAR)
    var pf = List[Int32](length=m + 1, fill=FAR)
    var cd = List[Int32](length=m + 1, fill=FAR)
    var ce = List[Int32](length=m + 1, fill=FAR)
    var cf = List[Int32](length=m + 1, fill=FAR)
    alpha = List[Int32](length=n + 1, fill=FAR)
    beta = List[Int32](length=n + 1, fill=FAR)
    for j in range(lo[0], hi[0] + 1):
        if j == 0:
            pd[0] = 0
        else:
            pe[j] = p.open + p.extend * Int32(j)
    for i in range(n + 1):
        if i > 0:
            # The current rows still hold row i - 2.
            if i >= 2:
                for j in range(lo[i - 2], hi[i - 2] + 1):
                    cd[j] = FAR
                    ce[j] = FAR
                    cf[j] = FAR
            var x = a[i - 1]
            for j in range(lo[i], hi[i] + 1):
                if j > 0:
                    var hd = min(pd[j - 1], min(pe[j - 1], pf[j - 1]))
                    cd[j] = hd + (Int32(0) if x == b[j - 1] else p.mismatch)
                    ce[j] = min(min(cd[j - 1], cf[j - 1]) + go, ce[j - 1] + p.extend)
                cf[j] = min(min(pd[j], pe[j]) + go, pf[j] + p.extend)
            swap(pd, cd)
            swap(pe, ce)
            swap(pf, cf)
        # pd, pe, pf now hold row i.
        if hi[i] < m:
            var j = hi[i]
            var h = min(pd[j], min(pe[j], pf[j]))
            if h < FAR // 2:
                var dd = Int32(j - i) * p.extend
                alpha[i] = min(alpha[i], h - dd)
                beta[i] = min(beta[i], h + dd)
        if i < n:
            for j in range(lo[i], min(hi[i], lo[i + 1] - 1) + 1):
                var h = min(pd[j], min(pe[j], pf[j]))
                if h < FAR // 2:
                    var dd = Int32(j - i) * p.extend
                    alpha[i] = min(alpha[i], h - dd)
                    beta[i] = min(beta[i], h + dd)
    return min(pd[m], min(pe[m], pf[m]))


def _band_bound(
    a: List[UInt8], b: List[UInt8], lo: List[Int], hi: List[Int], p: Prices, alpha: List[Int32], beta: List[Int32]
) -> Int32:
    """A lower bound on every alignment that leaves the band: one that does leaves it first after a cell p and
    returns for the last time at a cell q (a row below p's): its cost is at least the band's best to p, plus
    `extend` per diagonal between p and q (a gap letter each), plus the band's best from q (backward, in any state).
    The least over such p and q, from `_band_forward`'s alpha and beta and a backward pass over the band."""
    var n = len(a)
    var m = len(b)
    var go = p.open + p.extend
    # Prefix least of alpha and beta over rows before i.
    var pa = List[Int32](length=n + 1, fill=FAR)
    var pb = List[Int32](length=n + 1, fill=FAR)
    for i in range(1, n + 1):
        pa[i] = min(pa[i - 1], alpha[i - 1])
        pb[i] = min(pb[i - 1], beta[i - 1])
    var nh = List[Int32](length=m + 2, fill=FAR)
    var nf = List[Int32](length=m + 2, fill=FAR)
    var ch = List[Int32](length=m + 2, fill=FAR)
    var ce = List[Int32](length=m + 2, fill=FAR)
    var cf = List[Int32](length=m + 2, fill=FAR)
    var bound = FAR
    for i in reversed(range(n + 1)):
        # ce still holds row i + 1 (it runs along a row), ch and cf row i + 2.
        if i + 1 <= n:
            for j in range(lo[i + 1], hi[i + 1] + 1):
                ce[j] = FAR
        if i + 2 <= n:
            for j in range(lo[i + 2], hi[i + 2] + 1):
                ch[j] = FAR
                cf[j] = FAR
        for j in reversed(range(lo[i], hi[i] + 1)):
            var h: Int32
            if i == n and j == m:
                h = 0
            else:
                h = FAR
                if i < n and j < m:
                    h = nh[j + 1] + (Int32(0) if a[i] == b[j] else p.mismatch)
                if j < m:
                    h = min(h, ce[j + 1] + go)
                if i < n:
                    h = min(h, nf[j] + go)
            ch[j] = h
            ce[j] = min(h, ce[j + 1] + p.extend) if j < m else h
            cf[j] = min(h, nf[j] + p.extend) if i < n else h
        # Re-entries at row i: from above, (i, hi[i - 1] + 1 .. hi[i]); from below, (i, lo[i]).
        if i > 0 and pa[i] < FAR // 2:
            for j in range(hi[i - 1] + 1, hi[i] + 1):
                var g = min(ce[j], cf[j])
                var dd = Int32(j - i) * p.extend
                bound = min(bound, max(pa[i] + dd, pb[i] - dd) + g)
            if lo[i] > 0:
                var j = lo[i]
                var g = min(ce[j], cf[j])
                var dd = Int32(j - i) * p.extend
                bound = min(bound, max(pa[i] + dd, pb[i] - dd) + g)
        swap(nh, ch)
        swap(nf, cf)
    return bound


def _corridor(
    a: List[UInt8], b: List[UInt8], anchors: List[Anchor], p: Prices, best: Int32, prove: Bool
) raises -> String:
    """Corridors of half-width B = 27, 54, ... 1728 around the anchors' path: each one's best alignment over the
    optimum. With `prove`, it goes on doubling until the bound (`_band_bound`) proves the corridor's best optimal,
    and gives that B and the cells computed by then (forward and backward, all widths tried) against the full DP's.
    """
    var n = len(a)
    var m = len(b)
    var jlo = List[Int]()
    var jhi = List[Int]()
    _path_rows(n, m, anchors, jlo, jhi)
    var lo = List[Int](length=n + 1, fill=0)
    var hi = List[Int](length=n + 1, fill=0)
    var alpha = List[Int32]()
    var beta = List[Int32]()
    var band = 27
    var cells = 0
    var text = String("excess")
    var proved = -1
    while True:
        var width = 0
        var full = True
        for i in range(n + 1):
            lo[i] = max(0, jlo[i] - band)
            hi[i] = min(m, jhi[i] + band)
            width += hi[i] - lo[i] + 1
            if lo[i] > 0 or hi[i] < m:
                full = False
        var cost = _band_forward(a, b, lo, hi, p, alpha, beta)
        cells += width
        if band <= 1728:
            text += " " + String(band) + ":" + String(Int(cost - best))
        if prove and proved < 0:
            cells += width
            if cost <= _band_bound(a, b, lo, hi, p, alpha, beta):
                proved = band
        if full or (band >= 1728 and (not prove or proved >= 0)):
            break
        band *= 2
    if prove:
        text += (
            " | proved at "
            + String(proved)
            + ", cells "
            + String(Int(100.0 * Float64(cells) / Float64((n + 1) * (m + 1))))
            + "% of the full DP"
        )
    return text


def _segment(
    a: List[UInt8], b: List[UInt8], i0: Int, i1: Int, j0: Int, j1: Int, p: Prices, mut partner: List[Int]
) raises -> Int:
    """The best alignment of a[i0:i1] to b[j0:j1] under `p` (Gotoh, traced back); records a's letters' partners."""
    var n = i1 - i0
    var m = j1 - j0
    if n == 0 and m == 0:
        return 0
    if n == 0:
        return Int(p.open + p.extend * Int32(m))
    if m == 0:
        return Int(p.open + p.extend * Int32(n))
    var w = m + 1
    var go = p.open + p.extend
    # Per cell: the state each of D, E and F came from (0 D, 1 E, 2 F), two bits each.
    var moves = List[UInt8](length=(n + 1) * w, fill=0)
    var pd = List[Int32](length=w, fill=FAR)
    var pe = List[Int32](length=w, fill=FAR)
    var pf = List[Int32](length=w, fill=FAR)
    var cd = List[Int32](length=w, fill=FAR)
    var ce = List[Int32](length=w, fill=FAR)
    var cf = List[Int32](length=w, fill=FAR)
    pd[0] = 0
    for j in range(1, w):
        pe[j] = p.open + p.extend * Int32(j)
        moves[j] = UInt8(1 << 2) if j > 1 else UInt8(0)
    for i in range(1, n + 1):
        var x = a[i0 + i - 1]
        cd[0] = FAR
        ce[0] = FAR
        cf[0] = p.open + p.extend * Int32(i)
        moves[i * w] = UInt8(2 << 4) if i > 1 else UInt8(0)
        for j in range(1, w):
            var dh = UInt8(0)
            var hd = pd[j - 1]
            if pe[j - 1] < hd:
                hd = pe[j - 1]
                dh = 1
            if pf[j - 1] < hd:
                hd = pf[j - 1]
                dh = 2
            cd[j] = hd + (Int32(0) if x == b[j0 + j - 1] else p.mismatch)
            var eh = UInt8(1)
            var ev = ce[j - 1] + p.extend
            if cd[j - 1] + go < ev:
                ev = cd[j - 1] + go
                eh = 0
            if cf[j - 1] + go < ev:
                ev = cf[j - 1] + go
                eh = 2
            ce[j] = ev
            var fh = UInt8(2)
            var fv = pf[j] + p.extend
            if pd[j] + go < fv:
                fv = pd[j] + go
                fh = 0
            if pe[j] + go < fv:
                fv = pe[j] + go
                fh = 1
            cf[j] = fv
            moves[i * w + j] = dh | (eh << 2) | (fh << 4)
        swap(pd, cd)
        swap(pe, ce)
        swap(pf, cf)
    var state = 0
    var cost = pd[m]
    if pe[m] < cost:
        cost = pe[m]
        state = 1
    if pf[m] < cost:
        cost = pf[m]
        state = 2
    var i = n
    var j = m
    while i > 0 or j > 0:
        var came = Int((moves[i * w + j] >> UInt8(2 * state)) & 3)
        if state == 0:
            partner[i0 + i - 1] = j0 + j - 1
            i -= 1
            j -= 1
        elif state == 1:
            j -= 1
        else:
            i -= 1
        state = came
    return Int(cost)


def _chain_report(
    a: List[UInt8], b: List[UInt8], anchors: List[Anchor], t: Truth, p: Prices, prove: Bool = False
) raises -> String:
    """A chain's anchors (those on an optimal alignment), letters, forced letters it pairs right, the excess of the
    alignment through it (each gap aligned optimally) and its corridors (`_corridor`)."""
    var total = 0
    var pa = 0
    var pb = 0
    var zero = 0
    var letters = 0
    var right = 0
    var partner = List[Int](length=t.n, fill=-1)
    for x in anchors:
        total += _segment(a, b, pa, x.ref_pos, pb, x.query_pos, p, partner)
        if t.excess(x.ref_pos, x.query_pos, x.length) == 0:
            zero += 1
        letters += x.length
        for k in range(x.length):
            partner[x.ref_pos + k] = x.query_pos + k
        pa = x.ref_pos + x.length
        pb = x.query_pos + x.length
    total += _segment(a, b, pa, t.n, pb, t.m, p, partner)
    for i in range(t.n):
        if t.forced[i] and partner[i] == t.lo_j[i]:
            right += 1
    return (
        String(len(anchors))
        + " ("
        + String(zero)
        + ") | "
        + String(letters)
        + " | "
        + String(right)
        + " of "
        + String(t.forced_letters)
        + " | "
        + String(total - Int(t.best))
        + " | corridor "
        + _corridor(a, b, anchors, p, t.best, prove)
    )


def _trace_plegma(a: List[UInt8], b: List[UInt8], t: Truth, p: Prices, ambiguity: Float64) raises:
    """Every anchor `split_anchors` places (stage `flank` when the flanks alone place it too, else `split`): its
    excess; its letters an optimal alignment pairs exactly as it does, and the wrong ones at its start and end; the
    excess of its middle run of right letters alone; and the optimal diagonals at its start, middle and end."""
    var flanks = split_anchors(Span(a), Span(b), max_depth=0, ambiguity=ambiguity)
    var placed = split_anchors(Span(a), Span(b), ambiguity=ambiguity)
    var bad = List[Anchor]()
    for x in placed:
        var stage = "split"
        for f in flanks:
            if f == x:
                stage = "flank"
                break
        var good = List[Bool](length=x.length, fill=False)
        var right = 0
        for k in range(x.length):
            var i = x.ref_pos + k
            var j = x.query_pos + k
            if t.h(i, j) + t.hb[(i + 1) * t.w + j + 1] == t.best:
                good[k] = True
                right += 1
        var lead = 0
        while lead < x.length and not good[lead]:
            lead += 1
        var trail = 0
        while trail < x.length - lead and not good[x.length - 1 - trail]:
            trail += 1
        var core = x.length - lead - trail
        var trimmed = t.excess(x.ref_pos + lead, x.query_pos + lead, core) if core > 0 else -1
        var mid = x.ref_pos + x.length // 2
        var last = x.ref_pos + x.length
        if t.excess(x.ref_pos, x.query_pos, x.length) > 0:
            bad.append(x)
        print(
            "anchor",
            stage,
            "a",
            x.ref_pos,
            "b",
            x.query_pos,
            "len",
            x.length,
            "excess",
            t.excess(x.ref_pos, x.query_pos, x.length),
            "right",
            right,
            "lead",
            lead,
            "trail",
            trail,
            "trimmed",
            trimmed,
            "diag",
            x.query_pos - x.ref_pos,
            "optimal",
            t.lo_j[x.ref_pos] - x.ref_pos,
            t.hi_j[x.ref_pos] - x.ref_pos,
            t.lo_j[mid] - mid,
            t.hi_j[mid] - mid,
            t.lo_j[last] - last,
            t.hi_j[last] - last,
            "of",
            len(a),
            len(b),
        )
    # The alignment through all of them (each gap aligned optimally) over the optimum.
    var through = List[Int](length=t.n, fill=-1)
    var total = 0
    var pa = 0
    var pb = 0
    for x in placed:
        total += _segment(a, b, pa, x.ref_pos, pb, x.query_pos, p, through)
        pa = x.ref_pos + x.length
        pb = x.query_pos + x.length
    total += _segment(a, b, pa, t.n, pb, t.m, p, through)
    print("joint", total - Int(t.best), "optimum", t.best)
    if len(bad) == 0:
        return
    # One optimal alignment's gaps (a's deletions and b's insertions) within 1000 letters of a of a wrong anchor.
    var partner = List[Int](length=t.n, fill=-1)
    _ = _segment(a, b, 0, t.n, 0, t.m, p, partner)
    var i = 0
    var j = 0
    var near = List[Bool](length=t.n + 1, fill=False)
    for x in bad:
        for k in range(max(0, x.ref_pos - 1000), min(t.n, x.ref_pos + x.length + 1000) + 1):
            near[k] = True
    var text = String("optimal gaps near wrong anchors (a position, b position, size; + b letters, - a letters):")
    while i < t.n:
        if partner[i] < 0:
            var start = i
            while i < t.n and partner[i] < 0:
                i += 1
            if near[start]:
                text += " " + String(start) + "/" + String(j) + ":-" + String(i - start)
            continue
        if partner[i] > j:
            if near[i]:
                text += " " + String(i) + "/" + String(j) + ":+" + String(partner[i] - j)
        j = partner[i] + 1
        i += 1
    if t.m > j and near[t.n]:
        text += " " + String(t.n) + "/" + String(j) + ":+" + String(t.m - j)
    print(text)


def _induced_cost(a: List[UInt8], b: List[UInt8], on: List[Int], p: Prices) -> Int:
    """The cost of the alignment of a and b that pairs letter j of b with letter on[j] of a (-1: none; increasing
    otherwise), as an MSA would spell it: mismatches, and between pairs a gap in each sequence of its letters."""
    var cost = 0
    var pi = -1
    var pj = -1
    for j in range(len(b) + 1):
        var i = len(a) if j == len(b) else on[j]
        if j < len(b) and (i < 0 or i <= pi):
            continue
        var ga = i - pi - 1
        var gb = j - pj - 1
        if ga > 0:
            cost += Int(p.open + p.extend * Int32(ga))
        if gb > 0:
            cost += Int(p.open + p.extend * Int32(gb))
        if j < len(b) and a[i] != b[j]:
            cost += Int(p.mismatch)
        pi = i
        pj = j
    return cost


def _seed_audit(args: List[String]) raises:
    """`--seed-audit <dump> [--prices x,o,e]`: the anchors `plegma -S` used (PLEGMA_SEED_DUMP: per sequence, the
    consensus path it was anchored to, its letters and its anchors) against the sequence's optimal alignments to
    that path: per sequence, its anchors off an optimum and the excess of the alignment through all of them."""
    var prices = Prices(4, 6, 2)
    for k in range(len(args)):
        if args[k] == "--prices" and k + 1 < len(args):
            var parts = args[k + 1].split(",")
            prices = Prices(Int32(Int(String(parts[0]))), Int32(Int(String(parts[1]))), Int32(Int(String(parts[2]))))
    var text = String()
    with open(args[2], "r") as f:
        text = f.read()
    var lines = text.split("\n")
    var record = 0
    var total = 0
    var made_total = 0
    var optima = 0
    var k = 0
    while k + 5 < len(lines) + 1 and k < len(lines):
        if not lines[k].startswith(">path"):
            k += 1
            continue
        var a = List[UInt8](String(lines[k + 1]).as_bytes())
        var b = List[UInt8](String(lines[k + 3]).as_bytes())
        var anchors = List[Anchor]()
        var fields = String(lines[k + 4]).split(" ")
        for u in range(1, len(fields)):
            var xyz = fields[u].split(",")
            if len(xyz) == 3:
                anchors.append(Anchor(Int(String(xyz[0])), Int(String(xyz[1])), Int(String(xyz[2]))))
        k += 5
        # The alignment plegma made (`>aligned`: each letter's place on the path) and its cost.
        var made = -1
        if k < len(lines) and lines[k].startswith(">aligned"):
            var places = String(lines[k]).split(" ")
            var on = List[Int]()
            for u in range(1, len(places)):
                on.append(Int(String(places[u])))
            k += 1
            made = _induced_cost(a, b, on, prices)
        var t = Truth(a, b, prices)
        var off = 0
        var off_text = String()
        for x in anchors:
            var e = t.excess(x.ref_pos, x.query_pos, x.length)
            if e > 0:
                off += 1
                var mid = x.ref_pos + x.length // 2
                off_text += (
                    " ["
                    + String(x.ref_pos)
                    + ","
                    + String(x.query_pos)
                    + ","
                    + String(x.length)
                    + " ex "
                    + String(e)
                    + " diag "
                    + String(x.query_pos - x.ref_pos)
                    + " opt "
                    + String(t.lo_j[mid] - mid)
                    + ".."
                    + String(t.hi_j[mid] - mid)
                    + "]"
                )
        var through = List[Int](length=t.n, fill=-1)
        var cost = 0
        var pa = 0
        var pb = 0
        for x in anchors:
            cost += _segment(a, b, pa, x.ref_pos, pb, x.query_pos, prices, through)
            pa = x.ref_pos + x.length
            pb = x.query_pos + x.length
        cost += _segment(a, b, pa, t.n, pb, t.m, prices, through)
        var joint = cost - Int(t.best)
        total += joint
        optima += Int(t.best)
        if made >= 0:
            made_total += made - Int(t.best)
        print(
            "sequence",
            record,
            "path",
            t.n,
            "letters",
            t.m,
            "anchors",
            len(anchors),
            "off",
            off,
            "joint excess",
            joint,
            "made excess",
            made - Int(t.best) if made >= 0 else -1,
            "optimum",
            t.best,
            off_text,
        )
        record += 1
    print("total joint excess", total, "made excess", made_total, "of optima", optima)


def _coopt(args: List[String]) raises:
    """`--coopt <pairs.tsv> [--two-piece x,o1,e1,o2,e2]`: per line `a b row_a row_b` (two sequences and an alignment
    of them as two gapped rows), under two-piece gap costs (a mismatch x, a gap of k letters min(o1 + e1 k, o2 + e2
    k); by default pgrealign's abPOA scores as a cost, 48 and min(32 + 24k, 192 + 10k)): the optimum, the alignment's
    cost, its aligned letter pairs and how many of them some optimal alignment aligns (forward and backward DP).
    Prints `best cost pairs on_optimum flags` per line, `flags` one character per letter of a: `1` aligned as some
    optimal alignment aligns it, `0` aligned as none does, `.` not aligned."""
    var x = Int32(48)
    var o1 = Int32(32)
    var e1 = Int32(24)
    var o2 = Int32(192)
    var e2 = Int32(10)
    for k in range(len(args)):
        if args[k] == "--two-piece" and k + 1 < len(args):
            var parts = args[k + 1].split(",")
            x = Int32(Int(String(parts[0])))
            o1 = Int32(Int(String(parts[1])))
            e1 = Int32(Int(String(parts[2])))
            o2 = Int32(Int(String(parts[3])))
            e2 = Int32(Int(String(parts[4])))
    var text = String()
    with open(args[2], "r") as f:
        text = f.read()
    for line in text.split("\n"):
        var fields = String(line).split("\t")
        if len(fields) < 4:
            continue
        var a = List[UInt8](String(fields[0]).as_bytes())
        var b = List[UInt8](String(fields[1]).as_bytes())
        var ra = List[UInt8](String(fields[2]).as_bytes())
        var rb = List[UInt8](String(fields[3]).as_bytes())
        var n = len(a)
        var m = len(b)
        # The alignment's pairs (each letter of a: its partner in b or -1) and its cost.
        var partner = List[Int](length=n, fill=-1)
        var cost = 0
        var i = 0
        var j = 0
        var run = 0
        var side = 0
        for c in range(len(ra)):
            var gap_a = ra[c] == UInt8(ord("-"))
            var gap_b = rb[c] == UInt8(ord("-"))
            if gap_a and gap_b:
                continue
            var now = 1 if gap_a else (2 if gap_b else 0)
            if now != side and run > 0:
                cost += Int(min(o1 + e1 * Int32(run), o2 + e2 * Int32(run)))
                run = 0
            side = now
            if now != 0:
                run += 1
            else:
                partner[i] = j
                if a[i] != b[j]:
                    cost += Int(x)
            if not gap_a:
                i += 1
            if not gap_b:
                j += 1
        if run > 0:
            cost += Int(min(o1 + e1 * Int32(run), o2 + e2 * Int32(run)))
        # Forward: H at each pair's start cell (i, partner[i]); rows rolling, D E1 E2 F1 F2.
        var w = m + 1
        var h_prev = List[Int32](length=w, fill=FAR)
        var h_row = List[Int32](length=w, fill=FAR)
        var f1 = List[Int32](length=w, fill=FAR)
        var f2 = List[Int32](length=w, fill=FAR)
        var at_start = List[Int32](length=n, fill=FAR)
        h_prev[0] = 0
        for jj in range(1, w):
            h_prev[jj] = min(o1 + e1 * Int32(jj), o2 + e2 * Int32(jj))
        if n > 0 and partner[0] >= 0:
            at_start[0] = h_prev[partner[0]]
        for ii in range(1, n + 1):
            var ca = a[ii - 1]
            var g1 = FAR
            var g2 = FAR
            for jj in range(w):
                f1[jj] = min(h_prev[jj] + o1 + e1, f1[jj] + e1)
                f2[jj] = min(h_prev[jj] + o2 + e2, f2[jj] + e2)
                var h = min(f1[jj], f2[jj])
                if jj > 0:
                    g1 = min(h_row[jj - 1] + o1 + e1, g1 + e1)
                    g2 = min(h_row[jj - 1] + o2 + e2, g2 + e2)
                    var d = h_prev[jj - 1] + (Int32(0) if ca == b[jj - 1] else x)
                    h = min(h, min(d, min(g1, g2)))
                h_row[jj] = h
            if ii < n and partner[ii] >= 0:
                at_start[ii] = h_row[partner[ii]]
            swap(h_prev, h_row)
        var best = h_prev[m]
        # Backward: H' at each pair's end cell (i + 1, partner[i] + 1).
        var hb_next = List[Int32](length=w, fill=FAR)
        var hb_row = List[Int32](length=w, fill=FAR)
        var fb1 = List[Int32](length=w, fill=FAR)
        var fb2 = List[Int32](length=w, fill=FAR)
        var at_end = List[Int32](length=n, fill=FAR)
        for ii in reversed(range(n + 1)):
            var eb1 = FAR
            var eb2 = FAR
            for jj in reversed(range(w)):
                var h: Int32
                if ii == n and jj == m:
                    h = 0
                else:
                    h = FAR
                    if ii < n and jj < m:
                        h = hb_next[jj + 1] + (Int32(0) if a[ii] == b[jj] else x)
                    if jj < m:
                        h = min(h, min(eb1 + o1 + e1, eb2 + o2 + e2))
                    if ii < n:
                        h = min(h, min(fb1[jj] + o1 + e1, fb2[jj] + o2 + e2))
                hb_row[jj] = h
                eb1 = min(h, eb1 + e1) if jj < m else h
                eb2 = min(h, eb2 + e2) if jj < m else h
                fb1[jj] = min(h, fb1[jj] + e1) if ii < n else h
                fb2[jj] = min(h, fb2[jj] + e2) if ii < n else h
            swap(hb_next, hb_row)
            # hb_next now holds row ii: the end cells of letters ii - 1.
            if ii >= 1 and partner[ii - 1] >= 0:
                at_end[ii - 1] = hb_next[partner[ii - 1] + 1]
        var pairs = 0
        var on = 0
        var flags = List[UInt8](length=n, fill=UInt8(ord(".")))
        for ii in range(n):
            if partner[ii] < 0:
                continue
            pairs += 1
            var sub = Int32(0) if a[ii] == b[partner[ii]] else x
            if at_start[ii] + sub + at_end[ii] == best:
                on += 1
                flags[ii] = UInt8(ord("1"))
            else:
                flags[ii] = UInt8(ord("0"))
        if n > 0:
            print(best, cost, pairs, on, StringSlice(unsafe_from_utf8=Span(flags)))
        else:
            print(best, cost, pairs, on, "-")


def _band_of(cigar: String, n: Int, m: Int, mut lo: List[Int], mut hi: List[Int]) raises:
    """Widen [lo[i], hi[i]] (rows 0..n of the DP over a's letters) to every column the alignment `cigar` of a (the
    reference) and b visits in row i."""
    var runs = cigar_runs(cigar)
    var i = 0
    var j = 0
    lo[0] = min(lo[0], 0)
    hi[0] = max(hi[0], 0)
    for x in range(len(runs.ops)):
        var op = runs.ops[x]
        for _ in range(runs.counts[x]):
            if op == DELETION:
                i += 1
            elif op == INSERTION:
                j += 1
            else:
                i += 1
                j += 1
            lo[i] = min(lo[i], j)
            hi[i] = max(hi[i], j)


def _coopt_set(args: List[String]) raises:
    """`--coopt-set <pairs.tsv> [--two-piece x,o1,e1,o2,e2]`: per line `a b` (two sequences), under two-piece costs
    (as --coopt), every letter pair some optimal alignment aligns. The forward and backward DP run only between the
    pair's leftmost and rightmost optimal alignments (karan, Ties.LEFT and Ties.RIGHT) and MARGIN columns either
    side, MARGIN being 64 letters and the extremes' longest gap (an optimal alignment that splits it may shift by up
    to its length between them); the band is widened (the margin four times over) and the pair done again until its optimum is karan's and a
    band four times wider finds the same co-optimal letter pairs (with two-piece gaps some optimal alignments stray
    past the two extremes, and one that leaves a band is seen by neither pass inside it). ANCHORS_COOPT_MARGIN sets
    the first margin (a margin past both lengths is the full DP). Prints per line `best n partners`, `partners` one token per letter of a: `.` for none, else
    its partners in b as intervals `j` or `j1-j2` joined by `:`."""
    var x = Int32(48)
    var o1 = Int32(32)
    var e1 = Int32(24)
    var o2 = Int32(192)
    var e2 = Int32(10)
    for k in range(len(args)):
        if args[k] == "--two-piece" and k + 1 < len(args):
            var parts = args[k + 1].split(",")
            x = Int32(Int(String(parts[0])))
            o1 = Int32(Int(String(parts[1])))
            e1 = Int32(Int(String(parts[2])))
            o2 = Int32(Int(String(parts[3])))
            e2 = Int32(Int(String(parts[4])))
    var costs = Costs.two_piece(Int(x), Int(o1), Int(e1), Int(o2), Int(e2))
    var text = String()
    with open(args[2], "r") as f:
        text = f.read()
    var aligner = Aligner()
    for line in text.split("\n"):
        var fields = String(line).split("\t")
        if len(fields) < 2:
            continue
        var a = List[UInt8](String(fields[0]).as_bytes())
        var b = List[UInt8](String(fields[1]).as_bytes())
        var n = len(a)
        var m = len(b)
        var left = aligner.align(Span(a), Span(b), costs, Mode.GLOBAL, ties=Ties.LEFT, max_memory=4 << 30)
        var right = aligner.align(Span(a), Span(b), costs, Mode.GLOBAL, ties=Ties.RIGHT, max_memory=4 << 30)
        var optimum = Int32(left.cost)
        var path_lo = List[Int](length=n + 1, fill=m + 1)
        var path_hi = List[Int](length=n + 1, fill=-1)
        _band_of(left.cigar, n, m, path_lo, path_hi)
        _band_of(right.cigar, n, m, path_lo, path_hi)
        # An optimal alignment may split a gap of the extremes and shift by up to its length between them.
        var longest = 0
        for side in range(2):
            var runs = cigar_runs(left.cigar if side == 0 else right.cigar)
            for r in range(len(runs.ops)):
                if runs.ops[r] == DELETION or runs.ops[r] == INSERTION:
                    longest = max(longest, runs.counts[r])
        var margin = 64 + longest
        var last_cells = -1
        var forced = getenv("ANCHORS_COOPT_MARGIN")
        if forced != "":
            margin = Int(forced)
        while True:
            # The band: the two paths' columns per row and `margin` either side, kept monotone so each row's band
            # reaches the next's.
            var lo = List[Int](length=n + 1, fill=0)
            var hi = List[Int](length=n + 1, fill=m)
            for i in range(n + 1):
                lo[i] = max(0, path_lo[i] - margin)
                hi[i] = min(m, path_hi[i] + margin)
            for i in range(1, n + 1):
                if lo[i] < lo[i - 1]:
                    lo[i] = lo[i - 1]
            for i in reversed(range(n)):
                if hi[i] > hi[i + 1]:
                    hi[i] = hi[i + 1]
            var start = List[Int](length=n + 2, fill=0)
            for i in range(n + 1):
                start[i + 1] = start[i] + (hi[i] - lo[i] + 1)
            # Forward, every row's H kept (the band's cells); F1 F2 (gaps in b, down a column) rolling.
            var fwd = List[Int32](length=start[n + 1], fill=FAR)
            var f1 = List[Int32](length=m + 1, fill=FAR)
            var f2 = List[Int32](length=m + 1, fill=FAR)
            for jj in range(lo[0], hi[0] + 1):
                fwd[start[0] + jj - lo[0]] = Int32(0) if jj == 0 else min(o1 + e1 * Int32(jj), o2 + e2 * Int32(jj))
            for ii in range(1, n + 1):
                var ca = a[ii - 1]
                var g1 = FAR
                var g2 = FAR
                var row = start[ii] - lo[ii]
                var prev = start[ii - 1] - lo[ii - 1]
                for jj in range(lo[ii], hi[ii] + 1):
                    var up = fwd[prev + jj] if jj >= lo[ii - 1] and jj <= hi[ii - 1] else FAR
                    var diag = fwd[prev + jj - 1] if jj - 1 >= lo[ii - 1] and jj - 1 <= hi[ii - 1] else FAR
                    f1[jj] = min(up + o1 + e1, f1[jj] + e1)
                    f2[jj] = min(up + o2 + e2, f2[jj] + e2)
                    var h = min(f1[jj], f2[jj])
                    if jj > 0:
                        var west = fwd[row + jj - 1] if jj - 1 >= lo[ii] else FAR
                        g1 = min(west + o1 + e1, g1 + e1)
                        g2 = min(west + o2 + e2, g2 + e2)
                        if diag < FAR:
                            h = min(h, diag + (Int32(0) if ca == b[jj - 1] else x))
                        h = min(h, min(g1, g2))
                    fwd[row + jj] = h
                # A column leaving the band upward resets its gap states (they are FAR above the band).
                if ii < n:
                    for jj in range(lo[ii], lo[ii + 1]):
                        f1[jj] = FAR
                        f2[jj] = FAR
            var best = fwd[start[n] + m - lo[n]] if lo[n] <= m and m <= hi[n] else FAR
            if best != optimum and margin < n + m:
                margin *= 4
                continue
            # Backward, rolling; each match cell (ii, jj) checked against the forward H at (ii - 1, jj - 1).
            var hb_next = List[Int32](length=m + 2, fill=FAR)
            var hb_row = List[Int32](length=m + 2, fill=FAR)
            var fb1 = List[Int32](length=m + 1, fill=FAR)
            var fb2 = List[Int32](length=m + 1, fill=FAR)
            var partners = List[List[Int]]()
            for _ in range(n):
                partners.append(List[Int]())
            for ii in reversed(range(n + 1)):
                var eb1 = FAR
                var eb2 = FAR
                for jj in reversed(range(lo[ii], hi[ii] + 1)):
                    var h: Int32
                    if ii == n and jj == m:
                        h = 0
                    else:
                        h = FAR
                        if ii < n and jj < m and jj + 1 >= lo[ii + 1] and jj + 1 <= hi[ii + 1]:
                            h = hb_next[jj + 1] + (Int32(0) if a[ii] == b[jj] else x)
                        if jj < m:
                            h = min(h, min(eb1 + o1 + e1, eb2 + o2 + e2))
                        if ii < n:
                            h = min(h, min(fb1[jj] + o1 + e1, fb2[jj] + o2 + e2))
                    hb_row[jj] = h
                    eb1 = min(h, eb1 + e1) if jj < m else h
                    eb2 = min(h, eb2 + e2) if jj < m else h
                    fb1[jj] = min(h, fb1[jj] + e1) if ii < n else h
                    fb2[jj] = min(h, fb2[jj] + e2) if ii < n else h
                    if ii >= 1 and jj >= 1 and jj - 1 >= lo[ii - 1] and jj - 1 <= hi[ii - 1]:
                        var f = fwd[start[ii - 1] + jj - 1 - lo[ii - 1]]
                        if f < FAR and h < FAR and f + (Int32(0) if a[ii - 1] == b[jj - 1] else x) + h == best:
                            partners[ii - 1].append(jj - 1)
                # Columns of this row outside the next row's band read as FAR there.
                if ii > 0:
                    for jj in range(hi[ii - 1] + 1, hi[ii] + 1):
                        fb1[jj] = FAR
                        fb2[jj] = FAR
                for jj in range(lo[ii], hi[ii] + 1):
                    hb_next[jj] = hb_row[jj]
                if ii > 0 and lo[ii] > lo[ii - 1]:
                    for jj in range(lo[ii - 1], lo[ii]):
                        hb_next[jj] = FAR
            # Stable across a band four times wider: some optimal alignments stray past the two extremes (a long gap
            # split in two, its pieces before and after theirs), and one that leaves the band is seen by neither
            # pass inside it; so the set is kept once a band four times wider finds the same cells.
            var cells = 0
            for i in range(n):
                cells += len(partners[i])
            if getenv("ANCHORS_COOPT_STABLE") == "1" and margin < n + m and cells != last_cells:
                last_cells = cells
                margin *= 4
                continue
            var out = String(best, " ", n, " ")
            for i in range(n):
                if i > 0:
                    out += ","
                ref ps = partners[i]
                if len(ps) == 0:
                    out += "."
                    continue
                # Descending as the backward pass found them: intervals ascending.
                var k = len(ps) - 1
                var first = True
                while k >= 0:
                    var s0 = ps[k]
                    var s1 = s0
                    while k - 1 >= 0 and ps[k - 1] == s1 + 1:
                        k -= 1
                        s1 = ps[k]
                    if not first:
                        out += ":"
                    out += String(s0) if s0 == s1 else String(s0, "-", s1)
                    first = False
                    k -= 1
            print(out)
            break


def _brute(args: List[String]) raises:
    """`--brute <sequences.fa> <a> <b> <period> [bin] [--prices x,o,e]`: every exact match of 21 letters or more
    between sequences a and b (every MEM, each occurrence) checked against every optimal alignment of the pair under
    the prices (unit cost by default; a gap of k letters costs o + e * k), from the full forward and backward DP.

    Prints the pinches (boundaries of a every optimal alignment crosses at one letter of b: the placements the
    optimum forces) and the letters of a every optimal alignment pairs with the same letter of b; per kind (MUM,
    SMEM, other MEM) and length, how many matches lie on an optimal alignment and their forced excess (the best
    alignment through the match, over the optimum); per MUM and long SMEM, where it is, its diagonal against the
    optimal diagonals there (in periods) and its excess; colinear chains of the matches (all, MUMs and SMEMs, MUMs) under
    three gap penalties, with the excess of the alignment through each chain and the corridor around it that provably
    holds an optimum (`_corridor`); and along a, per `bin` letters (500), the optimal
    diagonals, the pinches and the matches."""
    var records = read_fasta(args[2])
    var seqs = List[List[UInt8]]()
    for r in range(len(records)):
        var letters = List[UInt8](records.seq(r))
        for k in range(len(letters)):
            if letters[k] >= UInt8(ord("a")) and letters[k] <= UInt8(ord("z")):
                letters[k] -= 32
        var dup = False
        for x in seqs:
            if x == letters:
                dup = True
                break
        if not dup:
            seqs.append(letters^)
    ref a = seqs[Int(args[3])]
    ref b = seqs[Int(args[4])]
    var period = Int(args[5])
    var bin = Int(args[6]) if len(args) > 6 and not args[6].startswith("--") else 500
    var prices = Prices(1, 0, 1)
    for k in range(len(args)):
        if args[k] == "--prices" and k + 1 < len(args):
            var parts = args[k + 1].split(",")
            prices = Prices(Int32(Int(String(parts[0]))), Int32(Int(String(parts[1]))), Int32(Int(String(parts[2]))))
    var n = len(a)
    var m = len(b)
    var t = Truth(a, b, prices)
    var best = t.best
    ref lo_j = t.lo_j
    ref hi_j = t.hi_j
    ref count = t.count
    ref forced = t.forced
    var forced_letters = t.forced_letters
    var pinches = 0
    for i in range(n + 1):
        if count[i] == 1:
            pinches += 1
    print(
        "pair",
        args[3],
        args[4],
        "lengths",
        n,
        m,
        "optimum",
        best,
        "prices",
        prices.mismatch,
        prices.open,
        prices.extend,
    )
    print("pinches", pinches, "of", n + 1, "boundaries of a; letters of a forced to one partner:", forced_letters)
    # Runs of pinches (the stretches of a whose placement the optimum fixes).
    var runs = String("pinch runs (a's letters, b's letter at the start):")
    var start = -1
    for i in range(n + 2):
        var p = i <= n and count[i] == 1
        if p and start < 0:
            start = i
        elif not p and start >= 0:
            if i - start >= 10:
                runs += " " + String(start) + "-" + String(i - 1) + "@" + String(lo_j[start])
            start = -1
    print(runs)
    for k in range(len(args)):
        if args[k] == "--anchors":
            var ambiguity = 0.0
            for u in range(len(args)):
                if args[u] == "--ambiguity" and u + 1 < len(args):
                    ambiguity = Float64(args[u + 1])
            _trace_plegma(a, b, t, prices, ambiguity)
            return
    var found = find_anchors(Span(a), Span(b), 21, 1000)
    var smem_set = smems(find_anchors(Span(a), Span(b), 21, 1000))
    comptime KINDS = 3
    comptime BUCKETS = 6
    # Per kind and length bucket: matches, on an optimal alignment (excess 0), excess 1-9, 10-99, 100+.
    var tally = List[Int](length=KINDS * BUCKETS * 5, fill=0)
    var bins = (n + bin - 1) // bin
    var bin_right = List[Int](length=bins, fill=0)
    var bin_near = List[Int](length=bins, fill=0)
    var bin_far = List[Int](length=bins, fill=0)
    var bin_mum_right = List[Int](length=bins, fill=0)
    var bin_mum = List[Int](length=bins, fill=0)
    var covered = List[Bool](length=n, fill=False)
    var kinds = ["MUM", "SMEM (not unique)", "other MEM"]
    var buckets = ["21-50", "51-100", "101-200", "201-400", "401-1000", ">1000"]
    var detail = List[String]()
    for x in found:
        var L = x.length
        var excess = t.excess(x.ref_pos, x.query_pos, L)
        var kind = 2
        if x.ref_copies == 1 and x.query_copies == 1:
            kind = 0
        else:
            for y in smem_set:
                if y == x:
                    kind = 1
                    break
        var bucket = 0 if L <= 50 else (
            1 if L <= 100 else (2 if L <= 200 else (3 if L <= 400 else (4 if L <= 1000 else 5)))
        )
        var base = (kind * BUCKETS + bucket) * 5
        tally[base] += 1
        tally[base + (1 if excess == 0 else (2 if excess < 10 else (3 if excess < 100 else 4)))] += 1
        var mid = x.ref_pos + L // 2
        var k = mid // bin
        # Diagonal against the optimal cells' diagonals at the match's middle, in periods.
        var d = x.query_pos - x.ref_pos
        var dlo = lo_j[mid] - mid
        var dhi = hi_j[mid] - mid
        var off = 0
        if d < dlo:
            off = -((dlo - d + period // 2) // period)
        elif d > dhi:
            off = (d - dhi + period // 2) // period
        if excess == 0:
            bin_right[k] += 1
            for t_ in range(L):
                covered[x.ref_pos + t_] = True
        elif off == 0:
            bin_near[k] += 1
        else:
            bin_far[k] += 1
        if kind == 0:
            bin_mum[k] += 1
            if excess == 0:
                bin_mum_right[k] += 1
        if kind == 0:
            var f = 0
            # The longest run of a's letters in the match that an optimal alignment matches exactly, one after
            # another (each matched letter's diagonal step stays on an optimum and continues the run).
            var run = 0
            var longest = 0
            for t_ in range(L):
                if forced[x.ref_pos + t_]:
                    f += 1
                run = run + 1 if t.matched(a, b, x.ref_pos + t_) else 0
                longest = max(longest, run)
            detail.append(
                String("  ")
                + kinds[kind]
                + " a "
                + String(x.ref_pos)
                + " b "
                + String(x.query_pos)
                + " len "
                + String(L)
                + " copies "
                + String(x.ref_copies)
                + "/"
                + String(x.query_copies)
                + " diag "
                + String(d)
                + " optimal "
                + String(dlo)
                + ".."
                + String(dhi)
                + " periods off "
                + String(off)
                + " excess "
                + String(excess)
                + " forced letters "
                + String(f)
                + " optimal exact run "
                + String(longest)
            )
    var cover = 0
    var opt_letters = 0
    for i in range(n):
        # A letter of a some optimal alignment matches: an optimal cell (i, j) with a[i] == b[j] whose diagonal
        # step stays optimal.
        if t.matched(a, b, i):
            opt_letters += 1
            if covered[i]:
                cover += 1
    print(
        "matches",
        len(found),
        "| a's letters some optimal alignment matches:",
        opt_letters,
        "| covered by matches with excess 0:",
        cover,
    )
    print("kind               length     all  excess 0  1-9  10-99  100+")
    for kind in range(KINDS):
        for bucket in range(BUCKETS):
            var base = (kind * BUCKETS + bucket) * 5
            if tally[base] > 0:
                print(
                    "  ",
                    kinds[kind],
                    " " * (17 - kinds[kind].byte_length()),
                    buckets[bucket],
                    " " * (9 - buckets[bucket].byte_length()),
                    tally[base],
                    tally[base + 1],
                    tally[base + 2],
                    tally[base + 3],
                    tally[base + 4],
                )
    print("MUMs:")
    for line in detail:
        print(line)
    # Colinear chains: matches in order of a, each after its predecessor in both sequences; score the letters
    # matched minus a gap penalty of |gap in a - gap in b| + slope * min(gap in a, gap in b).
    var by_start = List[List[Int]]()
    for _ in range(n):
        by_start.append(List[Int]())
    for k in range(len(found)):
        by_start[found[k].ref_pos].append(k)
    var order = List[Int]()
    for i in range(n):
        for k in by_start[i]:
            order.append(k)
    var kind_of = List[Int](length=len(found), fill=2)
    for k in range(len(found)):
        ref x = found[k]
        if x.ref_copies == 1 and x.query_copies == 1:
            kind_of[k] = 0
        else:
            for y in smem_set:
                if y == x:
                    kind_of[k] = 1
                    break
    var set_names = ["all MEMs", "MUMs + SMEMs", "MUMs"]
    var slopes = [0.0, 0.5, 1.0]
    print(
        "chains: set, slope | anchors (excess 0) | letters | forced letters paired right | excess through the chain |"
        " corridor: excess by half-width B (and the B the bound proves optimal)"
    )
    for set in range(3):
        for slope in slopes:
            var pick = List[Int]()
            for k in order:
                if kind_of[k] <= 2 - set:
                    pick.append(k)
            var score = List[Float64](length=len(pick), fill=0.0)
            var back = List[Int](length=len(pick), fill=-1)
            var top = -1
            for u in range(len(pick)):
                ref x = found[pick[u]]
                var sx = (
                    Float64(x.length)
                    - Float64(abs(x.ref_pos - x.query_pos))
                    - slope * Float64(min(x.ref_pos, x.query_pos))
                )
                var bu = -1
                for v in range(u):
                    ref y = found[pick[v]]
                    var ga = x.ref_pos - (y.ref_pos + y.length)
                    var gb = x.query_pos - (y.query_pos + y.length)
                    if ga < 0 or gb < 0:
                        continue
                    var c = score[v] + Float64(x.length) - Float64(abs(ga - gb)) - slope * Float64(min(ga, gb))
                    if c > sx:
                        sx = c
                        bu = v
                score[u] = sx
                back[u] = bu
            var top_score = -1.0e18
            for u in range(len(pick)):
                ref x = found[pick[u]]
                var ea = n - (x.ref_pos + x.length)
                var eb = m - (x.query_pos + x.length)
                var f = score[u] - Float64(abs(ea - eb)) - slope * Float64(min(ea, eb))
                if f > top_score:
                    top_score = f
                    top = u
            var chain = List[Int]()
            var u = top
            while u >= 0:
                chain.append(pick[u])
                u = back[u]
            chain.reverse()
            var anchors = List[Anchor]()
            for k in chain:
                anchors.append(found[k])
            print("  ", set_names[set], slope, "|", _chain_report(a, b, anchors, t, prices))
    # ChainX's (and llchain's) chaining: links cost the longer side of their gap, overlaps their diagonal
    # difference; its cost is the anchored edit distance (every match outside the anchors an edit).
    for set in range(2):
        var pool = List[Anchor]()
        for k in range(len(found)):
            if set == 0 or kind_of[k] == 0:
                pool.append(found[k])
        var cx = chain(pool, n, m, ChainCosts.chainx())
        print(
            "  ChainX "
            + ("all MEMs" if set == 0 else "MUMs")
            + " (anchored edit distance "
            + String(Int(cx.cost))
            + ") |",
            _chain_report(a, b, without_overlaps(cx.anchors), t, prices),
        )
    print("  plegma split_anchors |", _chain_report(a, b, split_anchors(Span(a), Span(b)), t, prices, True))
    print(
        "  plegma flanks only   |", _chain_report(a, b, split_anchors(Span(a), Span(b), max_depth=0), t, prices, True)
    )
    print("  straight corridor (no seeds) |", _corridor(a, b, List[Anchor](), prices, best, True))
    print(
        "along a (per",
        bin,
        (
            "letters): optimal diagonals (widest), pinches | matches with excess 0, off within the optimal diagonals,"
            " off by periods | MUMs right/all"
        ),
    )
    for k in range(bins):
        var dlo = Int.MAX
        var dhi = Int.MIN
        var widest = 0
        var p = 0
        for i in range(k * bin, min(n, (k + 1) * bin)):
            dlo = min(dlo, lo_j[i] - i)
            dhi = max(dhi, hi_j[i] - i)
            widest = max(widest, hi_j[i] - lo_j[i])
            if count[i] == 1:
                p += 1
        print(
            "  ",
            k * bin,
            ": diagonals",
            dlo,
            "..",
            dhi,
            "(widest",
            widest,
            ") pinches",
            p,
            "|",
            bin_right[k],
            bin_near[k],
            bin_far[k],
            "| MUMs",
            bin_mum_right[k],
            "/",
            bin_mum[k],
        )


def main() raises:
    var args = argv()
    if len(args) > 1 and String(args[1]) == "--coopt-set":
        var list = List[String]()
        for k in range(len(args)):
            list.append(String(args[k]))
        _coopt_set(list)
        return
    if len(args) > 1 and String(args[1]) == "--coopt":
        var list = List[String]()
        for k in range(len(args)):
            list.append(String(args[k]))
        _coopt(list)
        return
    if len(args) > 1 and String(args[1]) == "--seed-audit":
        var list = List[String]()
        for k in range(len(args)):
            list.append(String(args[k]))
        _seed_audit(list)
        return
    if len(args) > 1 and String(args[1]) == "--brute":
        var list = List[String]()
        for k in range(len(args)):
            list.append(String(args[k]))
        _brute(list)
        return
    if len(args) > 1 and String(args[1]) == "--compose2":
        var list = List[String]()
        for k in range(len(args)):
            list.append(String(args[k]))
        _compose2(list)
        return
    if len(args) > 1 and String(args[1]) == "--compose":
        var list = List[String]()
        for k in range(len(args)):
            list.append(String(args[k]))
        _compose(list)
        return
    if len(args) < 6:
        print("usage: anchors <sequences.fa> <period> <max pairs> <out.tsv> <plan> [<plan> ...]")
        return
    var records = read_fasta(String(args[1]))
    var period = Int(String(args[2]))
    var max_pairs = Int(String(args[3]))
    var plans = List[Plan]()
    for k in range(5, len(args)):
        plans.append(Plan(String(args[k])))
    var seqs = List[List[UInt8]]()
    for r in range(len(records)):
        var letters = List[UInt8](records.seq(r))
        for k in range(len(letters)):
            if letters[k] >= UInt8(ord("a")) and letters[k] <= UInt8(ord("z")):
                letters[k] -= 32
        var dup = False
        for x in seqs:
            if x == letters:
                dup = True
                break
        if not dup:
            seqs.append(letters^)
    # MSA plans: each distinct sequence's row.
    for pi in range(len(plans)):
        if plans[pi].kind != "msa":
            continue
        var msa = read_fasta(plans[pi].msa_path)
        var rows = List[List[UInt8]]()
        for _ in range(len(seqs)):
            rows.append(List[UInt8]())
        for r in range(len(msa)):
            var row = List[UInt8](msa.seq(r))
            var spelled = List[UInt8]()
            for k in range(len(row)):
                if row[k] >= UInt8(ord("a")) and row[k] <= UInt8(ord("z")):
                    row[k] -= 32
                if row[k] == UInt8(ord(".")):
                    row[k] = UInt8(ord("-"))
                if row[k] != UInt8(ord("-")):
                    spelled.append(row[k])
            for d in range(len(seqs)):
                if len(rows[d]) == 0 and seqs[d] == spelled:
                    rows[d] = row.copy()
        for d in range(len(seqs)):
            if len(rows[d]) == 0:
                raise Error(t"{plans[pi].msa_path}: no row spells sequence {d}")
        plans[pi].rows = rows^
    # The region's unit: the longest sequence's array at the period, its consensus.
    var unit = List[UInt8]()
    if period > 0:
        var longest = 0
        for k in range(1, len(seqs)):
            if len(seqs[k]) > len(seqs[longest]):
                longest = k
        var top = repeat_at(seqs[longest], period)
        if top.period > 0:
            unit = consensus_unit(seqs[longest], top.start, top.end, period)
    var pairs = List[Tuple[Int, Int]]()
    for x in range(len(seqs)):
        for y in range(x + 1, len(seqs)):
            pairs.append((x, y) if len(seqs[x]) >= len(seqs[y]) else (y, x))
    var chosen = List[Tuple[Int, Int]]()
    var step = max(1, len(pairs) // max(1, max_pairs))
    var t = 0
    while t < len(pairs) and len(chosen) < max_pairs:
        chosen.append(pairs[t])
        t += step
    var n_plans = len(plans)
    comptime DEPTHS = 8
    # Per plan and depth (the last one gathering deeper ones): anchors, off anchors, letters, off letters.
    var tally = List[Int](length=n_plans * DEPTHS * 4, fill=0)
    var by_class = List[Int](length=n_plans * 3, fill=0)
    var forced = List[Int](length=n_plans, fill=0)
    var optimum = 0
    var aligner = Aligner()
    var text = String()
    text.write(
        "plan\ta\tb\tlen_a\tlen_b\tdepth\tlevel\tk\tlength\tref_copies\tquery_copies\twin_ref\twin_query\tin_repeat",
        "\tsiblings\tdiag_start\tdiag_end\tletters\toff\toff_near\toff_copies\toff_other\tref_pos\tquery_pos\n",
    )
    var done = 0
    for pq in chosen:
        ref a = seqs[pq[0]]
        ref b = seqs[pq[1]]
        var n = len(a)
        var m = len(b)
        if (n + 1) * (m + 1) > 150_000_000:
            continue
        done += 1
        var w = m + 1
        var fwd = List[Int32](length=(n + 1) * w, fill=0)
        for j in range(w):
            fwd[j] = Int32(j)
        for i in range(1, n + 1):
            fwd[i * w] = Int32(i)
            var x = a[i - 1]
            for j in range(1, w):
                var d = fwd[(i - 1) * w + j - 1] + (Int32(0) if x == b[j - 1] else Int32(1))
                fwd[i * w + j] = min(d, min(fwd[(i - 1) * w + j] + 1, fwd[i * w + j - 1] + 1))
        var best = fwd[n * w + m]
        optimum += Int(best)
        for pi in range(n_plans):
            ref plan = plans[pi]
            var repeats = RepeatMask(Span(a), plan.mask_k, plan.mask_distance)
            var traced = List[Traced]()
            if plan.kind == "seeds":
                _trace(a, b, 0, n, 0, m, plan, 0, 0, repeats, traced)
            elif plan.kind == "msa":
                _stacked(plan.rows[pq[0]], plan.rows[pq[1]], traced)
            elif plan.kind == "libsplit":
                for x in split_anchors(Span(a), Span(b)):
                    traced.append(Traced(x.ref_pos, x.query_pos, x.length, 0, 0, 1, 1, n, m, False, 0, 0, 0))
            elif plan.kind == "split_arrays":
                # The arrays' boundaries as seeds (each sequence's array at the period, as longdust or TRF would
                # find it): the left flanks, the arrays and the right flanks split each on its own.
                var ra = repeat_at(a, period) if period > 0 else Repeat(0, 0, 0, 0.0)
                var rb = repeat_at(b, period) if period > 0 else Repeat(0, 0, 0, 0.0)
                if ra.period == 0 or rb.period == 0:
                    _split(a, b, 0, n, 0, m, plan.min_run, plan.ends, plan.len_d, plan.max_depth, 0, traced)
                else:
                    _split(
                        a, b, 0, ra.start, 0, rb.start, plan.min_run, plan.ends, plan.len_d, plan.max_depth, 1, traced
                    )
                    _split(
                        a,
                        b,
                        ra.start,
                        ra.end,
                        rb.start,
                        rb.end,
                        plan.min_run,
                        plan.ends,
                        plan.len_d,
                        plan.max_depth,
                        1,
                        traced,
                    )
                    _split(a, b, ra.end, n, rb.end, m, plan.min_run, plan.ends, plan.len_d, plan.max_depth, 1, traced)
            elif plan.kind == "split":
                if plan.mask_k > 0:
                    # The flanks first, then splitting every gap they leave.
                    var flanks = List[Traced]()
                    _trace(a, b, 0, n, 0, m, plan, 0, 0, repeats, flanks)
                    var pa = 0
                    var pb = 0
                    for f in flanks:
                        _split(
                            a,
                            b,
                            pa,
                            f.ref_pos,
                            pb,
                            f.query_pos,
                            plan.min_run,
                            plan.ends,
                            plan.len_d,
                            plan.max_depth,
                            1,
                            traced,
                        )
                        traced.append(f.copy())
                        pa = f.ref_pos + f.length
                        pb = f.query_pos + f.length
                    _split(a, b, pa, n, pb, m, plan.min_run, plan.ends, plan.len_d, plan.max_depth, 1, traced)
                else:
                    _split(a, b, 0, n, 0, m, plan.min_run, plan.ends, plan.len_d, plan.max_depth, 0, traced)
            elif plan.kind == "greedy":
                _greedy(a, b, period * plan.min_run if period > 0 else 21 * plan.min_run, plan.stay, traced)
            elif plan.kind == "copymums":
                _trace(a, b, 0, n, 0, m, plan, 0, 0, repeats, traced)
                if period > 0:
                    _copy_mums(a, b, period, unit, traced)
            elif plan.kind == "copies":
                _trace(a, b, 0, n, 0, m, plan, 0, 0, repeats, traced)
                if period > 0:
                    _copy_anchors(a, b, period, unit, plan.min_run, traced)
            else:
                _aligned(a, b, plan.kind == "edit", traced)
            # Each anchored letter pair, by partner letter: (query letter, anchor).
            var at = List[List[Int]]()
            var owner = List[List[Int]]()
            for _ in range(n):
                at.append(List[Int]())
                owner.append(List[Int]())
            var letters = List[Int](length=len(traced), fill=0)
            var off = List[Int](length=len(traced), fill=0)
            var near = List[Int](length=len(traced), fill=0)
            var copies = List[Int](length=len(traced), fill=0)
            var other = List[Int](length=len(traced), fill=0)
            for k in range(len(traced)):
                ref an = traced[k]
                var lo = min(plan.trim, an.length // 2)
                for q in range(lo, an.length - lo):
                    at[an.ref_pos + q].append(an.query_pos + q)
                    owner[an.ref_pos + q].append(k)
            var next_row = List[Int32](length=w, fill=0)
            var row = List[Int32](length=w, fill=0)
            for j in range(w):
                next_row[j] = Int32(m - j)
            for i in reversed(range(n)):
                row[m] = Int32(n - i)
                var x = a[i]
                for j in reversed(range(m)):
                    var d = next_row[j + 1] + (Int32(0) if x == b[j] else Int32(1))
                    row[j] = min(d, min(next_row[j] + 1, row[j + 1] + 1))
                for e in range(len(at[i])):
                    var j = at[i][e]
                    var k = owner[i][e]
                    letters[k] += 1
                    if fwd[i * w + j] + (Int32(0) if x == b[j] else Int32(1)) + next_row[j + 1] == best:
                        continue
                    off[k] += 1
                    var offset = m
                    for dj in range(1, m):
                        var hit = False
                        for jj in [j - dj, j + dj]:
                            if jj >= 0 and jj < m:
                                if fwd[i * w + jj] + (Int32(0) if x == b[jj] else Int32(1)) + next_row[jj + 1] == best:
                                    hit = True
                        if hit:
                            offset = dj
                            break
                    if offset <= 2:
                        near[k] += 1
                    elif period > 0 and abs(offset - ((offset + period // 2) // period) * period) <= 2:
                        copies[k] += 1
                    else:
                        other[k] += 1
                next_row = row.copy()
            # What the anchors force: the pieces between consecutive anchored letters (kept while they rise on both
            # sequences, by partner letter) aligned optimally.
            var keys = List[UInt64]()
            for k in range(len(traced)):
                keys.append((UInt64(traced[k].ref_pos) << 32) | UInt64(k))
            sort(keys)
            var cost = 0
            var pi_at = 0
            var pj_at = 0
            for key in keys:
                ref an = traced[Int(key & 0xFFFFFFFF)]
                var lo = min(plan.trim, an.length // 2)
                var i0 = an.ref_pos + lo
                var j0 = an.query_pos + lo
                var len_kept = an.length - 2 * lo
                if len_kept <= 0 or i0 < pi_at or j0 < pj_at:
                    continue
                cost += aligner.distance(Span(a)[pi_at:i0], Span(b)[pj_at:j0])
                pi_at = i0 + len_kept
                pj_at = j0 + len_kept
            cost += aligner.distance(Span(a)[pi_at:n], Span(b)[pj_at:m])
            forced[pi] += cost - Int(best)
            for k in range(len(traced)):
                ref an = traced[k]
                var d = min(an.depth, DEPTHS - 1)
                var base = (pi * DEPTHS + d) * 4
                tally[base] += 1
                tally[base + 1] += 1 if off[k] > 0 else 0
                tally[base + 2] += letters[k]
                tally[base + 3] += off[k]
                by_class[pi * 3] += near[k]
                by_class[pi * 3 + 1] += copies[k]
                by_class[pi * 3 + 2] += other[k]
                text.write(
                    plan.text,
                    "\t",
                    pq[0],
                    "\t",
                    pq[1],
                    "\t",
                    n,
                    "\t",
                    m,
                    "\t",
                    an.depth,
                    "\t",
                    an.level,
                    "\t",
                    plan.levels[an.level].min_len,
                    "\t",
                    an.length,
                    "\t",
                    an.ref_copies,
                    "\t",
                    an.query_copies,
                    "\t",
                    an.win_ref,
                    "\t",
                    an.win_query,
                    "\t",
                    1 if an.in_repeat else 0,
                    "\t",
                    an.siblings,
                    "\t",
                    an.diag_start,
                    "\t",
                    an.diag_end,
                    "\t",
                    letters[k],
                    "\t",
                    off[k],
                    "\t",
                    near[k],
                    "\t",
                    copies[k],
                    "\t",
                    other[k],
                    "\t",
                    an.ref_pos,
                    "\t",
                    an.query_pos,
                    "\n",
                )
    with open(String(args[4]), "w") as f:
        f.write(text)
    print("pairs", done, "of", len(pairs), "period", period)
    for pi in range(n_plans):
        var anchors = 0
        var off_anchors = 0
        var letters = 0
        var off_letters = 0
        for d in range(DEPTHS):
            var base = (pi * DEPTHS + d) * 4
            anchors += tally[base]
            off_anchors += tally[base + 1]
            letters += tally[base + 2]
            off_letters += tally[base + 3]
        print(
            plans[pi].text,
            "| anchors",
            anchors,
            "off",
            off_anchors,
            "| letters",
            letters,
            "off",
            off_letters,
            "| off by: near",
            by_class[pi * 3],
            "copies",
            by_class[pi * 3 + 1],
            "other",
            by_class[pi * 3 + 2],
            "| forced excess",
            forced[pi],
            "of optima",
            optimum,
        )
        for d in range(DEPTHS):
            var base = (pi * DEPTHS + d) * 4
            if tally[base] > 0:
                print(
                    "    depth",
                    d,
                    "anchors",
                    tally[base],
                    "off",
                    tally[base + 1],
                    "letters",
                    tally[base + 2],
                    "off",
                    tally[base + 3],
                )
