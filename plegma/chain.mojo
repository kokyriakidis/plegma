"""Colinear chaining: the anchors an alignment of two sequences is pinned to.

A chain runs from the start of both sequences to their ends through anchors in order on both. It scores what its
anchors earn less what each link costs (`ChainCosts`). With anchors earning nothing, a link costing the longer side
of its gap and overlaps allowed, this is ChainX's chaining (Jain, Gibney and Thankachan 2022), whose optimum is the
edit distance with every match outside the anchors counted as an edit.
"""

from .anchors import Anchor


struct ChainCosts(ImplicitlyCopyable):
    """What anchors earn and links cost.

    An anchor earns `gain` per letter it adds past its predecessor, times (ref copies x query copies)^(-1/2) with
    `rarity`. A link with gaps of `g_r` partner and `g_q` query letters costs `substitution` per letter of the
    shorter, and for its indel `d = g_r - g_q` the least of `open1 + extend1 |d|` and `open2 + extend2 |d|` (0 when
    `d` is 0), plus `off_period` per letter `|d|` lies from a multiple of `period` (none with `period` 0).
    """

    var gain: Float64
    var rarity: Bool
    var substitution: Float64
    var open1: Float64
    var extend1: Float64
    var open2: Float64
    var extend2: Float64
    var period: Int
    var off_period: Float64
    var overlaps: Bool
    """Whether an anchor may overlap its predecessor (consecutive repeat copies do)."""

    def __init__(
        out self,
        *,
        gain: Float64 = 0.0,
        rarity: Bool = False,
        substitution: Float64 = 1.0,
        open1: Float64 = 0.0,
        extend1: Float64 = 1.0,
        open2: Float64 = 0.0,
        extend2: Float64 = 1.0,
        period: Int = 0,
        off_period: Float64 = 0.0,
        overlaps: Bool = True,
    ):
        self.gain = gain
        self.rarity = rarity
        self.substitution = substitution
        self.open1 = open1
        self.extend1 = extend1
        self.open2 = open2
        self.extend2 = extend2
        self.period = period
        self.off_period = off_period
        self.overlaps = overlaps

    @staticmethod
    def heaviest() -> Self:
        """MUMmer's chain: the most anchored letters, links free, no overlaps."""
        return Self(gain=1.0, substitution=0.0, extend1=0.0, extend2=0.0, overlaps=False)

    @staticmethod
    def chainx() -> Self:
        """ChainX's: links cost the longer side of their gap, anchors earn nothing."""
        return Self()

    def indel(self, d: Int) -> Float64:
        """What a link's indel of `d` letters costs."""
        if d == 0:
            return 0.0
        var size = Float64(abs(d))
        var cost = min(self.open1 + self.extend1 * size, self.open2 + self.extend2 * size)
        if self.period > 0:
            var off = abs(d) % self.period
            cost += self.off_period * Float64(min(off, self.period - off))
        return cost

    def link(self, before: Anchor, after: Anchor) -> Float64:
        """What linking `before` to `after` costs."""
        var query_gap = after.query_pos - before.query_end()
        var ref_gap = after.ref_pos - before.ref_end()
        return self.substitution * Float64(max(0, min(query_gap, ref_gap))) + self.indel(ref_gap - query_gap)

    def earns(self, before: Anchor, after: Anchor) -> Float64:
        """What `after` earns following `before`."""
        if self.gain == 0.0:
            return 0.0
        var overlap = max(0, before.query_end() - after.query_pos, before.ref_end() - after.ref_pos)
        var weight = self.gain
        if self.rarity:
            weight /= Float64(after.ref_copies * after.query_copies) ** 0.5
        return weight * Float64(max(0, after.length - overlap))

    def precedes(self, before: Anchor, after: Anchor) -> Bool:
        """Whether `after` may follow `before` in a chain."""
        if self.overlaps:
            return (
                before.query_pos <= after.query_pos
                and before.query_end() <= after.query_end()
                and before.ref_pos <= after.ref_pos
                and before.ref_end() <= after.ref_end()
                and before != after
            )
        return before.query_end() <= after.query_pos and before.ref_end() <= after.ref_pos


struct Chain(Movable):
    """An optimal chain: its anchors in order, its score and what its links cost."""

    var anchors: List[Anchor]
    var score: Float64
    var cost: Float64

    def __init__(out self, var anchors: List[Anchor], score: Float64, cost: Float64):
        self.anchors = anchors^
        self.score = score
        self.cost = cost


def chain(anchors: List[Anchor], ref_len: Int, query_len: Int, costs: ChainCosts) -> Chain:
    """The best chain from the start of both sequences to their ends through `anchors` (in query order, as
    `find_anchors` gives them). Of equal chains, the one taking the earliest predecessor at each anchor. Quadratic
    in the anchors at worst; under ChainX-like costs each anchor looks back only as far as its best link can
    reach."""
    var n = len(anchors)
    var start = Anchor(0, 0, 0)
    var end = Anchor(ref_len, query_len, 0)
    var score = List[Float64](length=n, fill=0.0)
    var paid = List[Float64](length=n, fill=0.0)
    var back = List[Int](length=n, fill=-1)
    # With anchors earning nothing and indels costing at least `substitution` a letter, a link costs at least
    # `substitution` per letter of its query gap, and every score is at most 0: predecessors further back than
    # the best so far can pay for are never taken, and the scan back stops there.
    var bounded = (
        costs.gain == 0.0
        and costs.substitution > 0.0
        and costs.period == 0
        and costs.open1 >= 0.0
        and costs.open2 >= 0.0
        and min(costs.extend1, costs.extend2) >= costs.substitution
    )
    var longest = 0
    for a in anchors:
        longest = max(longest, a.length)
    for i in range(n):
        var a = anchors[i]
        var best = costs.earns(start, a) - costs.link(start, a)
        var best_paid = costs.link(start, a)
        var from_index = -1
        # Back from the nearest, the earliest of equal predecessors kept (and the start before any of them).
        for j in reversed(range(i)):
            var b = anchors[j]
            if bounded and costs.substitution * Float64(a.query_pos - b.query_pos - longest) > -best:
                break
            if not costs.precedes(b, a):
                continue
            var link = costs.link(b, a)
            var value = score[j] + costs.earns(b, a) - link
            if value > best or (value == best and from_index >= 0):
                best = value
                best_paid = paid[j] + link
                from_index = j
        score[i] = best
        paid[i] = best_paid
        back[i] = from_index
    var best = -costs.link(start, end)
    var best_paid = costs.link(start, end)
    var last = -1
    for i in range(n):
        var a = anchors[i]
        if a.query_end() > query_len or a.ref_end() > ref_len:
            continue
        var link = costs.link(a, end)
        if score[i] - link > best:
            best = score[i] - link
            best_paid = paid[i] + link
            last = i
    var picked = List[Anchor]()
    var at = last
    while at >= 0:
        picked.append(anchors[at])
        at = back[at]
    picked.reverse()
    return Chain(picked^, best, best_paid)


def linked_cost(anchors: List[Anchor], ref_len: Int, query_len: Int, costs: ChainCosts) -> Float64:
    """What linking `anchors` in their order costs, from the start of both sequences to their ends."""
    var before = Anchor(0, 0, 0)
    var total = 0.0
    for a in anchors:
        total += costs.link(before, a)
        before = a
    return total + costs.link(before, Anchor(ref_len, query_len, 0))


def without_overlaps(anchors: List[Anchor]) -> List[Anchor]:
    """A chain's anchors, each cut at its start until it clears its predecessor on both sequences; an anchor cut
    away whole is dropped."""
    var out = List[Anchor]()
    for a in anchors:
        var cut = 0
        if len(out) > 0:
            var before = out[len(out) - 1]
            cut = max(0, before.query_end() - a.query_pos, before.ref_end() - a.ref_pos)
        if cut < a.length:
            out.append(Anchor(a.ref_pos + cut, a.query_pos + cut, a.length - cut, a.ref_copies, a.query_copies))
    return out^
