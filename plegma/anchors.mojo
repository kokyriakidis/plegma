"""Exact-match anchors between a query and a partner sequence, from sporos.

Every MEM of a pair comes with how often its string occurs in each sequence, so MUMs (once in both), rare MEMs
(a few times) and SMEMs (no other match's query interval contains it) are views of one list.
"""

from sporos.mems import MemIndex, MemParams, find_mems, occurrence_counts
from sporos.text import SeqText


struct Anchor(Equatable, TrivialRegisterPassable):
    """An exact match: `length` letters from `ref_pos` in the partner and `query_pos` in the query."""

    var ref_pos: Int
    var query_pos: Int
    var length: Int
    var ref_copies: Int
    """How often the anchor's string occurs in the partner, at most one past the cap it was counted to."""
    var query_copies: Int
    """How often it occurs in the query, likewise."""

    def __init__(out self, ref_pos: Int, query_pos: Int, length: Int, ref_copies: Int = 1, query_copies: Int = 1):
        self.ref_pos = ref_pos
        self.query_pos = query_pos
        self.length = length
        self.ref_copies = ref_copies
        self.query_copies = query_copies

    def __eq__(self, other: Self) -> Bool:
        return self.ref_pos == other.ref_pos and self.query_pos == other.query_pos and self.length == other.length

    @inline(.always)
    def ref_end(self) -> Int:
        return self.ref_pos + self.length

    @inline(.always)
    def query_end(self) -> Int:
        return self.query_pos + self.length

    @inline(.always)
    def diagonal(self) -> Int:
        """The partner's offset less the query's, as karan counts diagonals."""
        return self.ref_pos - self.query_pos

    @inline(.always)
    def is_mum(self) -> Bool:
        return self.ref_copies == 1 and self.query_copies == 1

    def shifted(self, ref_offset: Int, query_offset: Int) -> Self:
        """The same anchor, its positions moved by the offsets."""
        return Self(
            self.ref_pos + ref_offset, self.query_pos + query_offset, self.length, self.ref_copies, self.query_copies
        )


def find_anchors(partner: ImmSpan[UInt8, _], query: ImmSpan[UInt8, _], min_len: Int, cap: Int) raises -> List[Anchor]:
    """Every MEM of at least `min_len` letters between `partner` and `query` (forward strands), in query order,
    each with its copies in both counted to `cap`."""
    var anchors = List[Anchor]()
    if len(partner) < min_len or len(query) < min_len:
        return anchors^
    var refs = List[type_of(partner)]()
    refs.append(partner)
    var queries = List[type_of(query)]()
    queries.append(query)
    var rtext = SeqText(refs)
    var qtext = SeqText(queries)
    var index = MemIndex(rtext, MemParams.choose(min_len, len(rtext), len(qtext)))
    var mems = find_mems(index, rtext, qtext)
    var counts = occurrence_counts(mems, qtext, cap, cap)
    var ref_start = rtext.starts[0]
    var query_start = qtext.starts[0]
    anchors.reserve(len(mems))
    for i in range(len(mems)):
        var m = mems[i]
        anchors.append(
            Anchor(
                Int(m.ref_pos) - ref_start,
                Int(m.query_pos) - query_start,
                Int(m.length),
                Int(counts[0][i]),
                Int(counts[1][i]),
            )
        )
    return anchors^


def mums(anchors: List[Anchor]) -> List[Anchor]:
    """The anchors whose string occurs once in each sequence."""
    var out = List[Anchor]()
    for a in anchors:
        if a.is_mum():
            out.append(a)
    return out^


def rare(anchors: List[Anchor], max_copies: Int) -> List[Anchor]:
    """The anchors whose string occurs at most `max_copies` times in each sequence."""
    var out = List[Anchor]()
    for a in anchors:
        if a.ref_copies <= max_copies and a.query_copies <= max_copies:
            out.append(a)
    return out^


def smems(anchors: List[Anchor]) -> List[Anchor]:
    """The anchors whose query interval no other anchor's strictly contains: a string's every partner
    occurrence is kept. `anchors` must be in query order, as `find_anchors` gives them."""
    var out = List[Anchor]()
    var reach = -1  # the furthest query end of the anchors starting before the current start
    var i = 0
    while i < len(anchors):
        var start = anchors[i].query_pos
        var j = i
        var longest = 0
        while j < len(anchors) and anchors[j].query_pos == start:
            longest = max(longest, anchors[j].length)
            j += 1
        if start + longest > reach:
            for k in range(i, j):
                if anchors[k].length == longest:
                    out.append(anchors[k])
        reach = max(reach, start + longest)
        i = j
    return out^
