"""The order sequences join the graph: abPOA's guide order, from the similarity of their minimizers.

Each sequence's minimizers (sporos, canonical k-mers) form a multiset; two sequences' similarity is the multiset
Jaccard index of theirs. The most similar pair goes first, then each time the sequence most similar in total to
those already placed.
"""

from sporos import minimizer_positions

from .alphabet import encode_base


def _minimizer_keys(seq: Span[UInt8, _], k: Int, w: Int) raises -> List[UInt64]:
    """The 2-bit forward k-mers at the sequence's minimizer positions, sorted (a multiset)."""
    var keys = List[UInt64]()
    if len(seq) < k + w:
        return keys^
    var positions = minimizer_positions(seq, k=k, w=w)
    var mask = (UInt64(1) << UInt64(2 * k)) - 1
    for p in positions:
        var key = UInt64(0)
        var ok = True
        for i in range(k):
            var code = encode_base(seq[Int(p) + i])
            if code > 3:
                ok = False
                break
            key = ((key << 2) | UInt64(code)) & mask
        if ok:
            keys.append(key)
    sort(keys)
    return keys^


def _shared(a: List[UInt64], b: List[UInt64]) -> Int:
    """The size of the multiset intersection of two sorted lists."""
    var i = 0
    var j = 0
    var count = 0
    while i < len(a) and j < len(b):
        if a[i] == b[j]:
            count += 1
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
    return count


def _before(a: List[UInt8], b: List[UInt8]) -> Bool:
    """A content order on sequences, for ties: the longer first, then the first in letter order (identical ones are
    interchangeable)."""
    if len(a) != len(b):
        return len(a) > len(b)
    for i in range(len(a)):
        if a[i] != b[i]:
            return a[i] < b[i]
    return False


def guide_order(seqs: List[List[UInt8]], k: Int = 19, w: Int = 10) raises -> List[Int]:
    """abPOA's progressive order of `seqs`: the most similar pair by minimizer multiset Jaccard first, then repeatedly
    the sequence with the highest summed similarity to those placed. Of the first pair, and in ties, the sequences
    decide, not their places in `seqs` (`_before`: the longer, then the first in letter order), so the order is the
    same however the reads come: abPOA's (the earlier first, ties to the earlier) made a read set's consensus
    depend on its input order. (The pair's more central read first instead cost a simulated set a letter.)"""
    var n = len(seqs)
    var order = List[Int]()
    if n <= 2:
        for i in range(n):
            order.append(i)
        if n == 2 and _before(seqs[1], seqs[0]):
            order.reverse()
        return order^
    var keys = List[List[UInt64]]()
    for s in seqs:
        keys.append(_minimizer_keys(Span(s), k, w))
    var similar = List[Float64](length=n * n, fill=0.0)
    for i in range(n):
        for j in range(i + 1, n):
            var both = _shared(keys[i], keys[j])
            var either = len(keys[i]) + len(keys[j]) - both
            var value = Float64(both) / Float64(either) if either > 0 else 0.0
            similar[i * n + j] = value
            similar[j * n + i] = value

    def ahead(i: Int, j: Int) capturing -> Bool:
        """Whether `i` goes before `j` of a pair."""
        return _before(seqs[i], seqs[j])

    var best = -1.0
    var first = 0
    var second = 1
    for i in range(n):
        for j in range(i + 1, n):
            var value = similar[i * n + j]
            var a = i if ahead(i, j) else j
            var b = j if a == i else i
            # As similar: the pair whose first goes before the other pair's first, then whose second does.
            if value > best or (
                value == best and (ahead(a, first) and a != first or (a == first and b != second and ahead(b, second)))
            ):
                best = value
                first = a
                second = b
    var placed = List[Bool](length=n, fill=False)
    var total = List[Float64](length=n, fill=0.0)
    for x in [first, second]:
        order.append(x)
        placed[x] = True
        for i in range(n):
            total[i] += similar[x * n + i]
    while len(order) < n:
        var pick = -1
        for i in range(n):
            if placed[i]:
                continue
            if pick < 0 or total[i] > total[pick] or (total[i] == total[pick] and _before(seqs[i], seqs[pick])):
                pick = i
        order.append(pick)
        placed[pick] = True
        for i in range(n):
            total[i] += similar[pick * n + i]
    return order^
