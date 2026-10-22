from std.testing import assert_equal, assert_true, TestSuite

from karan import distance

from plegma.anchors import Anchor, find_anchors, mums, smems
from plegma.backbone import Level, MUMS, ALL, RepeatMask, backbone, windows
from plegma.chain import ChainCosts, chain, linked_cost, without_overlaps


def _random(n: Int, seed: Int) -> List[UInt8]:
    var out = List[UInt8](capacity=n)
    var state = UInt64(seed) * 0x9E3779B97F4A7C15 + 1
    for _ in range(n):
        state ^= state >> 12
        state ^= state << 25
        state ^= state >> 27
        out.append(UInt8(ord("ACGT"[byte=Int((state * 0x2545F4914F6CDD1D) >> 62)])))
    return out^


def _vntr(copies: Int, unit: StringSlice, seed: Int) -> List[UInt8]:
    """`copies` of `unit` between random flanks, one copy in four changed at one letter."""
    var out = _random(60, seed)
    var state = UInt64(seed) * 0x2545F4914F6CDD1D + 7
    for _ in range(copies):
        state ^= state >> 12
        state ^= state << 25
        state ^= state >> 27
        var at = len(out)
        out.extend(unit.as_bytes())
        if state % 4 == 0:
            out[at + Int((state >> 20) % UInt64(unit.byte_length()))] = UInt8(ord("ACGT"[byte=Int(state >> 62)]))
    out.extend(_random(60, seed + 1))
    return out^


def _mutated(src: List[UInt8], every: Int, seed: Int) -> List[UInt8]:
    """A copy of `src` with a substitution, deletion or insertion about once per `every` letters."""
    var out = List[UInt8]()
    var state = UInt64(seed) * 0x9E3779B97F4A7C15 + 3
    for c in src:
        state ^= state >> 12
        state ^= state << 25
        state ^= state >> 27
        var r = state * 0x2545F4914F6CDD1D
        if Int(r % UInt64(every)) == 0:
            var kind = (r >> 40) % 3
            if kind == 0:
                out.append(UInt8(ord("ACGT"[byte=Int(r >> 62)])))
            elif kind == 2:
                out.append(c)
                out.append(UInt8(ord("ACGT"[byte=Int(r >> 62)])))
            continue
        out.append(c)
    return out^


def _anchored_distance(r: List[UInt8], q: List[UInt8], anchors: List[Anchor]) -> Int:
    """Edit distance where a match counts only on an anchor's diagonal inside it (ChainX's anchored distance)."""
    var allowed = List[Bool](length=(len(r) + 1) * (len(q) + 1), fill=False)
    for a in anchors:
        for k in range(a.length):
            allowed[(a.ref_pos + k + 1) * (len(q) + 1) + a.query_pos + k + 1] = True
    var prev = List[Int](capacity=len(q) + 1)
    for j in range(len(q) + 1):
        prev.append(j)
    for i in range(1, len(r) + 1):
        var cur = List[Int](length=len(q) + 1, fill=i)
        for j in range(1, len(q) + 1):
            var diagonal = prev[j - 1] + (0 if r[i - 1] == q[j - 1] and allowed[i * (len(q) + 1) + j] else 1)
            cur[j] = min(diagonal, prev[j] + 1, cur[j - 1] + 1)
        prev = cur^
    return prev[len(q)]


def _text(bytes: List[UInt8]) -> String:
    var out = String()
    for c in bytes:
        out += chr(Int(c))
    return out^


def test_chainx_cost_is_the_anchored_distance() raises:
    for case_id in range(8):
        var r = _vntr(6 + case_id, "GGTGGTTAGGGTCACGGT", case_id) if case_id % 2 == 0 else _random(300, case_id)
        var q = _mutated(r, 25, case_id + 50)
        for min_len in [8, 12]:
            var anchors = find_anchors(Span(r), Span(q), min_len, 1000)
            var found = chain(anchors, len(r), len(q), ChainCosts.chainx())
            assert_equal(Int(found.cost), _anchored_distance(r, q, anchors))


def test_filled_backbone_is_never_below_the_optimum() raises:
    for case_id in range(6):
        var r = _vntr(10 + 3 * case_id, "TGATGGTGGTGGTTGGTAGTATTACTGTTGGTGG", case_id)
        var q = _mutated(_vntr(7 + 2 * case_id, "TGATGGTGGTGGTTGGTAGTATTACTGTTGGTGG", case_id), 60, case_id + 9)
        var best = distance(_text(r), _text(q))
        var levels: List[Level] = [Level(MUMS, 31), Level(MUMS, 15)]
        for repeat in [False, True]:
            var picked = backbone(Span(r), Span(q), levels, ChainCosts.chainx(), 0, repeat)
            var ref_at = 0
            var query_at = 0
            for a in picked:
                assert_true(a.ref_pos >= ref_at and a.query_pos >= query_at)
                for k in range(a.length):
                    assert_equal(r[a.ref_pos + k], q[a.query_pos + k])
                ref_at = a.ref_end()
                query_at = a.query_end()
            var filled = 0
            for w in windows(picked, len(r), len(q)):
                var rs = List[UInt8](r[w.ref_start : w.ref_end])
                var qs = List[UInt8](q[w.query_start : w.query_end])
                filled += distance(_text(rs), _text(qs))
            assert_true(filled >= best)
            var cost = linked_cost(picked, len(r), len(q), ChainCosts.chainx())
            assert_true(Float64(filled) <= cost)


def test_smems_keep_only_query_maximal_matches() raises:
    var r = _vntr(12, "ACGTTGCAAGT", 3)
    var q = _vntr(9, "ACGTTGCAAGT", 4)
    var every = find_anchors(Span(r), Span(q), 10, 1000)
    var kept = smems(every)
    for a in kept:
        for b in every:
            var contains = b.query_pos <= a.query_pos and b.query_end() >= a.query_end()
            assert_true(not contains or b.length == a.length)
    for b in every:
        var covered = False
        for a in kept:
            covered = covered or (a.query_pos <= b.query_pos and a.query_end() >= b.query_end())
        assert_true(covered)


def test_without_overlaps_clears_each_predecessor() raises:
    var anchors: List[Anchor] = [Anchor(0, 0, 10), Anchor(6, 5, 10), Anchor(12, 12, 3), Anchor(20, 20, 4)]
    var cut = without_overlaps(anchors)
    assert_equal(len(cut), 3)
    assert_equal(cut[1].ref_pos, 11)
    assert_equal(cut[1].query_pos, 10)
    assert_equal(cut[1].length, 5)
    assert_equal(cut[2].ref_pos, 20)


def test_repeat_mask_marks_tandem_copies() raises:
    var flank = _random(200, 77)
    var text = List[UInt8](flank)
    for _ in range(20):
        text.extend("ACGTTGCAGGTACCATG".as_bytes())
    text.extend(_random(200, 78))
    var mask = RepeatMask(Span(text), 12, 500)
    assert_true(mask.covers(220, 260))
    assert_true(not mask.covers(10, 150))
    assert_true(not mask.covers(len(text) - 150, len(text) - 10))


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
