#!/usr/bin/env python3
"""A port of minimap2's pinned alignment (2.31, -x asm5 as dipcall runs it: -z200000,10000 -r2k), forward strand,
one reference: what makes minimap2 (and so dipcall, and so HG002's truth) write a copy-number change in a tandem
array as several records where the exact global alignment writes one.

Ported line by line from sketch.c (mm_sketch), seed.c (mm_seed_mz_flt, mm_seed_collect_all, mm_seed_select,
mm_collect_matches), index.c (mm_idx_cal_max_occ), lchain.c (mg_lchain_rmq, comput_sc_simple, mg_chain_backtrack),
map.c (the long-join re-chain) and align.c (mm_fix_bad_ends, mm_filter_bad_seeds, mm_filter_bad_seeds_alt and
mm_align1's gap filling: a pin at each kept seed's middle once both sequences advanced MIN_KSW_LEN). The RMQ's ties
fall to the later anchor (minimap2's follow its tree's shape). The segments between pins are aligned globally by
karan at -A's costs (asm5's scores made costs), ties left; the ends too (minimap2 extends them).

  pins(ref, query) -> [(ref, query)] the points the alignment passes through
  chained_rows(pairs) -> {name: (ref row, query row)}, pairs {name: (ref, query)}: every segment in one karan run.
"""
import bisect
import math
import struct
import os
import re
import subprocess
import tempfile

K = 19
W = 19
MID_OCC_FRAC = 2e-4
MIN_MID_OCC = 50
MAX_MID_OCC = 500
MAX_MAX_OCC = 4095
OCC_DIST = 500
Q_OCC_FRAC = 0.01
MAX_GAP = 10000
RMQ_INNER_DIST = 1000
BW = 2000  # dipcall's -r2k
BW_LONG = 100000
MAX_CHAIN_SKIP = 25
MIN_CNT = 3
MIN_CHAIN_SCORE = 40
CHAIN_GAP_SCALE = 0.8
RMQ_RESCUE_SIZE = 1000
RMQ_RESCUE_RATIO = 0.1
MIN_KSW_LEN = 200
LONG_JOIN = 1 << 40
IGNORE = 1 << 41
TANDEM = 1 << 42
COSTS = "two-piece:40,78,7,162,3"
NT4 = {"A": 0, "C": 1, "G": 2, "T": 3}
U64 = (1 << 64) - 1


def hash64(key, mask):
    key = (~key + (key << 21)) & mask
    key = key ^ (key >> 24)
    key = ((key + (key << 3)) + (key << 8)) & mask
    key = key ^ (key >> 14)
    key = ((key + (key << 2)) + (key << 4)) & mask
    key = key ^ (key >> 28)
    key = (key + (key << 31)) & mask
    return key


def sketch(s, w=W, k=K):
    """mm_sketch: [(x = hash << 8 | span, y = pos << 1 | strand)]."""
    shift1, mask = 2 * (k - 1), (1 << (2 * k)) - 1
    kmer = [0, 0]
    buf = [(U64, U64)] * w
    mn = (U64, U64)
    out = []
    l = buf_pos = min_pos = 0
    for i, ch in enumerate(s):
        c = NT4.get(ch, 4)
        info = (U64, U64)
        if c < 4:
            span = l + 1 if l + 1 < k else k
            kmer[0] = ((kmer[0] << 2) | c) & mask
            kmer[1] = (kmer[1] >> 2) | ((3 ^ c) << shift1)
            if kmer[0] == kmer[1]:
                continue
            z = 0 if kmer[0] < kmer[1] else 1
            l += 1
            if l >= k and span < 256:
                info = (hash64(kmer[z], mask) << 8 | span, i << 1 | z)
        else:
            l = 0
        buf[buf_pos] = info
        if l == w + k - 1 and mn[0] != U64:
            for j in list(range(buf_pos + 1, w)) + list(range(0, buf_pos)):
                if mn[0] == buf[j][0] and buf[j][1] != mn[1]:
                    out.append(buf[j])
        if info[0] <= mn[0]:
            if l >= w + k and mn[0] != U64:
                out.append(mn)
            mn, min_pos = info, buf_pos
        elif buf_pos == min_pos:
            if l >= w + k - 1 and mn[0] != U64:
                out.append(mn)
            mn = (U64, U64)
            for j in list(range(buf_pos + 1, w)) + list(range(0, buf_pos + 1)):
                if mn[0] >= buf[j][0]:
                    mn, min_pos = buf[j], j
            if l >= w + k - 1 and mn[0] != U64:
                for j in list(range(buf_pos + 1, w)) + list(range(0, buf_pos + 1)):
                    if mn[0] == buf[j][0] and mn[1] != buf[j][1]:
                        out.append(buf[j])
        buf_pos += 1
        if buf_pos == w:
            buf_pos = 0
    if mn[0] != U64:
        out.append(mn)
    return out


def mid_occ(index):
    counts = sorted(len(v) for v in index.values())
    if not counts:
        return MIN_MID_OCC
    thres = counts[min(int((1.0 - MID_OCC_FRAC) * len(counts)), len(counts) - 1)] + 1
    return min(max(thres, MIN_MID_OCC), MAX_MID_OCC)


def seed_select(seeds, qlen, max_occ):
    """mm_seed_select: in each streak of high-occurrence seeds, up to its length / OCC_DIST of the rarest kept."""
    n = len(seeds)
    if n <= 1 or not any(s["n"] > max_occ for s in seeds):
        return
    last0 = -1
    for i in range(n + 1):
        if i == n or seeds[i]["n"] <= max_occ:
            if i - last0 > 1:
                ps = 0 if last0 < 0 else seeds[last0]["q_pos"] >> 1
                pe = qlen if i == n else seeds[i]["q_pos"] >> 1
                st, en = last0 + 1, i
                max_high = int((pe - ps) / OCC_DIST + 0.499)
                if max_high > 0:
                    chosen = sorted(range(st, en), key=lambda j: (seeds[j]["n"], j))[: min(max_high, 128)]
                    for j in chosen:
                        seeds[j]["flt"] = 1
                for j in range(st, en):
                    seeds[j]["flt"] ^= 1
                for j in range(st, en):
                    if seeds[j]["n"] > MAX_MAX_OCC:
                        seeds[j]["flt"] = 1
            last0 = i


def anchors(ref, query):
    """The forward-strand anchors (x = ref end, y = flags | span << 32 | query end), sorted as minimap2 sorts them."""
    index = {}
    for x, y in sketch(ref):
        index.setdefault(x >> 8, []).append(y)
    max_occ = mid_occ(index)
    mv = sketch(query)
    # mm_seed_mz_flt: a minimizer the query holds too often.
    if len(mv) > max_occ:
        counts = {}
        for x, _ in mv:
            counts[x] = counts.get(x, 0) + 1
        mv = [(x, y) for x, y in mv if not (counts[x] > max_occ and counts[x] > len(mv) * Q_OCC_FRAC)]
    seeds = []
    for i, (x, y) in enumerate(mv):
        hits = index.get(x >> 8)
        if not hits:
            continue
        tandem = (i > 0 and mv[i - 1][0] >> 8 == x >> 8) or (i + 1 < len(mv) and mv[i + 1][0] >> 8 == x >> 8)
        seeds.append({"q_pos": y, "q_span": x & 0xFF, "hits": hits, "n": len(hits), "tandem": tandem, "flt": 0})
    if MAX_MAX_OCC > max_occ:
        seed_select(seeds, len(query), max_occ)
    else:
        for s in seeds:
            s["flt"] = int(s["n"] > max_occ)
    out = []
    for s in seeds:
        if s["flt"]:
            continue
        for r in s["hits"]:
            if (r & 1) == (s["q_pos"] & 1):
                out.append((r >> 1, (TANDEM if s["tandem"] else 0) | s["q_span"] << 32 | s["q_pos"] >> 1))
    out.sort()
    return out


def qpos(a):
    return a[1] & 0xFFFFFFFF


def span(a):
    return (a[1] >> 32) & 0xFF


def f32(x):
    """`x` rounded to a C float."""
    return struct.unpack("f", struct.pack("f", x))[0]


def mg_log2(x):
    """minimap2's mg_log2: log2 from the float's exponent and a polynomial in its mantissa (not exact), in floats."""
    i = struct.unpack("I", struct.pack("f", float(x)))[0]
    log_2 = f32(((i >> 23) & 255) - 128)
    i = (i & ~(255 << 23) & 0xFFFFFFFF) + (127 << 23)
    z = struct.unpack("f", struct.pack("I", i))[0]
    return f32(log_2 + f32(f32(f32(f32(-0.34484843 * z) + 2.02466578) * z) - 0.67487759))


def comput_sc_simple(ai, aj, pen_gap, pen_skip):
    dq = qpos(ai) - qpos(aj)
    dr = ai[0] - aj[0]
    dd = abs(dr - dq)
    dg = min(dr, dq)
    q_span = span(aj)
    sc = min(q_span, dg)
    exact = dd == 0 and dg <= q_span
    if dd or dq > q_span:
        # In C floats: lin_pen = chn_pen_gap * dd + chn_pen_skip * dg; sc -= (int)(lin_pen + .5f * log_pen).
        lin = f32(f32(f32(pen_gap) * f32(dd)) + f32(f32(pen_skip) * f32(dg)))
        log_pen = mg_log2(dd + 1) if dd >= 1 else 0.0
        sc -= int(f32(lin + f32(0.5 * log_pen)))
    return sc, exact, dd


def lchain_rmq(a, bw, max_dist=MAX_GAP, max_dist_inner=RMQ_INNER_DIST):
    """mg_lchain_rmq and mg_chain_backtrack: the chains, best first, each its anchors in order."""
    n = len(a)
    if n == 0:
        return []
    pen_gap = CHAIN_GAP_SCALE * 0.01 * K
    max_dist = max(max_dist, bw)
    max_dist_inner = min(max(max_dist_inner, 0), max_dist)
    f, p, v, t = [0] * n, [-1] * n, [0] * n, [0] * n
    # The trees: anchors keyed (query position, index); the outer one a segment tree of each key's priority
    # (f + gap penalty x (ref + query) / 2: the RMQ's), the inner one a sorted list walked down from the anchor.
    keys = sorted((qpos(a[j]), j) for j in range(n))
    leaf = {key[1]: m for m, key in enumerate(keys)}
    size = 1
    while size < n:
        size <<= 1
    NONE = (-math.inf, -1, -1)
    tree = [NONE] * (2 * size)

    def tree_set(j, value):
        m = leaf[j] + size
        tree[m] = value
        m >>= 1
        while m:
            tree[m] = max(tree[2 * m], tree[2 * m + 1])
            m >>= 1

    def tree_max(lo, hi):
        best = NONE
        lo += size
        hi += size
        while lo < hi:
            if lo & 1:
                best = max(best, tree[lo])
                lo += 1
            if hi & 1:
                hi -= 1
                best = max(best, tree[hi])
            lo >>= 1
            hi >>= 1
        return best

    inner = []
    i0 = st = st_inner = 0
    for i in range(n):
        max_j, max_f = -1, span(a[i])
        if i0 < i and a[i0][0] != a[i][0]:
            for j in range(i0, i):
                tree_set(j, (f[j] + 0.5 * pen_gap * (a[j][0] + qpos(a[j])), qpos(a[j]), j))
                if max_dist_inner > 0:
                    bisect.insort(inner, (qpos(a[j]), j))
            i0 = i
        while st < i and a[i][0] > a[st][0] + max_dist:
            if st < i0:
                tree_set(st, NONE)
            st += 1
        if max_dist_inner > 0:
            while st_inner < i and a[i][0] > a[st_inner][0] + max_dist_inner:
                m = bisect.bisect_left(inner, (qpos(a[st_inner]), st_inner))
                if m < len(inner) and inner[m] == (qpos(a[st_inner]), st_inner):
                    del inner[m]
                st_inner += 1
        yi = qpos(a[i])
        best = tree_max(bisect.bisect_right(keys, (yi - max_dist, n)), bisect.bisect_left(keys, (yi, -1)))
        if best[2] >= 0:
            j = best[2]
            sc, exact, width = comput_sc_simple(a[i], a[j], pen_gap, 0.0)
            sc += f[j]
            if width <= bw and sc > max_f:
                max_f, max_j = sc, j
            if not exact and inner and yi > 0:
                n_skip = 0
                for m in range(bisect.bisect_right(inner, (yi - 1, n)) - 1, -1, -1):
                    j = inner[m][1]
                    if qpos(a[j]) < yi - max_dist_inner:
                        break
                    sc, _, width = comput_sc_simple(a[i], a[j], pen_gap, 0.0)
                    sc += f[j]
                    if width <= bw:
                        if sc > max_f:
                            max_f, max_j = sc, j
                            if n_skip > 0:
                                n_skip -= 1
                        elif t[j] == i:
                            n_skip += 1
                            if n_skip > MAX_CHAIN_SKIP:
                                break
                        if p[j] >= 0:
                            t[p[j]] = i
        f[i], p[i] = max_f, max_j
        v[i] = v[max_j] if max_j >= 0 and v[max_j] > max_f else max_f
    # mg_chain_backtrack, max_drop = bw.
    z = sorted((f[i], i) for i in range(n) if f[i] >= MIN_CHAIN_SCORE)
    t = [0] * n
    chains = []
    for k in range(len(z) - 1, -1, -1):
        zf, zi = z[k]
        if t[zi]:
            continue
        # mg_chain_bk_end
        i, end_i, max_i, max_s = zi, -1, zi, 0
        if t[i] == 0:
            while True:
                t[i] = 2
                end_i = i = p[i]
                s = zf if i < 0 else zf - f[i]
                if s > max_s:
                    max_s, max_i = s, i
                elif max_s - s > bw:
                    break
                if not (i >= 0 and t[i] == 0):
                    break
            i = zi
            while i >= 0 and i != end_i:
                t[i] = 0
                i = p[i]
        end_i = max_i
        members = []
        i = zi
        while i != end_i:
            members.append(i)
            t[i] = 1
            i = p[i]
        sc = zf if i < 0 else zf - f[i]
        if sc >= MIN_CHAIN_SCORE and len(members) >= MIN_CNT:
            chains.append((sc, [a[m] for m in reversed(members)]))
    chains.sort(key=lambda c: -c[0])
    return chains


def chain(ref, query):
    """The primary chain, as map.c gets it (re-chained with BW_LONG when there are several)."""
    a = anchors(ref, query)
    chains = lchain_rmq(a, BW)
    if len(chains) > 1:
        # map.c tests the chain compact_a puts first: chains sorted by their first anchor's reference position (not
        # by score; "u[] may not be sorted by score here").
        first = min(chains, key=lambda c: c[1][0][0])[1]
        st, en = qpos(first[0]), qpos(first[-1])
        if len(query) - (en - st) > RMQ_RESCUE_SIZE or en - st > len(query) * RMQ_RESCUE_RATIO:
            pooled = sorted(x for _, c in chains for x in c)
            chains = lchain_rmq(pooled, BW_LONG)
    return list(chains[0][1]) if chains else []


def fix_bad_ends(a, bw, min_match):
    """mm_fix_bad_ends: (first, count) of the chain's anchors kept."""
    cnt = len(a)
    if cnt < 3:
        return 0, cnt
    mlen = span(a[0])
    for i in range(1, cnt):
        tl, ql = a[i][0] - a[i - 1][0], qpos(a[i]) - qpos(a[i - 1])
        s = span(a[i])
        mlen += s if tl > s and ql > s else min(tl, ql)
    as_ = 0
    m = l = span(a[0])
    for i in range(1, cnt - 1):
        if a[i][1] & LONG_JOIN:
            break
        lr, lq = a[i][0] - a[i - 1][0], qpos(a[i]) - qpos(a[i - 1])
        mn, mx = min(lr, lq), max(lr, lq)
        if mx - mn > l >> 1:
            as_ = i
        l += mn
        m += min(mn, span(a[i]))
        if l >= bw << 1 or (m >= min_match and m >= bw) or m >= mlen >> 1:
            break
    count = cnt - as_
    m = l = span(a[cnt - 1])
    for i in range(cnt - 2, as_, -1):
        if a[i + 1][1] & LONG_JOIN:
            break
        lr, lq = a[i + 1][0] - a[i][0], qpos(a[i + 1]) - qpos(a[i])
        mn, mx = min(lr, lq), max(lr, lq)
        if mx - mn > l >> 1:
            count = i + 1 - as_
        l += mn
        m += min(mn, span(a[i + 1]))
        if l >= bw << 1 or (m >= min_match and m >= bw) or m >= mlen >> 1:
            break
    return as_, count


def gap_of(a, i):
    return (qpos(a[i]) - qpos(a[i - 1])) - (a[i][0] - a[i - 1][0])


def long_gaps(a, min_gap):
    out = [i for i in range(1, len(a)) if abs(gap_of(a, i)) > min_gap]
    return out if len(out) > 1 else []


def flag(a, i, bit):
    a[i] = (a[i][0], a[i][1] | bit)


def filter_bad_seeds(a, min_gap=10, diff_thres=40, max_ext_len=MAX_GAP >> 1, max_ext_cnt=10):
    ks = long_gaps(a, min_gap)
    if not ks:
        return
    n = len(ks)
    best, max_st, max_en = 0, -1, -1
    k = 0
    while True:
        if k == n or k >= max_en:
            if max_en > 0:
                for i in range(ks[max_st], ks[max_en]):
                    flag(a, i, IGNORE)
            best, max_st, max_en = 0, -1, -1
            if k == n:
                break
        i = ks[k]
        g = gap_of(a, i)
        n_ins, n_del = max(g, 0), max(-g, 0)
        qs, rs = qpos(a[i - 1]), a[i - 1][0]
        max_diff, max_diff_l = 0, -1
        for l in range(k + 1, min(n, k + max_ext_cnt + 1)):
            j = ks[l]
            if qpos(a[j]) - qs > max_ext_len or a[j][0] - rs > max_ext_len:
                break
            g = gap_of(a, j)
            if g > 0:
                n_ins += g
            else:
                n_del += -g
            diff = n_ins + n_del - abs(n_ins - n_del)
            if max_diff < diff:
                max_diff, max_diff_l = diff, l
        if max_diff > diff_thres and max_diff > best:
            best, max_st, max_en = max_diff, k, max_diff_l
        k += 1


def filter_bad_seeds_alt(a, min_gap=30, max_ext=MAX_GAP >> 1):
    ks = long_gaps(a, min_gap)
    if not ks:
        return
    n = len(ks)
    k = 0
    while k < n:
        i = ks[k]
        gap1 = abs(gap_of(a, i))
        re1, qe1 = a[i][0], qpos(a[i])
        l = k + 1
        while l < n:
            j = ks[l]
            if qpos(a[j]) - qe1 > max_ext or a[j][0] - re1 > max_ext:
                break
            gap2 = abs(gap_of(a, j))
            pre = span(a[j - 1])
            rs2, qs2 = a[j - 1][0] + pre, qpos(a[j - 1]) + pre
            if min(rs2 - re1, qs2 - qe1) > gap1 + gap2:
                break
            re1, qe1 = a[j][0], qpos(a[j])
            gap1 = gap2
            l += 1
        if l > k + 1:
            end = ks[l - 1]
            for j in range(ks[k], end):
                flag(a, j, IGNORE)
            flag(a, end, LONG_JOIN)
        k = l


def pins(ref, query):
    a = chain(ref, query)
    if not a:
        return []
    as_, cnt = fix_bad_ends(a, BW, MIN_CHAIN_SCORE * 2)
    a = a[as_ : as_ + cnt]
    filter_bad_seeds(a)
    filter_bad_seeds_alt(a)
    half = K >> 1
    rs, qs = a[0][0] - half, qpos(a[0]) - half
    out = [(rs, qs)]
    for i in range(1, len(a)):
        if (a[i][1] & (IGNORE | TANDEM)) and i != len(a) - 1:
            continue
        re_, qe = a[i][0] - half, qpos(a[i]) - half
        if i == len(a) - 1 or (a[i][1] & LONG_JOIN) or (qe - qs >= MIN_KSW_LEN and re_ - rs >= MIN_KSW_LEN):
            out.append((re_, qe))
            rs, qs = re_, qe
    return [(r, q) for r, q in out if 0 < r < len(ref) and 0 < q < len(query)]


def pins_in_context(ref, query, left, right):
    """The pins of `left + ref + right` against `left + query + right` (the flanks only context for the seeds and the
    chain) that fall inside `ref` and `query`, in their own coordinates."""
    n = len(left)
    return [(r - n, q - n) for r, q in pins(left + ref + right, left + query + right)
            if 0 < r - n < len(ref) and 0 < q - n < len(query)]


def chained_rows(pairs, context=None):
    """`context` {name: (left, right)}: the flanks the pins are found with."""
    segments = {}
    for name, (ref, query) in pairs.items():
        found = pins_in_context(ref, query, *context[name]) if context else pins(ref, query)
        points = [(0, 0)] + found + [(len(ref), len(query))]
        segments[name] = [(p0, p1) for p0, p1 in zip(points, points[1:]) if p1 != p0]
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "pairs.tsv")
        with open(path, "w") as f:
            for name, segs in segments.items():
                ref, query = pairs[name]
                for n, ((r0, q0), (r1, q1)) in enumerate(segs):
                    if r1 > r0 and q1 > q0:
                        f.write("%s:%d\t%s\t%s\n" % (name, n, ref[r0:r1], query[q0:q1]))
        out = subprocess.run([os.environ.get("KARAN", "karan"), "--pairs", path, "--costs", COSTS, "--ties", "left",
                              "--max-memory", str(8 << 30)], capture_output=True, text=True, check=True).stdout
    cigars = {}
    for line in out.splitlines()[1:]:
        w = line.split("\t")
        cigars[w[0]] = w[9]
    rows = {}
    for name, segs in segments.items():
        ref, query = pairs[name]
        rr, qq = [], []
        for n, ((r0, q0), (r1, q1)) in enumerate(segs):
            cigar = cigars.get("%s:%d" % (name, n))
            if cigar is None:
                rr.append(ref[r0:r1] + "-" * (q1 - q0))
                qq.append("-" * (r1 - r0) + query[q0:q1])
                continue
            i, j = r0, q0
            for count, op in re.findall(r"(\d+)([=XMID])", cigar):
                count = int(count)
                if op in "=XM":
                    rr.append(ref[i : i + count])
                    qq.append(query[j : j + count])
                    i += count
                    j += count
                elif op == "I":
                    rr.append("-" * count)
                    qq.append(query[j : j + count])
                    j += count
                else:
                    rr.append(ref[i : i + count])
                    qq.append("-" * count)
                    i += count
        rows[name] = ("".join(rr), "".join(qq))
    return rows
