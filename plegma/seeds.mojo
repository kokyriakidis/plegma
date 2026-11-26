"""minimap2 2.31 minimizer sketching, reference indexing, seed collection, and LSD radix sort.

Extracted from `pins.mojo` so sketching, reference index reuse (across many queries in `on_reference`), seed
filtering (`mm_seed_mz_flt`, `mm_seed_select`), and radix sorting can be tested and reused independently.
"""

from .alphabet import base_code

comptime K = 19
comptime W = 19
comptime MID_OCC_FRAC = 2e-4
comptime MIN_MID_OCC = 50
comptime MAX_MID_OCC = 500
comptime MAX_MAX_OCC = 4095
comptime OCC_DIST = 500
comptime Q_OCC_FRAC = 0.01
comptime LONG_JOIN = 1
comptime IGNORE = 2
comptime TANDEM = 4
comptime NONE = UInt64.MAX
comptime POSITION_BITS = 25
"""Sequences past 2^25 letters are not pinned (the sort key's room)."""


@fieldwise_init
struct Seed(ImplicitlyCopyable, TrivialRegisterPassable):
    """A chained seed: its reference and query end letters, its span and minimap2's flags (`TANDEM`, `IGNORE`,
    `LONG_JOIN`)."""

    var x: Int
    var q: Int
    var span: Int
    var flags: Int


struct Sketch(Movable):
    """mm_sketch's minimizers: `x` hash << 8 | span, `y` end << 1 | strand."""

    var x: List[UInt64]
    var y: List[UInt64]

    def __init__(out self):
        self.x = List[UInt64]()
        self.y = List[UInt64]()

    def push(mut self, x: UInt64, y: UInt64):
        self.x.append(x)
        self.y.append(y)


def _hash64(key: UInt64, mask: UInt64) -> UInt64:
    var k = (~key + (key << 21)) & mask
    k = k ^ (k >> 24)
    k = ((k + (k << 3)) + (k << 8)) & mask
    k = k ^ (k >> 14)
    k = ((k + (k << 2)) + (k << 4)) & mask
    k = k ^ (k >> 28)
    k = (k + (k << 31)) & mask
    return k


def sketch(s: List[UInt8]) -> Sketch:
    """mm_sketch: the (W, K) minimizers, every k-mer tied for a window's least kept, symmetric k-mers skipped."""
    var out = Sketch()
    var shift1 = UInt64(2 * (K - 1))
    var mask = (UInt64(1) << UInt64(2 * K)) - 1
    var fwd = UInt64(0)
    var rev = UInt64(0)
    var buf_x = Array[UInt64, W](fill=NONE)
    var buf_y = Array[UInt64, W](fill=NONE)
    var min_x = NONE
    var min_y = NONE
    var l = 0
    var buf_pos = 0
    var min_pos = 0
    for i in range(len(s)):
        var c = base_code(s[i])
        var info_x = NONE
        var info_y = NONE
        if c < 4:
            var span = l + 1 if l + 1 < K else K
            fwd = ((fwd << 2) | UInt64(c)) & mask
            rev = (rev >> 2) | (UInt64(3 ^ c) << shift1)
            if fwd == rev:
                continue
            var z = 0 if fwd < rev else 1
            l += 1
            if l >= K and span < 256:
                info_x = (_hash64(fwd if z == 0 else rev, mask) << 8) | UInt64(span)
                info_y = (UInt64(i) << 1) | UInt64(z)
        else:
            l = 0
        buf_x[buf_pos] = info_x
        buf_y[buf_pos] = info_y
        if l == W + K - 1 and min_x != NONE:
            for j in range(buf_pos + 1, W):
                if min_x == buf_x[j] and buf_y[j] != min_y:
                    out.push(buf_x[j], buf_y[j])
            for j in range(buf_pos):
                if min_x == buf_x[j] and buf_y[j] != min_y:
                    out.push(buf_x[j], buf_y[j])
        if info_x <= min_x:
            if l >= W + K and min_x != NONE:
                out.push(min_x, min_y)
            min_x = info_x
            min_y = info_y
            min_pos = buf_pos
        elif buf_pos == min_pos:
            if l >= W + K - 1 and min_x != NONE:
                out.push(min_x, min_y)
            min_x = NONE
            for j in range(buf_pos + 1, W):
                if min_x >= buf_x[j]:
                    min_x = buf_x[j]
                    min_y = buf_y[j]
                    min_pos = j
            for j in range(buf_pos + 1):
                if min_x >= buf_x[j]:
                    min_x = buf_x[j]
                    min_y = buf_y[j]
                    min_pos = j
            if l >= W + K - 1 and min_x != NONE:
                for j in range(buf_pos + 1, W):
                    if min_x == buf_x[j] and min_y != buf_y[j]:
                        out.push(buf_x[j], buf_y[j])
                for j in range(buf_pos + 1):
                    if min_x == buf_x[j] and min_y != buf_y[j]:
                        out.push(buf_x[j], buf_y[j])
        buf_pos += 1
        if buf_pos == W:
            buf_pos = 0
    if min_x != NONE:
        out.push(min_x, min_y)
    return out^


struct ReferenceIndex(Movable):
    """Pre-built minimizer index of a reference sequence, reusable read-only across many queries and threads."""

    var length: Int
    var slot_of: Dict[UInt64, Int]
    var slots: List[List[UInt64]]
    var max_occ: Int

    def __init__(out self):
        self.length = 0
        self.slot_of = Dict[UInt64, Int]()
        self.slots = List[List[UInt64]]()
        self.max_occ = MIN_MID_OCC

    def __init__(out self, reference: List[UInt8]):
        self.length = len(reference)
        var r = sketch(reference)
        self.slot_of = Dict[UInt64, Int]()
        self.slots = List[List[UInt64]]()
        for m in range(len(r.x)):
            var h = r.x[m] >> 8
            var at = self.slot_of.get(h)
            if at:
                self.slots[at.value()].append(r.y[m])
            else:
                self.slot_of[h] = len(self.slots)
                var hits: List[UInt64] = [r.y[m]]
                self.slots.append(hits^)
        var counts = List[Int](capacity=len(self.slots))
        for s in self.slots:
            counts.append(len(s))
        sort(counts)
        var max_occ = MIN_MID_OCC
        if len(counts) > 0:
            max_occ = counts[min(Int((1.0 - MID_OCC_FRAC) * Float64(len(counts))), len(counts) - 1)] + 1
        self.max_occ = min(max(max_occ, MIN_MID_OCC), MAX_MID_OCC)


comptime _DIGIT_BITS = 11
comptime _INSERTION = 32
"""Fewer keys than this are insertion sorted."""


def radix_sort(mut keys: List[UInt64]):
    """`keys` sorted (all keys distinct here, so any sort gives this order)."""
    var scratch = List[UInt64](unsafe_uninit_length=len(keys))
    radix_sort(keys.unsafe_ptr(), scratch.unsafe_ptr(), len(keys))


def radix_sort(keys: MutPointer[UInt64, _], scratch: MutPointer[UInt64, _], n: Int):
    """LSD radix sort of `keys[0:n]` (sporos's); `scratch` holds `n` keys.
    Digits every key shares are skipped."""
    if n <= 1:
        return
    if n <= _INSERTION:
        for i in range(1, n):
            var k = keys[unsafe_offset=i]
            var at = i
            while at > 0 and keys[unsafe_offset=at - 1] > k:
                keys[unsafe_offset=at] = keys[unsafe_offset=at - 1]
                at -= 1
            keys[unsafe_offset=at] = k
        return
    var differ = UInt64(0)
    for i in range(1, n):
        differ |= keys[unsafe_offset=i] ^ keys[unsafe_offset=0]
    comptime digits = 1 << _DIGIT_BITS
    var counts = Array[Int, digits](fill=0)
    var src = keys.unsafe_origin_cast[MutUntrackedOrigin]()
    var dst = scratch.unsafe_origin_cast[MutUntrackedOrigin]()
    var in_scratch = False
    for shift in range(0, 64, _DIGIT_BITS):
        var bits = min(_DIGIT_BITS, 64 - shift)
        var mask = UInt64((1 << bits) - 1)
        if (differ >> UInt64(shift)) & mask == 0:
            continue
        for d in range(1 << bits):
            counts[d] = 0
        for i in range(n):
            counts[Int((src[unsafe_offset=i] >> UInt64(shift)) & mask)] += 1
        var at = 0
        for d in range(1 << bits):
            var c = counts[d]
            counts[d] = at
            at += c
        for i in range(n):
            var k = src[unsafe_offset=i]
            var d = Int((k >> UInt64(shift)) & mask)
            dst[unsafe_offset=counts[d]] = k
            counts[d] += 1
        var t = src
        src = dst
        dst = t
        in_scratch = not in_scratch
    if in_scratch:
        for i in range(n):
            dst[unsafe_offset=i] = src[unsafe_offset=i]


def _key(a: Seed) -> UInt64:
    """The sort key of minimap2's (x, y) order: reference end, then the tandem flag, the span, the query end."""
    var tandem = UInt64(1) if a.flags & TANDEM else UInt64(0)
    return (UInt64(a.x) << 34) | (tandem << 33) | (UInt64(a.span) << UInt64(POSITION_BITS)) | UInt64(a.q)


def _seed_of(key: UInt64) -> Seed:
    var mask = (UInt64(1) << UInt64(POSITION_BITS)) - 1
    return Seed(
        Int(key >> 34),
        Int(key & mask),
        Int((key >> UInt64(POSITION_BITS)) & 0xFF),
        TANDEM if (key >> 33) & 1 else 0,
    )


def sorted_seeds(seeds: List[Seed]) -> List[Seed]:
    var keys = List[UInt64](capacity=len(seeds))
    for a in seeds:
        keys.append(_key(a))
    radix_sort(keys)
    var out = List[Seed](capacity=len(keys))
    for key in keys:
        out.append(_seed_of(key))
    return out^


@fieldwise_init
struct _Hits(Copyable, Movable):
    """A query minimizer's reference hits: its end and strand (`q_pos`), span, the index of its hit slot, the
    tandem flag and mm_seed_select's filter."""

    var q_pos: UInt64
    var span: Int
    var slot: Int
    var count: Int
    var tandem: Bool
    var flt: Bool


def _seed_select(mut seeds: List[_Hits], qlen: Int, max_occ: Int):
    """mm_seed_select: in each streak of seeds past `max_occ`, up to its length / OCC_DIST of the rarest kept."""
    var n = len(seeds)
    if n <= 1:
        return
    var high = False
    for s in seeds:
        if s.count > max_occ:
            high = True
    if not high:
        return
    var last0 = -1
    for i in range(n + 1):
        if i == n or seeds[i].count <= max_occ:
            if i - last0 > 1:
                var ps = 0 if last0 < 0 else Int(seeds[last0].q_pos >> 1)
                var pe = qlen if i == n else Int(seeds[i].q_pos >> 1)
                var st = last0 + 1
                var en = i
                var max_high = Int(Float64(pe - ps) / Float64(OCC_DIST) + 0.499)
                if max_high > 0:
                    var order = List[UInt64]()
                    for j in range(st, en):
                        order.append((UInt64(seeds[j].count) << 32) | UInt64(j))
                    radix_sort(order)
                    for m in range(min(min(max_high, 128), len(order))):
                        seeds[Int(order[m] & 0xFFFFFFFF)].flt = True
                for j in range(st, en):
                    seeds[j].flt = not seeds[j].flt
                for j in range(st, en):
                    if seeds[j].count > MAX_MAX_OCC:
                        seeds[j].flt = True
            last0 = i


def collect_anchors(index: ReferenceIndex, query: List[UInt8]) -> List[Seed]:
    """The forward-strand seeds of `query` against pre-indexed `index`, minimap2's filters applied, in its order."""
    var max_occ = index.max_occ
    var mv = sketch(query)
    if len(mv.x) > max_occ:
        var in_query = Dict[UInt64, Int]()
        for x in mv.x:
            in_query[x] = in_query.get(x, 0) + 1
        var kept = Sketch()
        for m in range(len(mv.x)):
            var c = in_query.get(mv.x[m], 0)
            if not (c > max_occ and Float64(c) > Float64(len(mv.x)) * Q_OCC_FRAC):
                kept.push(mv.x[m], mv.y[m])
        mv = kept^
    var seeds = List[_Hits]()
    for i in range(len(mv.x)):
        var h = mv.x[i] >> 8
        var at = index.slot_of.get(h)
        if not at:
            continue
        var slot = at.value()
        var tandem = (i > 0 and mv.x[i - 1] >> 8 == h) or (i + 1 < len(mv.x) and mv.x[i + 1] >> 8 == h)
        seeds.append(_Hits(mv.y[i], Int(mv.x[i] & 0xFF), slot, len(index.slots[slot]), tandem, False))
    if MAX_MAX_OCC > max_occ:
        _seed_select(seeds, len(query), max_occ)
    else:
        for ref s in seeds:
            s.flt = s.count > max_occ
    var out = List[Seed]()
    for s in seeds:
        if s.flt:
            continue
        for hit in index.slots[s.slot]:
            if (hit & 1) == (s.q_pos & 1):
                out.append(Seed(Int(hit >> 1), Int(s.q_pos >> 1), s.span, TANDEM if s.tandem else 0))
    return sorted_seeds(out)


def collect_anchors(reference: List[UInt8], query: List[UInt8]) -> List[Seed]:
    """The forward-strand seeds of `query` against `reference`, minimap2's filters applied, in its order."""
    var index = ReferenceIndex(reference)
    return collect_anchors(index, query)
