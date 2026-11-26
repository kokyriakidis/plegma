"""The two haplotypes' consensus sequences of reads already phased (`phased_consensus`): each made from its labelled
reads, then again with the unlabelled reads on the haplotype whose letters they carry where the two differ
(`_sides`, at the places `_differences` finds), and each homopolymer run they differ by decided on both haplotypes'
reads together (`runs`)."""

from std.atomic import Atomic

from karan import Aligner, Alignment, Costs, Mode, Ties
from karan.common import spread, thread_count

from .alphabet import encode_base
from .cigar import ALIGN_MATCH, INSERTION, MATCH, haplotype_at, query_at
from .consensus import ConsensusOptions, ConsensusScratch
from .homopolymer import STRAND_UNKNOWN, HomopolymerModel, reads_say_one
from .kernel import BOTH, INSIDE
from .placement import Clustered, _consensus_of, _placed, alignment_costs
from .runs import (
    _align_group_reads,
    _best_length,
    _best_shared_length,
    _flanks,
    _read_lengths_from,
    _run_pairs,
    _set_runs,
)
from .trace import PHASED, Trace

# Re-exported for callers/tests that import `_haplotype_at` from here.
comptime _haplotype_at = haplotype_at

comptime TO_BOTH = 2
"""`_sides`: a read for both haplotypes, covering no place where the two differ."""
comptime MIN_CARRIERS = 1
comptime MIN_CARRIER_SHARE = 0.25
comptime PHASED_CARRIERS = 3
comptime PHASED_AGREEMENT = 0.75
comptime SITE_FLANK = 8


@fieldwise_init
struct _Site(Copyable, Movable):
    """A place where two haplotypes differ, on the longer's letters `start ..< end` (`SITE_FLANK` either side
    included) and the shorter's `short_start ..< short_end`: the longer spells it `longer`, the shorter `shorter`."""

    var start: Int
    var end: Int
    var short_start: Int
    var short_end: Int
    var longer: List[UInt8]
    var shorter: List[UInt8]


def _differences(mut aligner: Aligner, longer: List[UInt8], shorter: List[UInt8], costs: Costs) raises -> List[_Site]:
    """Where `shorter`, placed inside `longer` (`Mode.INFIX`), spells it otherwise, each place widened by
    `SITE_FLANK` and those that then touch joined."""
    var out = List[_Site]()
    if len(longer) == 0 or len(shorter) == 0 or longer == shorter:
        return out^
    var placed = aligner.align(Span(longer), Span(shorter), costs, Mode.INFIX, ties=Ties.LEFT)
    var at = query_at(placed)
    var spans = List[Tuple[Int, Int]]()
    var r = placed.reference_start
    var run_n = 0
    for ch in placed.cigar.as_bytes():
        if ch >= 48 and ch <= 57:
            run_n = run_n * 10 + Int(ch - 48)
            continue
        var length = run_n if run_n > 0 else 1
        run_n = 0
        if ch == MATCH or ch == ALIGN_MATCH:
            r += length
            continue
        var consumed = 0 if ch == INSERTION else length
        var start = max(placed.reference_start, r - SITE_FLANK)
        var end = min(placed.reference_end, r + consumed + SITE_FLANK)
        if len(spans) > 0 and start <= spans[len(spans) - 1][1]:
            spans[len(spans) - 1] = (spans[len(spans) - 1][0], max(end, spans[len(spans) - 1][1]))
        else:
            spans.append((start, end))
        r += consumed
    for ref span in spans:
        var long_letters = List[UInt8](capacity=span[1] - span[0])
        for i in range(span[0], span[1]):
            long_letters.append(longer[i])
        var short_start = at[span[0] - placed.reference_start]
        var short_end = at[span[1] - placed.reference_start]
        var short_letters = List[UInt8](capacity=short_end - short_start)
        for i in range(short_start, short_end):
            short_letters.append(shorter[i])
        out.append(_Site(span[0], span[1], short_start, short_end, long_letters^, short_letters^))
    return out^


def _spelling(seq: List[UInt8], start: Int, end: Int, at: List[Int], lo: Int, hi: Int, site: _Site) -> Int:
    """Which spelling of `site` a read carries: 0 the longer's, 1 the shorter's, -1 neither."""
    if lo <= start or hi >= end:
        return -1
    var q0 = at[lo - start]
    var q1 = at[hi - start]
    var span_len = q1 - q0
    if span_len == len(site.longer):
        var eq = True
        for k in range(span_len):
            if seq[q0 + k] != site.longer[k]:
                eq = False
                break
        if eq:
            return 0
    if span_len == len(site.shorter):
        var eq = True
        for k in range(span_len):
            if seq[q0 + k] != site.shorter[k]:
                eq = False
                break
        if eq:
            return 1
    return -1


def _first(a: List[UInt8], b: List[UInt8]) -> Bool:
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            return a[i] < b[i]
    return len(a) <= len(b)


def _spell_one_builder(
    mut aligner: Aligner,
    made: List[List[UInt8]],
    seq: List[UInt8],
    h: Int,
    long_index: Int,
    short_index: Int,
    sites: List[_Site],
    costs: Costs,
) raises -> List[Int]:
    """Aligns one builder read to its target haplotype and scores its spelling at every differing site."""
    var on_shorter = h == short_index
    ref haplotype = made[short_index if on_shorter else long_index]
    var start: Int
    var end: Int
    var at: List[Int]
    if len(seq) > len(haplotype):
        var inside = aligner.align(Span(seq), Span(haplotype), costs, Mode.INFIX, ties=Ties.LEFT)
        at = haplotype_at(inside)
        start = 0
        end = len(haplotype)
    else:
        var placed = aligner.align(Span(haplotype), Span(seq), costs, Mode.INFIX, ties=Ties.LEFT)
        at = query_at(placed)
        start = placed.reference_start
        end = placed.reference_end
    var row = List[Int](length=len(sites), fill=-1)
    for k in range(len(sites)):
        ref site = sites[k]
        if on_shorter:
            row[k] = _spelling(seq, start, end, at, site.short_start, site.short_end, site)
        else:
            row[k] = _spelling(seq, start, end, at, site.start, site.end, site)
    return row^


def _sides(
    mut aligner: Aligner,
    made: List[List[UInt8]],
    seqs: List[List[UInt8]],
    builders: List[Int],
    reads: List[Int],
    costs: Costs,
    labels: List[Int] = List[Int](),
    threads: Int = 1,
) raises -> List[Int]:
    """For each of `reads`, the haplotype of `made` whose letters it carries where the two differ."""
    var out = List[Int](length=len(reads), fill=TO_BOTH)
    var long_index = (
        0 if len(made[0]) > len(made[1]) or (len(made[0]) == len(made[1]) and _first(made[0], made[1])) else 1
    )
    var short_index = 1 - long_index
    ref longer = made[long_index]
    var sites = _differences(aligner, longer, made[short_index], costs)
    if len(sites) == 0:
        return out^
    var n_builders = len(builders)
    var spelt = List[List[Int]](length=n_builders, fill=List[Int]())
    var row_of = List[Int](length=len(seqs), fill=-1)
    for b in range(n_builders):
        row_of[builders[b]] = b
    var workers = thread_count(threads, n_builders)
    if workers > 1 and n_builders >= 2:
        var order = List[UInt64](capacity=n_builders)
        for b in range(n_builders):
            order.append((UInt64(len(seqs[builders[b]])) << 32) | UInt64(b))
        sort(order)
        order.reverse()
        var failed = List[Bool](length=n_builders, fill=False)
        var spelt_p = spelt.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def spell_worker(
            slot: Int,
        ) {
            mut taken,
            imm order,
            imm builders,
            imm labels,
            imm made,
            imm seqs,
            imm long_index,
            imm short_index,
            imm sites,
            imm costs,
            imm spelt_p,
            imm failed_p,
        }:
            var local_aligner = Aligner(warm=False)
            while True:
                var idx = Int(taken.fetch_add(1))
                if idx >= len(order):
                    return
                var b = Int(order[idx] & 0xFFFFFFFF)
                var r = builders[b]
                var h = labels[r] if len(labels) > 0 else -1
                try:
                    spelt_p[unsafe_offset=b] = _spell_one_builder(
                        local_aligner, made, seqs[r], h, long_index, short_index, sites, costs
                    )
                except:
                    failed_p[unsafe_offset=b] = True

        spread(spell_worker, workers, workers)
        for b in range(n_builders):
            if failed[b]:
                var r = builders[b]
                var h = labels[r] if len(labels) > 0 else -1
                spelt[b] = _spell_one_builder(aligner, made, seqs[r], h, long_index, short_index, sites, costs)
    else:
        for b in range(n_builders):
            var r = builders[b]
            var h = labels[r] if len(labels) > 0 else -1
            spelt[b] = _spell_one_builder(aligner, made, seqs[r], h, long_index, short_index, sites, costs)

    var counted = List[Bool](length=len(sites), fill=False)
    for k in range(len(sites)):
        var carriers = [0, 0]
        var own = List[Int](length=4, fill=0)
        for b in range(len(spelt)):
            var c = spelt[b][k]
            if c < 0:
                continue
            carriers[c] += 1
            var h = labels[builders[b]] if len(labels) > 0 else -1
            if h == 0 or h == 1:
                own[2 * h + c] += 1
        var fewer = min(carriers[0], carriers[1])
        counted[k] = fewer >= MIN_CARRIERS and Float64(fewer) >= MIN_CARRIER_SHARE * Float64(carriers[0] + carriers[1])
        if not counted[k]:
            var mine = own[2 * long_index]
            var theirs = own[2 * short_index + 1]
            var mine_all = own[2 * long_index] + own[2 * long_index + 1]
            var theirs_all = own[2 * short_index] + own[2 * short_index + 1]
            counted[k] = (
                mine >= PHASED_CARRIERS
                and theirs >= PHASED_CARRIERS
                and Float64(mine) >= PHASED_AGREEMENT * Float64(mine_all)
                and Float64(theirs) >= PHASED_AGREEMENT * Float64(theirs_all)
            )
    for i in range(len(reads)):
        var votes = [0, 0]
        ref row = spelt[row_of[reads[i]]]
        for k in range(len(sites)):
            if counted[k] and row[k] >= 0:
                votes[long_index if row[k] == 0 else 1 - long_index] += 1
        if votes[0] > votes[1]:
            out[i] = 0
        elif votes[1] > votes[0]:
            out[i] = 1
        elif votes[0] > 0:
            out[i] = -1
    return out^


def _any_builds(anchor_of: List[Int]) -> Bool:
    for a in anchor_of:
        if a != INSIDE:
            return True
    return False


def _build_pair(
    mut scratch: ConsensusScratch,
    seqs: List[List[UInt8]],
    groups: List[List[Int]],
    build_mask: List[Bool],
    build_anchor: List[Int],
    strand_of: List[Int],
    options: ConsensusOptions,
    homopolymers: HomopolymerModel,
    mut made: List[List[UInt8]],
    mut joined: List[Bool],
) raises:
    """Builds whichever of the two haplotypes `build_mask` requests, concurrently when both run and `threads > 1`."""
    if options.poa.threads > 1 and build_mask[0] and build_mask[1]:
        var sub_options = options.copy()
        sub_options.poa.threads = max(1, options.poa.threads // 2)
        var built = List[List[UInt8]](length=2, fill=List[UInt8]())
        var built_joined = List[Bool](length=2, fill=True)
        var failed = List[Bool](length=2, fill=False)
        var built_p = built.unsafe_ptr()
        var joined_p = built_joined.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()

        def build_hap(
            k: Int,
        ) {
            imm seqs,
            imm groups,
            imm build_anchor,
            imm strand_of,
            imm sub_options,
            imm homopolymers,
            imm built_p,
            imm joined_p,
            imm failed_p,
        }:
            var local_scratch = ConsensusScratch()
            var whole = True
            try:
                built_p[unsafe_offset=k] = _consensus_of(
                    local_scratch, seqs, groups[k], build_anchor, strand_of, sub_options, homopolymers, whole
                )
                joined_p[unsafe_offset=k] = whole
            except:
                failed_p[unsafe_offset=k] = True

        spread(build_hap, 2, 2)
        for k in range(2):
            if failed[k]:
                made[k] = _consensus_of(
                    scratch, seqs, groups[k], build_anchor, strand_of, options, homopolymers, joined[k]
                )
            else:
                made[k] = built[k].copy()
                joined[k] = built_joined[k]
        return

    for k in range(2):
        if build_mask[k]:
            made[k] = _consensus_of(scratch, seqs, groups[k], build_anchor, strand_of, options, homopolymers, joined[k])


def phased_consensus(
    var seqs: List[List[UInt8]],
    haplotypes: List[Int],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> Clustered:
    """The two haplotypes' consensus sequences of reads already phased (`haplotypes`: 0 or 1 per read, -1 for a read
    of neither), each made as `consensus` makes one from its reads, and every read's alignment to its haplotype's."""
    var scratch = ConsensusScratch()
    return phased_consensus(scratch, seqs^, haplotypes, anchors, strands, options, homopolymers)


def phased_consensus(
    mut scratch: ConsensusScratch,
    var seqs: List[List[UInt8]],
    haplotypes: List[Int],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> Clustered:
    """As `phased_consensus`, reusing `scratch` across calls."""
    var n = len(seqs)
    if len(haplotypes) != n:
        raise Error(t"{len(haplotypes)} haplotypes for {n} sequences")
    if len(anchors) != 0 and len(anchors) != n:
        raise Error(t"{len(anchors)} anchors for {n} sequences")
    if len(strands) != 0 and len(strands) != n:
        raise Error(t"{len(strands)} strands for {n} sequences")
    var anchor_of = anchors.copy() if len(anchors) > 0 else List[Int](length=n, fill=BOTH)
    var strand_of = strands.copy() if len(strands) > 0 else List[Int](length=n, fill=STRAND_UNKNOWN)
    var build_anchor = anchor_of.copy()
    if INSIDE in anchor_of and not _any_builds(anchor_of):
        for r in range(n):
            build_anchor[r] = BOTH
    var groups = List[List[Int]](length=2, fill=List[Int]())
    var cluster_of = List[Int](length=n, fill=-1)
    for r in range(n):
        var h = haplotypes[r]
        if h != 0 and h != 1 and h != -1:
            raise Error(t"read {r}: haplotype {h} is not 0, 1 or -1")
        if h >= 0:
            cluster_of[r] = h
            if build_anchor[r] != INSIDE:
                groups[h].append(r)
    var labelled = [len(groups[0]) > 0, len(groups[1]) > 0]
    var unlabelled = List[Int]()
    for r in range(n):
        if haplotypes[r] < 0 and build_anchor[r] != INSIDE:
            unlabelled.append(r)
    if not labelled[0] and not labelled[1] and len(unlabelled) > 0:
        groups[0] = unlabelled.copy()
        groups[1] = unlabelled.copy()
        unlabelled.clear()
    var made = List[List[UInt8]](length=2, fill=List[UInt8]())
    var joined: List[Bool] = [True, True]
    var initial_mask: List[Bool] = [len(groups[0]) > 0, len(groups[1]) > 0]
    _build_pair(scratch, seqs, groups, initial_mask, build_anchor, strand_of, options, homopolymers, made, joined)
    if len(groups[1]) == 0:
        made[1] = made[0].copy()
        joined[1] = joined[0]
    if len(groups[0]) == 0:
        made[0] = made[1].copy()
        joined[0] = joined[1]
    var threads = options.poa.threads
    if len(unlabelled) > 0:
        var costs = alignment_costs(options.poa.scoring)
        var builders = groups[0].copy()
        builders.extend(groups[1].copy())
        builders.extend(unlabelled.copy())
        var sides = _sides(scratch.aligner, made, seqs, builders, unlabelled, costs, haplotypes, threads)
        var trace = Trace()
        if trace.on(PHASED):
            var count = [0, 0, 0, 0]
            for side in sides:
                count[min(max(side, 0), TO_BOTH) if side >= 0 else 3] += 1
            trace.line(
                PHASED,
                String(
                    len(groups[0]),
                    " and ",
                    len(groups[1]),
                    " labelled reads build; of ",
                    len(unlabelled),
                    " unlabelled: ",
                    count[0],
                    " to haplotype 0, ",
                    count[1],
                    " to 1, ",
                    count[2],
                    " to both, ",
                    count[3],
                    " to neither",
                ),
            )
        var grew: List[Bool] = [False, False]
        for i in range(len(unlabelled)):
            for k in range(2):
                if sides[i] == k or sides[i] == TO_BOTH:
                    groups[k].append(unlabelled[i])
                    grew[k] = True
        _build_pair(scratch, seqs, groups, grew, build_anchor, strand_of, options, homopolymers, made, joined)
    var runs0 = List[Int]()
    var runs1 = List[Int]()
    if homopolymers.loaded() and made[0] != made[1]:
        _ = _run_pairs(scratch.aligner, made[0], made[1], max(options.noise_run, 2), runs0, runs1)
    if len(runs0) > 0:
        var aligned0 = _align_group_reads(scratch.aligner, made[0], seqs, groups[0], anchor_of, threads)
        var aligned1 = _align_group_reads(scratch.aligner, made[1], seqs, groups[1], anchor_of, threads)
        var first = List[Int]()
        var second = List[Int]()
        var lengths0 = List[Int]()
        var lengths1 = List[Int]()
        for k in range(len(runs0)):
            var s = runs0[k]
            var e = s
            while e < len(made[0]) and made[0][e] == made[0][s]:
                e += 1
            var s1 = runs1[k]
            var e1 = s1
            while e1 < len(made[1]) and made[1][e1] == made[1][s1]:
                e1 += 1
            var base = Int(encode_base(made[0][s]))
            var flanks0 = _flanks(made[0], s, e)
            var flanks1 = _flanks(made[1], s1, e1)
            var low = max(2, min(e - s, e1 - s1) - 4)
            var high = max(e - s, e1 - s1) + 4
            if not homopolymers.covers(base, high):
                first.append(e - s)
                second.append(e1 - s1)
                continue
            _read_lengths_from(aligned0, made[0], s, e, seqs, groups[0], lengths0)
            _read_lengths_from(aligned1, made[1], s1, e1, seqs, groups[1], lengths1)
            var shared_score = 0.0
            var score0 = 0.0
            var score1 = 0.0
            var shared = _best_shared_length(
                homopolymers,
                base,
                flanks0,
                flanks1,
                lengths0,
                lengths1,
                strand_of,
                groups[0],
                groups[1],
                low,
                high,
                shared_score,
            )
            var own0 = _best_length(
                homopolymers, base, flanks0[0], flanks0[1], lengths0, strand_of, groups[0], low, high, score0
            )
            var own1 = _best_length(
                homopolymers, base, flanks1[0], flanks1[1], lengths1, strand_of, groups[1], low, high, score1
            )
            if shared < 0 or own0 < 0 or own1 < 0:
                first.append(e - s)
                second.append(e1 - s1)
                continue
            if score0 + score1 - shared_score > options.two_lengths_db:
                first.append(own0)
                second.append(own1)
            else:
                first.append(shared)
                second.append(shared)
            if reads_say_one(lengths0):
                first[k] = 1
            if reads_say_one(lengths1):
                second[k] = 1
        var a = made[0].copy()
        var b = made[1].copy()
        made[0] = _set_runs(a, runs0, first)
        made[1] = _set_runs(b, runs1, second)
    return _placed(
        scratch.aligner, made^, joined^, seqs, cluster_of, anchor_of, alignment_costs(options.poa.scoring), threads
    )
