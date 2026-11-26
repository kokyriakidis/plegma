"""Several consensus sequences of one read set: the reads clustered into haplotypes (abPOA's `-d`, as a diploid
caller asks for), each cluster's consensus polished against its own reads, and every read aligned to its cluster's."""

from std.atomic import Atomic

from karan import Aligner, Alignment, Costs, Mode
from karan.common import spread, thread_count

from .clusters import GAP, MISSING, cluster_rows
from .consensus import ConsensusOptions, ConsensusScratch
from .graph import Graph
from .homopolymer import STRAND_UNKNOWN, HomopolymerModel, reads_say_one
from .kernel import BOTH, INSIDE, LEFT, RIGHT
from .phased import phased_consensus
from .placement import Clustered, Placement, _aligned, _consensus_of, _placed, alignment_costs
from .poa import poa_with_scratch
from .runs import (
    _align_group_reads,
    _apart_by_runs_alone,
    _run_at,
    _run_pairs,
    _set_runs,
    _two_lengths_from,
)
from .trace import CLUSTERED, Trace

comptime MOVING_COSTS = Costs(5, 2, 1, -1, 0)
"""The refinement's costs of a read against a consensus: a one-letter gap 3, a substitution 5."""
comptime MAX_REFINEMENTS = 4
"""Rounds of moving reads to their nearest consensus, at most."""
comptime MIN_UNIT_SPLIT_FREQ = 0.3
"""Two clusters apart by one unit of a tandem repeat of 2 to 6 letters stay apart only while the smaller holds this
share of the reads."""


def _one_unit_apart(a: List[UInt8], b: List[UInt8]) -> Bool:
    """Whether `a` and `b` differ only by one unit, 2 to 6 letters not all one, of a tandem repeat."""
    var longer = a.copy() if len(a) > len(b) else b.copy()
    var shorter = b.copy() if len(a) > len(b) else a.copy()
    var p = len(longer) - len(shorter)
    if p < 2 or p > 6:
        return False
    var i = 0
    while i < len(shorter) and longer[i] == shorter[i]:
        i += 1
    for k in range(i, len(shorter)):
        if longer[k + p] != shorter[k]:
            return False
    var one_letter = True
    for k in range(i + 1, i + p):
        one_letter = one_letter and longer[k] == longer[i]
    if one_letter:
        return False
    var after = i + 2 * p <= len(longer)
    var before = i >= p
    for k in range(p):
        after = after and longer[i + p + k] == longer[i + k]
        before = before and longer[i - p + k] == longer[i + k]
    return after or before


def _first_split(
    mut scratch: ConsensusScratch,
    seqs: List[List[UInt8]],
    spanning: List[Int],
    anchor: Int,
    max_clusters: Int,
    min_freq: Float64,
    options: ConsensusOptions,
) raises -> List[List[Int]]:
    """The reads `spanning`, all anchored as `anchor` says, split as `cluster_rows` splits their MSA's rows."""
    var names = List[String](capacity=len(spanning))
    var mine = List[List[UInt8]](capacity=len(spanning))
    for r in spanning:
        names.append(String(t"r{r}"))
        mine.append(seqs[r].copy())
    var made = poa_with_scratch(
        names, mine, List[Int](length=len(spanning), fill=anchor), options.poa, Graph(), scratch.poa, scratch.aligner
    )
    var by_path = made.graph.msa()
    var rows = List[String](capacity=len(spanning))
    for k in range(len(spanning)):
        var row = List[UInt8](by_path[made.shown[k]].as_bytes())
        if anchor == LEFT:
            var c = len(row) - 1
            while c >= 0 and row[c] == GAP:
                row[c] = MISSING
                c -= 1
        elif anchor == RIGHT:
            var c = 0
            while c < len(row) and row[c] == GAP:
                row[c] = MISSING
                c += 1
        rows.append(String(unsafe_from_utf8=row))
    var out = List[List[Int]]()
    for part in cluster_rows(rows, max_clusters, min_freq, options.poa.threads):
        var reads = List[Int](capacity=len(part))
        for k in part:
            reads.append(spanning[k])
        out.append(reads^)
    return out^


def clustered_consensus(
    var seqs: List[List[UInt8]],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    max_clusters: Int = 2,
    min_freq: Float64 = 0.25,
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> Clustered:
    """At most `max_clusters` consensus sequences of `seqs`, each from the reads nearest it and made as `consensus`
    makes one, and every read's alignment to its cluster's consensus."""
    var scratch = ConsensusScratch()
    return clustered_consensus(scratch, seqs^, anchors, strands, max_clusters, min_freq, options, homopolymers)


def clustered_consensus(
    mut scratch: ConsensusScratch,
    var seqs: List[List[UInt8]],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    max_clusters: Int = 2,
    min_freq: Float64 = 0.25,
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> Clustered:
    """As `clustered_consensus`, reusing `scratch` across calls."""
    var n = len(seqs)
    if max_clusters < 1:
        raise Error(t"max_clusters {max_clusters}: must be 1 or more")
    if not (min_freq >= 0.0 and min_freq <= 1.0):
        raise Error(t"min_freq {min_freq}: must be within 0 and 1")
    if len(anchors) != 0 and len(anchors) != n:
        raise Error(t"{len(anchors)} anchors for {n} sequences")
    if len(strands) != 0 and len(strands) != n:
        raise Error(t"{len(strands)} strands for {n} sequences")
    var anchor_of = anchors.copy() if len(anchors) > 0 else List[Int](length=n, fill=BOTH)
    var strand_of = strands.copy() if len(strands) > 0 else List[Int](length=n, fill=STRAND_UNKNOWN)
    var out = Clustered()
    if n == 0:
        out.consensus.append(List[UInt8]())
        out.members.append(List[Int]())
        out.joined.append(True)
        return out^
    var builders = List[Int]()
    for r in range(n):
        if anchor_of[r] != INSIDE:
            builders.append(r)
    var build_anchor = anchor_of.copy()
    if len(builders) == 0:
        for r in range(n):
            builders.append(r)
            build_anchor[r] = BOTH
    var groups = List[List[Int]]()
    if max_clusters >= 2:
        var spanning = List[Int]()
        for r in builders:
            if build_anchor[r] == BOTH:
                spanning.append(r)
        if len(spanning) >= 2:
            groups = _first_split(scratch, seqs, spanning, BOTH, max_clusters, min_freq, options)
        elif _one_side(builders, build_anchor) != BOTH:
            groups = _first_split(
                scratch, seqs, builders, _one_side(builders, build_anchor), max_clusters, min_freq, options
            )
    if len(groups) <= 1:
        groups = [builders.copy()]
    var trace = Trace()
    if trace.on(CLUSTERED):
        trace.line(CLUSTERED, String(n, " reads, ", len(builders), " build; first split: ", _sizes(groups)))
    var costs = alignment_costs(options.poa.scoring)
    var moving = MOVING_COSTS
    var threads = options.poa.threads
    var joined = List[Bool]()
    var made = _made(scratch, seqs, groups, build_anchor, strand_of, options, homopolymers, joined)
    if len(groups) >= 2:
        var telling = List[Bool](length=n, fill=True)
        var floor = max(1, Int(min_freq * Float64(len(builders)) + 0.999999))
        var label = List[Int](length=n, fill=-1)
        for k in range(len(groups)):
            for r in groups[k]:
                label[r] = k
        var side = _one_side(builders, build_anchor)
        for refinement in range(MAX_REFINEMENTS):
            var moved = False
            var moves = 0
            var next = List[List[Int]](length=len(made), fill=List[Int]())
            var against = _cut(scratch.aligner, made, side)
            var assign_workers = thread_count(threads, len(builders))
            if assign_workers > 1 and len(builders) >= 2:
                var bests = List[Int](length=len(builders), fill=0)
                var tells = List[Bool](length=len(builders), fill=False)
                var failed = List[Bool](length=len(builders), fill=False)
                var bests_p = bests.unsafe_ptr()
                var tells_p = tells.unsafe_ptr()
                var failed_p = failed.unsafe_ptr()
                var taken = Atomic[Int64](0)
                var n_builders = len(builders)

                def assign_worker(
                    slot: Int,
                ) {
                    mut taken,
                    imm n_builders,
                    imm builders,
                    imm label,
                    imm against,
                    imm seqs,
                    imm build_anchor,
                    imm moving,
                    imm bests_p,
                    imm tells_p,
                    imm failed_p,
                }:
                    var local_aligner = Aligner(warm=False)
                    while True:
                        var idx = Int(taken.fetch_add(1))
                        if idx >= n_builders:
                            return
                        var r = builders[idx]
                        try:
                            var best = label[r] if label[r] >= 0 else 0
                            var best_cost = _aligned(
                                local_aligner, against[best], seqs[r], build_anchor[r], moving
                            ).cost
                            var is_telling = False
                            for k in range(len(against)):
                                if k == best:
                                    continue
                                var found = _aligned(local_aligner, against[k], seqs[r], build_anchor[r], moving)
                                var cost = found.cost
                                if cost != best_cost or not _short_of(found, against[k], against[best]):
                                    is_telling = True
                                if cost < best_cost:
                                    best = k
                                    best_cost = cost
                            bests_p[unsafe_offset=idx] = best
                            tells_p[unsafe_offset=idx] = is_telling
                        except:
                            failed_p[unsafe_offset=idx] = True

                spread(assign_worker, assign_workers, assign_workers)
                for idx in range(len(builders)):
                    var r = builders[idx]
                    var best: Int
                    var is_telling: Bool
                    if failed[idx]:
                        best = label[r] if label[r] >= 0 else 0
                        var best_cost = _aligned(scratch.aligner, against[best], seqs[r], build_anchor[r], moving).cost
                        is_telling = False
                        for k in range(len(made)):
                            if k == best:
                                continue
                            var found = _aligned(scratch.aligner, against[k], seqs[r], build_anchor[r], moving)
                            var cost = found.cost
                            if cost != best_cost or not _short_of(found, against[k], against[best]):
                                is_telling = True
                            if cost < best_cost:
                                best = k
                                best_cost = cost
                    else:
                        best = bests[idx]
                        is_telling = tells[idx]
                    telling[r] = is_telling
                    if best != label[r]:
                        moved = True
                        moves += 1
                        label[r] = best
                    next[best].append(r)
            else:
                for r in builders:
                    var best = label[r] if label[r] >= 0 else 0
                    var best_cost = _aligned(scratch.aligner, against[best], seqs[r], build_anchor[r], moving).cost
                    telling[r] = False
                    for k in range(len(made)):
                        if k == best:
                            continue
                        var found = _aligned(scratch.aligner, against[k], seqs[r], build_anchor[r], moving)
                        var cost = found.cost
                        if cost != best_cost or not _short_of(found, against[k], against[best]):
                            telling[r] = True
                        if cost < best_cost:
                            best = k
                            best_cost = cost
                    if best != label[r]:
                        moved = True
                        moves += 1
                        label[r] = best
                    next[best].append(r)
            if trace.on(CLUSTERED):
                trace.line(
                    CLUSTERED,
                    String(
                        "refinement ", refinement, ": ", moves, " reads moved, ", _sizes(next), " (floor ", floor, ")"
                    ),
                )
            floor = max(1, Int(min_freq * Float64(_telling(builders, telling)) + 0.999999))
            var kept = List[List[Int]]()
            for g in next:
                if _telling(g, telling) >= floor:
                    kept.append(g.copy())
            if len(kept) < len(next):
                groups = kept^
                if len(groups) <= 1:
                    break
                for r in builders:
                    label[r] = -1
                for k in range(len(groups)):
                    for r in groups[k]:
                        label[r] = k
                made = _made(scratch, seqs, groups, build_anchor, strand_of, options, homopolymers, joined)
                continue
            if not moved:
                break
            groups = next^
            made = _made(scratch, seqs, groups, build_anchor, strand_of, options, homopolymers, joined)
        # A read of a cluster dropped in the last round joins the consensus it differs from least.
        if len(groups) >= 2:
            var placed = List[Bool](length=n, fill=False)
            for g in groups:
                for r in g:
                    placed[r] = True
            var unplaced = List[Int]()
            for r in builders:
                if not placed[r]:
                    unplaced.append(r)
            if len(unplaced) > 0:
                var against = _cut(scratch.aligner, made, side)
                var drop_workers = thread_count(threads, len(unplaced))
                if drop_workers > 1 and len(unplaced) >= 2:
                    var bests = List[Int](length=len(unplaced), fill=0)
                    var failed = List[Bool](length=len(unplaced), fill=False)
                    var bests_p = bests.unsafe_ptr()
                    var failed_p = failed.unsafe_ptr()
                    var taken = Atomic[Int64](0)
                    var n_unplaced = len(unplaced)

                    def place_dropped(
                        slot: Int,
                    ) {
                        mut taken,
                        imm n_unplaced,
                        imm unplaced,
                        imm against,
                        imm seqs,
                        imm build_anchor,
                        imm moving,
                        imm bests_p,
                        imm failed_p,
                    }:
                        var local_aligner = Aligner(warm=False)
                        while True:
                            var idx = Int(taken.fetch_add(1))
                            if idx >= n_unplaced:
                                return
                            var r = unplaced[idx]
                            try:
                                var best = 0
                                var best_cost = _aligned(
                                    local_aligner, against[0], seqs[r], build_anchor[r], moving
                                ).cost
                                for k in range(1, len(against)):
                                    var cost = _aligned(
                                        local_aligner, against[k], seqs[r], build_anchor[r], moving
                                    ).cost
                                    if cost < best_cost:
                                        best = k
                                        best_cost = cost
                                bests_p[unsafe_offset=idx] = best
                            except:
                                failed_p[unsafe_offset=idx] = True

                    spread(place_dropped, drop_workers, drop_workers)
                    for idx in range(len(unplaced)):
                        var r = unplaced[idx]
                        var best = bests[idx]
                        if failed[idx]:
                            best = 0
                            var best_cost = _aligned(scratch.aligner, against[0], seqs[r], build_anchor[r], moving).cost
                            for k in range(1, len(made)):
                                var cost = _aligned(scratch.aligner, against[k], seqs[r], build_anchor[r], moving).cost
                                if cost < best_cost:
                                    best = k
                                    best_cost = cost
                        groups[best].append(r)
                else:
                    for r in unplaced:
                        var best = 0
                        var best_cost = _aligned(scratch.aligner, against[0], seqs[r], build_anchor[r], moving).cost
                        for k in range(1, len(made)):
                            var cost = _aligned(scratch.aligner, against[k], seqs[r], build_anchor[r], moving).cost
                            if cost < best_cost:
                                best = k
                                best_cost = cost
                        groups[best].append(r)
        var root = List[Int]()
        for k in range(len(groups)):
            root.append(k)
        var runs = List[Int]()
        var single = True
        for a in range(len(groups)):
            for b in range(a + 1, len(groups)):
                var cut = _cut(scratch.aligner, [made[a].copy(), made[b].copy()], side)
                if cut[0] == cut[1]:
                    _join(root, a, b)
                    trace.line(CLUSTERED, String("clusters ", a, " and ", b, " one: the same consensus"))
                elif _one_unit_apart(cut[0], cut[1]) and Float64(
                    min(_telling(groups[a], telling), _telling(groups[b], telling))
                ) < MIN_UNIT_SPLIT_FREQ * Float64(_telling(builders, telling)):
                    _join(root, a, b)
                    trace.line(
                        CLUSTERED,
                        String(
                            "clusters ",
                            a,
                            " and ",
                            b,
                            " one: a repeat unit apart, the smaller under ",
                            MIN_UNIT_SPLIT_FREQ,
                        ),
                    )
                else:
                    if not homopolymers.loaded():
                        if (
                            _apart_by_runs_alone(scratch.aligner, cut[0], cut[1], options.noise_run, runs, single)
                            and single
                        ):
                            _join(root, a, b)
                            trace.line(
                                CLUSTERED,
                                String(
                                    "clusters ",
                                    a,
                                    " and ",
                                    b,
                                    " one: apart by single letters in runs of ",
                                    options.noise_run,
                                    " or more",
                                ),
                            )
                        else:
                            trace.line(CLUSTERED, String("clusters ", a, " and ", b, " apart"))
                        continue
                    var runs_a = List[Int]()
                    var runs_b = List[Int]()
                    var only_runs = _run_pairs(scratch.aligner, cut[0], cut[1], options.noise_run, runs_a, runs_b)
                    var set_a = List[Int]()
                    var set_b = List[Int]()
                    var two = False
                    if len(runs_a) > 0:
                        var pair = List[List[UInt8]]()
                        pair.append(cut[0].copy())
                        pair.append(cut[1].copy())
                        var pair_groups = List[List[Int]]()
                        pair_groups.append(groups[a].copy())
                        pair_groups.append(groups[b].copy())
                        # Align each cluster's reads once against its consensus and reuse across all differing runs.
                        var aligned_a = _align_group_reads(
                            scratch.aligner, pair[0], seqs, pair_groups[0], build_anchor, threads
                        )
                        var aligned_b = _align_group_reads(
                            scratch.aligner, pair[1], seqs, pair_groups[1], build_anchor, threads
                        )
                        for k in range(len(runs_a)):
                            var len_a = _run_at(cut[0], runs_a[k])
                            var len_b = _run_at(cut[1], runs_b[k])
                            var n1 = 0
                            var n2 = 0
                            var one = -1
                            var said = List[List[Int]]()
                            var starts: List[Int] = [runs_a[k], runs_b[k]]
                            var gain = _two_lengths_from(
                                pair,
                                starts,
                                pair_groups,
                                aligned_a,
                                aligned_b,
                                seqs,
                                strand_of,
                                homopolymers,
                                n1,
                                n2,
                                one,
                                said,
                            )
                            if gain > options.two_lengths_db:
                                two = True
                                set_a.append(n1 if len_a < len_b else n2)
                                set_b.append(n2 if len_a < len_b else n1)
                            elif one > 0:
                                set_a.append(one)
                                set_b.append(one)
                            else:
                                set_a.append(len_a)
                                set_b.append(len_b)
                            if len(said) == 2 and reads_say_one(said[0]):
                                set_a[k] = 1
                            if len(said) == 2 and reads_say_one(said[1]):
                                set_b[k] = 1
                    made[a] = _set_runs(
                        made[a].copy(), _shifted(runs_a, len(made[a]) - len(cut[0]) if side == RIGHT else 0), set_a
                    )
                    made[b] = _set_runs(
                        made[b].copy(), _shifted(runs_b, len(made[b]) - len(cut[1]) if side == RIGHT else 0), set_b
                    )
                    cut = _cut(scratch.aligner, [made[a].copy(), made[b].copy()], side)
                    if (only_runs and not two) or cut[0] == cut[1]:
                        _join(root, a, b)
                    if trace.on(CLUSTERED):
                        trace.line(
                            CLUSTERED,
                            String(
                                "clusters ",
                                a,
                                " and ",
                                b,
                                ": ",
                                len(runs_a),
                                " runs apart, ",
                                "two lengths at one" if two else "one length at each",
                                "; one" if _root(root, a) == _root(root, b) else "; apart",
                            ),
                        )
        var merged = List[List[Int]]()
        var made_merged = List[List[UInt8]]()
        var joined_merged = List[Bool]()
        for k in range(len(groups)):
            if _root(root, k) != k:
                continue
            var members = List[Int]()
            var parts = 0
            for j in range(len(groups)):
                if _root(root, j) == k:
                    members.extend(groups[j].copy())
                    parts += 1
            var whole = joined[k]
            made_merged.append(
                made[k].copy() if parts
                == 1 else _consensus_of(scratch, seqs, members, build_anchor, strand_of, options, homopolymers, whole)
            )
            joined_merged.append(whole)
            merged.append(members^)
        if len(groups) <= 1 or len(merged) <= 1:
            groups = [builders.copy()]
            made = _made(scratch, seqs, groups, build_anchor, strand_of, options, homopolymers, joined)
        else:
            groups = merged^
            made = made_merged^
            joined = joined_merged^
    var cluster_of = List[Int](length=n, fill=-1)
    for k in range(len(groups)):
        for r in groups[k]:
            cluster_of[r] = k
    if trace.on(CLUSTERED):
        var lengths = String()
        for k in range(len(made)):
            lengths.write(", " if k > 0 else "", len(made[k]))
        trace.line(CLUSTERED, String(len(groups), " clusters of ", _sizes(groups), " reads, ", lengths, " letters"))
    return _placed(scratch.aligner, made^, joined^, seqs, cluster_of, anchor_of, costs, threads)


def _one_side(builders: List[Int], anchor_of: List[Int]) -> Int:
    """`LEFT` or `RIGHT` when every one of `builders` is anchored there alone, else `BOTH`."""
    var side = anchor_of[builders[0]] if len(builders) > 0 else BOTH
    for r in builders:
        if anchor_of[r] != side:
            return BOTH
    return side if side == LEFT or side == RIGHT else BOTH


def _cut(mut aligner: Aligner, made: List[List[UInt8]], side: Int) raises -> List[List[UInt8]]:
    """`made` as reads anchored at one `side` alone are weighed against it."""
    if side == BOTH or len(made) == 0:
        return made.copy()
    var shortest = 0
    for k in range(len(made)):
        if len(made[k]) < len(made[shortest]):
            shortest = k
    var out = List[List[UInt8]](capacity=len(made))
    for m in made:
        var found = aligner.align(
            Span(m), Span(made[shortest]), Costs.edit(), Mode.PREFIX if side == LEFT else Mode.SUFFIX
        )
        out.append(List[UInt8](Span(m)[0 : found.reference_end] if side == LEFT else Span(m)[found.reference_start :]))
    return out^


def _shifted(at: List[Int], by: Int) -> List[Int]:
    var out = at.copy()
    for k in range(len(out)):
        out[k] += by
    return out^


def _made(
    mut scratch: ConsensusScratch,
    seqs: List[List[UInt8]],
    groups: List[List[Int]],
    build_anchor: List[Int],
    strand_of: List[Int],
    options: ConsensusOptions,
    homopolymers: HomopolymerModel,
    mut joined: List[Bool],
) raises -> List[List[UInt8]]:
    """Each group's consensus (`_consensus_of`), in parallel across groups when `options.poa.threads > 1`."""
    var k_count = len(groups)
    var out = List[List[UInt8]](length=k_count, fill=List[UInt8]())
    joined = List[Bool](length=k_count, fill=True)
    var group_workers = thread_count(options.poa.threads, k_count)
    if group_workers > 1 and k_count >= 2:
        var sub_options = options.copy()
        sub_options.poa.threads = max(1, options.poa.threads // group_workers)
        var failed = List[Bool](length=k_count, fill=False)
        var out_p = out.unsafe_ptr()
        var joined_p = joined.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def build_cluster(
            slot: Int,
        ) {
            mut taken,
            imm k_count,
            imm seqs,
            imm groups,
            imm build_anchor,
            imm strand_of,
            imm sub_options,
            imm homopolymers,
            imm out_p,
            imm joined_p,
            imm failed_p,
        }:
            var local_scratch = ConsensusScratch()
            while True:
                var k = Int(taken.fetch_add(1))
                if k >= k_count:
                    return
                var whole = True
                try:
                    out_p[unsafe_offset=k] = _consensus_of(
                        local_scratch, seqs, groups[k], build_anchor, strand_of, sub_options, homopolymers, whole
                    )
                    joined_p[unsafe_offset=k] = whole
                except:
                    failed_p[unsafe_offset=k] = True

        spread(build_cluster, group_workers, group_workers)
        for k in range(k_count):
            if failed[k]:
                var whole = True
                out[k] = _consensus_of(scratch, seqs, groups[k], build_anchor, strand_of, options, homopolymers, whole)
                joined[k] = whole
        return out^

    for k in range(k_count):
        var whole = True
        out[k] = _consensus_of(scratch, seqs, groups[k], build_anchor, strand_of, options, homopolymers, whole)
        joined[k] = whole
    return out^


def _short_of(found: Alignment, consensus: List[UInt8], other: List[UInt8]) -> Bool:
    """Whether a read aligned to `consensus` (`found`) stops short of the letters where `other` differs from it."""
    var shortest = min(len(consensus), len(other))
    var prefix = 0
    while prefix < shortest and consensus[prefix] == other[prefix]:
        prefix += 1
    var suffix = 0
    while suffix < shortest and consensus[len(consensus) - 1 - suffix] == other[len(other) - 1 - suffix]:
        suffix += 1
    var after_suffix = len(consensus) - suffix
    return found.reference_start > min(prefix, after_suffix) or found.reference_end < max(prefix, after_suffix)


def _telling(reads: List[Int], telling: List[Bool]) -> Int:
    var n = 0
    for r in reads:
        if telling[r]:
            n += 1
    return n


def _sizes(groups: List[List[Int]]) -> String:
    var out = String("[")
    for k in range(len(groups)):
        out.write(", " if k > 0 else "", len(groups[k]))
    out.write("]")
    return out^


def _root(root: List[Int], k: Int) -> Int:
    var r = k
    while root[r] != r:
        r = root[r]
    return r


def _join(mut root: List[Int], a: Int, b: Int):
    var ra = _root(root, a)
    var rb = _root(root, b)
    if ra != rb:
        root[max(ra, rb)] = min(ra, rb)
