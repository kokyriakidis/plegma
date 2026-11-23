"""Several consensus sequences of one read set: the reads clustered into haplotypes (abPOA's `-d`, as a diploid
caller asks for), each cluster's consensus polished against its own reads, and every read aligned to its cluster's.

The clusters start as `cluster_rows` splits the spanning reads by the MSA columns where they disagree; then, as
k-means does, each read goes to the consensus it differs from least (`MOVING_COSTS`: a substitution weighs more than
a letter slipped in a run) and each cluster's consensus is made again
from its reads, until no read moves. Two clusters stay apart only while each holds at least `min_freq` of the reads
and their consensus sequences differ by more than homopolymer run lengths the reads fit as one (`_two_lengths`: the
reads' length errors there split one haplotype's reads as readily as two haplotypes'), nor only by one unit of a
short tandem repeat with the smaller under `MIN_UNIT_SPLIT_FREQ` of the reads. Every run two clusters differ by is
set so, whatever else they differ by: one length in both where the reads fit one. A read covering neither anchor
(`INSIDE`) builds no consensus; it is aligned inside the one it fits best.
"""

from karan import Aligner, Alignment, Costs

from .clusters import cluster_rows
from .consensus import ConsensusOptions, consensus
from .graph import Graph
from .homopolymer import HomopolymerModel, STRAND_UNKNOWN, reads_say_one
from .kernel import BOTH, INSIDE
from .placement import Clustered, _aligned, alignment_costs, _consensus_of, _placed
from .poa import poa
from .trace import CLUSTERED, Trace
from .runs import _apart_by_runs_alone, _run_at, _run_pairs, _set_runs, _two_lengths

# Re-exported for callers that import them from here.
from .phased import phased_consensus
from .placement import Clustered, Placement, alignment_costs


comptime MOVING_COSTS = Costs(5, 2, 1, -1, 0)
"""The refinement's costs of a read against a consensus: a one-letter gap 3, a substitution 5, under the 6 of the
insertion and deletion that would otherwise spell it (no costs weigh a substitution more than two one-letter gaps:
the alignment then writes it as those two)."""


comptime MAX_REFINEMENTS = 4
"""Rounds of moving reads to their nearest consensus, at most."""


comptime MIN_UNIT_SPLIT_FREQ = 0.3
"""Two clusters apart by one unit of a tandem repeat of 2 to 6 letters stay apart only while the smaller holds this
share of the reads: reads' slips in short tandem repeats split one haplotype's reads by a unit. On HG002 HiFi chr12,
18 and 20 (pgphase's regions no phase set held), such splits with the smaller under 30% of the reads were 1 true het
to 25 false; at 30% or more, 67 to 11."""


def _one_unit_apart(a: List[UInt8], b: List[UInt8]) -> Bool:
    """Whether `a` and `b` differ only by one unit, 2 to 6 letters not all one, of a tandem repeat: the longer one's
    extra letters, where the two first differ, repeat the unit beside them."""
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
    seqs: List[List[UInt8]], spanning: List[Int], max_clusters: Int, min_freq: Float64, options: ConsensusOptions
) raises -> List[List[Int]]:
    """The spanning reads split as `cluster_rows` splits their MSA's rows (one cluster when it finds no split)."""
    var names = List[String](capacity=len(spanning))
    var mine = List[List[UInt8]](capacity=len(spanning))
    for r in spanning:
        names.append(String(t"r{r}"))
        mine.append(seqs[r].copy())
    var made = poa(names, mine, List[Int](length=len(spanning), fill=BOTH), options.poa, Graph())
    var by_path = made.graph.msa()
    var rows = List[String](capacity=len(spanning))
    for k in range(len(spanning)):
        rows.append(by_path[made.shown[k]])
    var out = List[List[Int]]()
    for part in cluster_rows(rows, max_clusters, min_freq):
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
    """At most `max_clusters` consensus sequences of `seqs` (each anchored as `anchors` says: `BOTH`, `LEFT`,
    `RIGHT` or `INSIDE`; empty: all `BOTH`; on `strands`' strand for `homopolymers`), each from the reads nearest
    it and made as `consensus` makes one, a cluster kept apart only with at least `min_freq` of the reads that make
    consensus sequences; and every read's alignment to its cluster's consensus (at `options.poa.scoring`'s costs,
    gaps leftmost). With one cluster, or none found, the consensus of every read."""
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
        return out^
    # The reads that make consensus sequences; with every read inside, all of them as if spanning.
    var builders = List[Int]()
    for r in range(n):
        if anchor_of[r] != INSIDE:
            builders.append(r)
    var build_anchor = anchor_of.copy()
    if len(builders) == 0:
        for r in range(n):
            builders.append(r)
            build_anchor[r] = BOTH
    # `consensus` checks the anchors and strands as it is first called.
    var groups = List[List[Int]]()
    if max_clusters >= 2:
        var spanning = List[Int]()
        for r in builders:
            if build_anchor[r] == BOTH:
                spanning.append(r)
        if len(spanning) >= 2:
            groups = _first_split(seqs, spanning, max_clusters, min_freq, options)
    if len(groups) <= 1:
        groups = [builders.copy()]
    var trace = Trace()
    if trace.on(CLUSTERED):
        trace.line(CLUSTERED, String(n, " reads, ", len(builders), " build; first split: ", _sizes(groups)))
    var costs = alignment_costs(options.poa.scoring)
    # A read moves to the consensus it differs from least, a substitution weighing more than a letter slipped in a run
    # (reads slip in runs far more often than they misread a letter), as the first split weighs them. At the
    # scoring's costs a gap letter cost more than a mismatch, so reads moved by their run lengths, not by a SNP's
    # letter: HG002 chr20:56,892,630, G in a poly-T on 28 of 62 reads: the G cluster came to hold one T-run length's
    # reads instead, its consensus lost the G, and the two merged as apart by a run alone (written G on both).
    var moving = MOVING_COSTS
    var aligner = Aligner(warm=False)
    var made = List[List[UInt8]]()
    for g in groups:
        made.append(_consensus_of(seqs, g, build_anchor, strand_of, options, homopolymers))
    if len(groups) >= 2:
        # A read costing the same against every consensus and short of the letters where they differ tells the
        # clusters nothing (a partial read): where their sizes decide (`min_freq`, `MIN_UNIT_SPLIT_FREQ`) it counts
        # for none of them, nor in the reads they are shares of. Kept in the first input row's cluster, which takes
        # it, it let input order alone decide; counted evenly for each, it kept a cluster no read told apart. A read
        # over those letters at the same cost lies between the two (`_short_of`): it stays in its cluster and counts.
        var telling = List[Bool](length=n, fill=True)
        var floor = max(1, Int(min_freq * Float64(len(builders)) + 0.999999))
        var label = List[Int](length=n, fill=-1)
        for k in range(len(groups)):
            for r in groups[k]:
                label[r] = k
        for refinement in range(MAX_REFINEMENTS):
            var moved = False
            var moves = 0
            var next = List[List[Int]](length=len(made), fill=List[Int]())
            for r in builders:
                var best = label[r] if label[r] >= 0 else 0
                var best_cost = _aligned(aligner, made[best], seqs[r], build_anchor[r], moving).cost
                telling[r] = False
                for k in range(len(made)):
                    if k == best:
                        continue
                    var found = _aligned(aligner, made[k], seqs[r], build_anchor[r], moving)
                    var cost = found.cost
                    if cost != best_cost or not _short_of(found, made[k], made[best]):
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
                # A cluster under the floor: its reads are another's, so one cluster fewer.
                groups = kept^
                if len(groups) <= 1:
                    break
                # Every label again: the dropped cluster's reads kept its number, past the clusters left (or another's).
                for r in builders:
                    label[r] = -1
                for k in range(len(groups)):
                    for r in groups[k]:
                        label[r] = k
                made.clear()
                for g in groups:
                    made.append(_consensus_of(seqs, g, build_anchor, strand_of, options, homopolymers))
                continue
            if not moved:
                break
            groups = next^
            made.clear()
            for g in groups:
                made.append(_consensus_of(seqs, g, build_anchor, strand_of, options, homopolymers))
        # A read of a cluster dropped in the last round joins the consensus it differs from least.
        if len(groups) >= 2:
            var placed = List[Bool](length=n, fill=False)
            for g in groups:
                for r in g:
                    placed[r] = True
            for r in builders:
                if placed[r]:
                    continue
                var best = 0
                var best_cost = _aligned(aligner, made[0], seqs[r], build_anchor[r], moving).cost
                for k in range(1, len(made)):
                    var cost = _aligned(aligner, made[k], seqs[r], build_anchor[r], moving).cost
                    if cost < best_cost:
                        best = k
                        best_cost = cost
                groups[best].append(r)
        # Clusters spelling the same consensus are one; so are clusters apart only by homopolymer lengths, unless a
        # homopolymer model finds two lengths at one of those runs (`options.two_lengths_db`).
        # Each pair that is one joins (`root`: union of clusters); with two clusters, as one flag.
        var root = List[Int]()
        for k in range(len(groups)):
            root.append(k)
        var runs = List[Int]()
        var single = True
        # Over the clusters still standing: a refinement that leaves one stops with `made` as it was (the two
        # consensus sequences of the round before), and one cluster has nothing to merge.
        for a in range(len(groups)):
            for b in range(a + 1, len(groups)):
                if made[a] == made[b]:
                    _join(root, a, b)
                    trace.line(CLUSTERED, String("clusters ", a, " and ", b, " one: the same consensus"))
                elif _one_unit_apart(made[a], made[b]) and Float64(
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
                        # Without a model: apart unless every difference is a single letter of a run.
                        if _apart_by_runs_alone(aligner, made[a], made[b], options.noise_run, runs, single) and single:
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
                    # With a model, each run the two differ by, whatever else they differ by: two lengths where all
                    # the reads fit two better than one (`_two_lengths`), the cluster with the shorter run the
                    # shorter, else the one length in both. Clusters apart by runs alone stay apart only where a
                    # run has two. One run's two lengths no longer make every other run two (each was fitted to
                    # the pair it found), nor does another difference spare a run the test.
                    var runs_a = List[Int]()
                    var runs_b = List[Int]()
                    var only_runs = _run_pairs(aligner, made[a], made[b], options.noise_run, runs_a, runs_b)
                    var set_a = List[Int]()
                    var set_b = List[Int]()
                    var two = False
                    var pair = List[List[UInt8]]()
                    pair.append(made[a].copy())
                    pair.append(made[b].copy())
                    var pair_groups = List[List[Int]]()
                    pair_groups.append(groups[a].copy())
                    pair_groups.append(groups[b].copy())
                    for k in range(len(runs_a)):
                        var len_a = _run_at(made[a], runs_a[k])
                        var len_b = _run_at(made[b], runs_b[k])
                        var n1 = 0
                        var n2 = 0
                        var one = -1
                        var said = List[List[Int]]()
                        var starts = List[Int]()
                        starts.append(runs_a[k])
                        starts.append(runs_b[k])
                        var gain = _two_lengths(
                            aligner,
                            pair,
                            starts,
                            pair_groups,
                            seqs,
                            build_anchor,
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
                        # A cluster whose reads report the letter once there has no run, as in `phased_consensus`:
                        # the search starts at 2.
                        if len(said) == 2 and reads_say_one(said[0]):
                            set_a[k] = 1
                        if len(said) == 2 and reads_say_one(said[1]):
                            set_b[k] = 1
                    made[a] = _set_runs(made[a].copy(), runs_a, set_a)
                    made[b] = _set_runs(made[b].copy(), runs_b, set_b)
                    if (only_runs and not two) or made[a] == made[b]:
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
        # The clusters joined made one again from their reads; with every cluster one, the consensus of every read.
        # One flag for all pairs made three clusters one wherever two of them were (a third haplotype's lost).
        var merged = List[List[Int]]()
        var made_merged = List[List[UInt8]]()
        for k in range(len(groups)):
            if _root(root, k) != k:
                continue
            var members = List[Int]()
            var parts = 0
            for j in range(len(groups)):
                if _root(root, j) == k:
                    members.extend(groups[j].copy())
                    parts += 1
            made_merged.append(
                made[k].copy() if parts
                == 1 else _consensus_of(seqs, members, build_anchor, strand_of, options, homopolymers)
            )
            merged.append(members^)
        if len(groups) <= 1 or len(merged) <= 1:
            groups = [builders.copy()]
            made = [_consensus_of(seqs, builders, build_anchor, strand_of, options, homopolymers)]
        else:
            groups = merged^
            made = made_merged^
    var cluster_of = List[Int](length=n, fill=-1)
    for k in range(len(groups)):
        for r in groups[k]:
            cluster_of[r] = k
    if trace.on(CLUSTERED):
        var lengths = String()
        for k in range(len(made)):
            lengths.write(", " if k > 0 else "", len(made[k]))
        trace.line(CLUSTERED, String(len(groups), " clusters of ", _sizes(groups), " reads, ", lengths, " letters"))
    return _placed(aligner, made^, seqs, cluster_of, anchor_of, costs)


def _short_of(found: Alignment, consensus: List[UInt8], other: List[UInt8]) -> Bool:
    """Whether a read aligned to `consensus` (`found`) stops short of the letters where `other` differs from it: those
    between the first letter the two do not share at their start and the last at their end (inside a repeat the two
    overlap, a unit's insertion placeable anywhere in it, and the whole repeat counts). A read over them at the cost
    it has against `other` lies between the two, as a run length halfway between two clusters' does once a
    homopolymer model has moved one consensus's run: HG002 chr20:14,181,385, 19 of 24 reads of an A21 cluster
    between its consensus's A22 and the other's A20, counted for neither, and the A21 haplotype dropped under
    `min_freq`."""
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
    """How many of `reads` tell the clusters apart (`telling`)."""
    var n = 0
    for r in reads:
        if telling[r]:
            n += 1
    return n


def _sizes(groups: List[List[Int]]) -> String:
    """Clusters' sizes, for a trace line: `[12, 9]`."""
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
    """Clusters `a` and `b` one (`root`: each cluster's parent, a root its own)."""
    var ra = _root(root, a)
    var rb = _root(root, b)
    if ra != rb:
        root[max(ra, rb)] = min(ra, rb)
