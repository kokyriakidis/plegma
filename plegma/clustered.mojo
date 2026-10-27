"""Several consensus sequences of one read set: the reads clustered into haplotypes (abPOA's `-d`, as a diploid
caller asks for), each cluster's consensus polished against its own reads, and every read aligned to its cluster's.

The clusters start as `cluster_rows` splits the spanning reads by the MSA columns where they disagree; then, as
k-means does, each read goes to the consensus it aligns to more cheaply and each cluster's consensus is made again
from its reads, until no read moves. Two clusters stay apart only while each holds at least `min_freq` of the reads
and their consensus sequences differ by more than one letter in a homopolymer run (`_apart_by_runs_alone`: the
reads' length errors there split one haplotype's reads as readily as two haplotypes'). A read covering neither anchor (`INSIDE`) builds no consensus; it is aligned
inside the one it fits best.
"""

from std.math import log10

from karan import Aligner, Alignment, Costs, Mode, Ties

from .clusters import cluster_rows
from .consensus import ConsensusOptions, consensus
from .graph import Graph
from .alphabet import encode_base
from .homopolymer import FLANK_UNKNOWN, STRAND_UNKNOWN, HomopolymerModel
from .kernel import BOTH, INSIDE, LEFT, RIGHT, Scoring
from .poa import poa


comptime MAX_REFINEMENTS = 4
"""Rounds of moving reads to their nearest consensus, at most."""


@fieldwise_init
struct Placement(Copyable, Movable):
    """A read's alignment to its cluster's consensus: the `cluster`; the CIGAR (`=`, `X`, `I` a letter the read has
    and the consensus not, `D` one the consensus has and the read not), its `cost` under the scoring's costs, and
    the spans it covers, consensus letters `consensus_start ..< consensus_end` and read letters `read_start ..<
    read_end` (a read anchored at one end, or inside, covers part of the consensus)."""

    var cluster: Int
    var cigar: String
    var cost: Int
    var consensus_start: Int
    var consensus_end: Int
    var read_start: Int
    var read_end: Int


struct Clustered(Movable):
    """`clustered_consensus`'s outcome: per cluster its `consensus` and its `members` (read indices, in input
    order); per read its `placements` entry."""

    var consensus: List[List[UInt8]]
    var members: List[List[Int]]
    var placements: List[Placement]

    def __init__(out self):
        self.consensus = List[List[UInt8]]()
        self.members = List[List[Int]]()
        self.placements = List[Placement]()


def alignment_costs(scoring: Scoring) raises -> Costs:
    """`scoring` as costs, twice its scores so an odd match converts exactly: a mismatch costs 2 (mismatch +
    match), a gap letter 2 extend + match over 2 opens. Globally the least costly alignment is the best-scoring one."""
    var per_match = Int(scoring.match_score)
    var mismatch = 2 * (Int(scoring.mismatch) + per_match)
    if scoring.open2 > 0 and scoring.extend2 > 0:
        return Costs.two_piece(
            mismatch,
            2 * Int(scoring.open1),
            2 * Int(scoring.extend1) + per_match,
            2 * Int(scoring.open2),
            2 * Int(scoring.extend2) + per_match,
        )
    return Costs.affine(mismatch, 2 * Int(scoring.open1), 2 * Int(scoring.extend1) + per_match)


def _mode(anchor: Int) -> Mode:
    """How a read anchored so aligns to a consensus: whole (`BOTH`), from its start (`LEFT`, the consensus's end
    free), to its end (`RIGHT`), or anywhere inside it (`INSIDE`)."""
    if anchor == LEFT:
        return Mode.PREFIX
    if anchor == RIGHT:
        return Mode.SUFFIX
    if anchor == INSIDE:
        return Mode.INFIX
    return Mode.GLOBAL


def _aligned(
    mut aligner: Aligner, consensus: List[UInt8], seq: List[UInt8], anchor: Int, costs: Costs
) raises -> Alignment:
    return aligner.align(Span(consensus), Span(seq), costs, _mode(anchor), ties=Ties.LEFT)


def _consensus_of(
    seqs: List[List[UInt8]],
    reads: List[Int],
    anchors: List[Int],
    strands: List[Int],
    options: ConsensusOptions,
    homopolymers: HomopolymerModel,
) raises -> List[UInt8]:
    """The consensus of `reads` alone; of two sides not joined (no read spanning), the longer."""
    var mine = List[List[UInt8]](capacity=len(reads))
    var mine_anchors = List[Int](capacity=len(reads))
    var mine_strands = List[Int](capacity=len(reads))
    for r in reads:
        mine.append(seqs[r].copy())
        mine_anchors.append(anchors[r])
        mine_strands.append(strands[r])
    var parts = consensus(mine^, mine_anchors, mine_strands, options, homopolymers)
    if len(parts) > 1 and len(parts[1]) > len(parts[0]):
        return parts[1].copy()
    return parts[0].copy()


def _apart_by_runs_alone(
    mut aligner: Aligner, a: List[UInt8], b: List[UInt8], noise_run: Int, mut runs: List[Int], mut single: Bool
) raises -> Bool:
    """Whether `a` and `b` differ only by letters added to or dropped from runs of one letter at least `noise_run`
    long (never with `noise_run` 0), each such run's start in `a` then in `runs`, and `single` whether each by one
    letter alone: what reads' length errors in homopolymers make, the most common split of one haplotype's reads in
    two (on HG002 HiFi, 446 of 466 such splits where the truth is homozygous)."""
    runs.clear()
    single = True
    if noise_run <= 0:
        return False
    var found = aligner.align(Span(a), Span(b), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
    if found.cost == 0:
        return False
    var i = 0  # in a
    var j = 0  # in b
    var count = 0
    for c in found.cigar.as_bytes():
        if c >= UInt8(ord("0")) and c <= UInt8(ord("9")):
            count = count * 10 + Int(c - UInt8(ord("0")))
            continue
        if c == UInt8(ord("=")):
            i += count
            j += count
        elif c == UInt8(ord("X")):
            return False
        elif c == UInt8(ord("D")):
            # Letters of `a`'s run `b` lacks: all of the run's letter, the run long enough even in `b`.
            for k in range(count):
                if a[i + k] != a[i]:
                    return False
            if _run_at(a, i) - count < noise_run - 1 or _run_at(a, i) < noise_run:
                return False
            runs.append(_run_start(a, i))
            single = single and count == 1
            i += count
        else:
            for k in range(count):
                if b[j + k] != b[j]:
                    return False
            if _run_at(b, j) < noise_run:
                return False
            # The letters `b` has more: `a`'s run of them beside the gap.
            if i < len(a) and a[i] == b[j]:
                runs.append(_run_start(a, i))
            elif i > 0 and a[i - 1] == b[j]:
                runs.append(_run_start(a, i - 1))
            else:
                return False
            single = single and count == 1
            j += count
        count = 0
    return True


def _runs_in(a: List[UInt8], b: List[UInt8], runs: List[Int]) -> List[Int]:
    """The length in `b` of each run of `a` starting at `runs` (ascending), `a` and `b` apart only at those runs: a
    run of `b` starts where `a`'s does, shifted by the runs before it."""
    var out = List[Int]()
    var shift = 0
    for s in runs:
        var e = s
        while e < len(a) and a[e] == a[s]:
            e += 1
        var f = s + shift
        while f < len(b) and b[f] == a[s]:
            f += 1
        out.append(f - (s + shift))
        shift += out[len(out) - 1] - (e - s)
    return out^


def _with_runs(
    a: List[UInt8], b: List[UInt8], runs: List[Int], shorter: List[Int], longer: List[Int]
) raises -> List[List[UInt8]]:
    """`a` and `b`, apart only at the runs of `a` starting at `runs` (ascending), rebuilt from `a` with each such run
    `shorter[k]` long in the one whose run there is shorter and `longer[k]` in the other."""
    var in_b = _runs_in(a, b, runs)
    var b_longer = List[Bool]()
    for k in range(len(runs)):
        b_longer.append(in_b[k] > _run_at(a, runs[k]))
    var out = List[List[UInt8]]()
    out.append(List[UInt8]())
    out.append(List[UInt8]())
    for side in range(2):
        var at = 0
        for k in range(len(runs)):
            var s = runs[k]
            var e = s
            while e < len(a) and a[e] == a[s]:
                e += 1
            for i in range(at, s):
                out[side].append(a[i])
            var long_here = (side == 1) == b_longer[k]
            for _ in range(longer[k] if long_here else shorter[k]):
                out[side].append(a[s])
            at = e
        for i in range(at, len(a)):
            out[side].append(a[i])
    return out^


def _run_start(seq: List[UInt8], i: Int) -> Int:
    var b = i
    while b > 0 and seq[b - 1] == seq[i]:
        b -= 1
    return b


def _length_in(cigar: String, consensus_start: Int, letters: List[UInt8], s: Int, e: Int, letter: UInt8) -> Int:
    """The length a read reports for the consensus's run `s ..< e` of `letter`, from its CIGAR against the consensus
    (from `consensus_start`): its letters aligned within the run and inserted inside it or at its ends; -1 when it
    does not cover the letters either side of the run, or holds another letter inside it (no vote, as shasta2's)."""
    var i = consensus_start
    var j = 0
    var count = 0
    var length = 0
    var covered_left = False
    var covered_right = False
    for c in cigar.as_bytes():
        if c >= UInt8(ord("0")) and c <= UInt8(ord("9")):
            count = count * 10 + Int(c - UInt8(ord("0")))
            continue
        for _ in range(count):
            if c == UInt8(ord("I")):
                if i >= s and i <= e:
                    if letters[j] == letter:
                        length += 1
                    elif i > s and i < e:
                        return -1
                j += 1
            elif c == UInt8(ord("D")):
                covered_left = covered_left or i == s - 1
                covered_right = covered_right or i == e
                i += 1
            else:
                if i == s - 1:
                    covered_left = True
                if i == e:
                    covered_right = True
                if i >= s and i < e:
                    if letters[j] != letter:
                        return -1
                    length += 1
                i += 1
                j += 1
        count = 0
    return length if covered_left and covered_right else -1


def _two_lengths(
    mut aligner: Aligner,
    consensus: List[UInt8],
    s: Int,
    other: Int,
    seqs: List[List[UInt8]],
    reads: List[Int],
    anchors: List[Int],
    strands: List[Int],
    model: HomopolymerModel,
    mut shorter: Int,
    mut longer: Int,
) raises -> Float64:
    """How much better two true lengths of the run at `s` explain the reads' lengths there than one, under `model`:
    the log likelihood ratio in dB, each read under two from either length with probability a half; the two lengths
    in `shorter` and `longer`, searched from 4 under the shorter of the run's length and `other` (the other
    cluster's) to 4 over the longer. All the reads together, not each cluster's: clusters split by the run's length
    hold the reads biased apart."""
    var e = s
    while e < len(consensus) and consensus[e] == consensus[s]:
        e += 1
    var base = Int(encode_base(consensus[s]))
    if base > 3:
        return -Float64.MAX
    var left = Int(encode_base(consensus[s - 1])) if s > 0 else FLANK_UNKNOWN
    var right = Int(encode_base(consensus[e])) if e < len(consensus) else FLANK_UNKNOWN
    left = left if left <= 3 else FLANK_UNKNOWN
    right = right if right <= 3 else FLANK_UNKNOWN
    var lengths = List[Int]()
    var on = List[Int]()
    for r in reads:
        var found = aligner.align(Span(consensus), Span(seqs[r]), Costs.edit(), _mode(anchors[r]), ties=Ties.LEFT)
        var m = _length_in(found.cigar, found.reference_start, seqs[r], s, e, consensus[s])
        if m >= 0:
            lengths.append(m)
            on.append(strands[r])
    if len(lengths) < 2:
        return -Float64.MAX
    var n0 = e - s
    var low = max(2, min(n0, other) - 4)
    var high = max(n0, other) + 4
    var one = -Float64.MAX
    for n in range(low, high + 1):
        var total = 0.0
        for k in range(len(lengths)):
            total += model.log_probability(base, left, right, n, lengths[k], on[k])
        one = max(one, total)
    var two = -Float64.MAX
    shorter = n0
    longer = n0 + 1
    for n1 in range(low, high + 1):
        for n2 in range(n1 + 1, high + 1):
            var total = 0.0
            for k in range(len(lengths)):
                var x = model.log_probability(base, left, right, n1, lengths[k], on[k])
                var y = model.log_probability(base, left, right, n2, lengths[k], on[k])
                var top = max(x, y)
                total += top + 10.0 * log10(0.5 * 10.0 ** ((x - top) / 10.0) + 0.5 * 10.0 ** ((y - top) / 10.0))
            if total > two:
                two = total
                shorter = n1
                longer = n2
    return two - one


def _run_at(seq: List[UInt8], i: Int) -> Int:
    """The length of the run of one letter holding position `i` of `seq`."""
    var b = i
    var e = i
    while b > 0 and seq[b - 1] == seq[i]:
        b -= 1
    while e + 1 < len(seq) and seq[e + 1] == seq[i]:
        e += 1
    return e - b + 1


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
    var costs = alignment_costs(options.poa.scoring)
    var aligner = Aligner(warm=False)
    var made = List[List[UInt8]]()
    for g in groups:
        made.append(_consensus_of(seqs, g, build_anchor, strand_of, options, homopolymers))
    if len(groups) >= 2:
        var floor = max(1, Int(min_freq * Float64(len(builders)) + 0.999999))
        var label = List[Int](length=n, fill=-1)
        for k in range(len(groups)):
            for r in groups[k]:
                label[r] = k
        for _ in range(MAX_REFINEMENTS):
            var moved = False
            var next = List[List[Int]](length=len(made), fill=List[Int]())
            for r in builders:
                var best = label[r] if label[r] >= 0 else 0
                var best_cost = _aligned(aligner, made[best], seqs[r], build_anchor[r], costs).cost
                for k in range(len(made)):
                    if k == best:
                        continue
                    var cost = _aligned(aligner, made[k], seqs[r], build_anchor[r], costs).cost
                    if cost < best_cost:
                        best = k
                        best_cost = cost
                if best != label[r]:
                    moved = True
                    label[r] = best
                next[best].append(r)
            var kept = List[List[Int]]()
            for g in next:
                if len(g) >= floor:
                    kept.append(g.copy())
            if len(kept) < len(next):
                # A cluster under the floor: its reads are another's, so one cluster fewer.
                groups = kept^
                if len(groups) <= 1:
                    break
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
        # Clusters spelling the same consensus are one; so are clusters apart only by homopolymer lengths, unless a
        # homopolymer model finds two lengths at one of those runs (`options.two_lengths_db`).
        var distinct = True
        var runs = List[Int]()
        var single = True
        for a in range(len(made)):
            for b in range(a + 1, len(made)):
                if made[a] == made[b]:
                    distinct = False
                elif _apart_by_runs_alone(aligner, made[a], made[b], options.noise_run, runs, single):
                    # With a model: kept apart where it finds two lengths at one of the runs, then every differing run
                    # at the two lengths all the reads fit, the cluster with the shorter run the shorter. Without:
                    # apart unless every difference is a single letter.
                    var two = not single
                    if homopolymers.loaded():
                        two = False
                        var short_of = List[Int]()
                        var long_of = List[Int]()
                        var in_b = _runs_in(made[a], made[b], runs)
                        for k in range(len(runs)):
                            var n1 = 0
                            var n2 = 0
                            var gain = _two_lengths(
                                aligner,
                                made[a],
                                runs[k],
                                in_b[k],
                                seqs,
                                builders,
                                build_anchor,
                                strand_of,
                                homopolymers,
                                n1,
                                n2,
                            )
                            two = two or gain > options.two_lengths_db
                            short_of.append(n1)
                            long_of.append(n2)
                        if two:
                            var fitted = _with_runs(made[a], made[b], runs, short_of, long_of)
                            made[a] = fitted[0].copy()
                            made[b] = fitted[1].copy()
                    distinct = distinct and two
        if len(groups) <= 1 or not distinct:
            groups = [builders.copy()]
            made = [_consensus_of(seqs, builders, build_anchor, strand_of, options, homopolymers)]
    var cluster_of = List[Int](length=n, fill=-1)
    for k in range(len(groups)):
        for r in groups[k]:
            cluster_of[r] = k
    return _placed(aligner, made^, seqs, cluster_of, anchor_of, costs)


def _placed(
    mut aligner: Aligner,
    var made: List[List[UInt8]],
    seqs: List[List[UInt8]],
    cluster_of: List[Int],
    anchors: List[Int],
    costs: Costs,
) raises -> Clustered:
    """The clusters' consensus sequences `made` with every read placed: in its cluster (`cluster_of`), or with none
    (-1) in the one it aligns to most cheaply."""
    var out = Clustered()
    out.consensus = made^
    out.members = List[List[Int]](length=len(out.consensus), fill=List[Int]())
    for r in range(len(seqs)):
        var best = cluster_of[r]
        var found: Alignment
        if best >= 0:
            found = _aligned(aligner, out.consensus[best], seqs[r], anchors[r], costs)
        else:
            best = 0
            found = _aligned(aligner, out.consensus[0], seqs[r], anchors[r], costs)
            for k in range(1, len(out.consensus)):
                var other = _aligned(aligner, out.consensus[k], seqs[r], anchors[r], costs)
                if other.cost < found.cost:
                    best = k
                    found = other^
        out.members[best].append(r)
        out.placements.append(
            Placement(
                best,
                found.cigar,
                found.cost,
                found.reference_start,
                found.reference_end,
                found.query_start,
                found.query_end,
            )
        )
    return out^


def _set_runs(a: List[UInt8], runs: List[Int], lengths: List[Int]) -> List[UInt8]:
    """`a` with its runs starting at `runs` (ascending) `lengths` long."""
    var out = List[UInt8](capacity=len(a) + 8)
    var at = 0
    for k in range(len(runs)):
        var s = runs[k]
        var e = s
        while e < len(a) and a[e] == a[s]:
            e += 1
        for i in range(at, s):
            out.append(a[i])
        for _ in range(lengths[k]):
            out.append(a[s])
        at = e
    for i in range(at, len(a)):
        out.append(a[i])
    return out^


def _read_lengths(
    mut aligner: Aligner,
    consensus: List[UInt8],
    s: Int,
    e: Int,
    seqs: List[List[UInt8]],
    reads: List[Int],
    anchors: List[Int],
    mut lengths: List[Int],
) raises:
    """The length each of `reads` reports for `consensus`'s run `s ..< e` (-1 for none, see `_length_in`)."""
    lengths.clear()
    for r in reads:
        var found = aligner.align(Span(consensus), Span(seqs[r]), Costs.edit(), _mode(anchors[r]), ties=Ties.LEFT)
        lengths.append(_length_in(found.cigar, found.reference_start, seqs[r], s, e, consensus[s]))


def _best_length(
    model: HomopolymerModel,
    base: Int,
    left: Int,
    right: Int,
    lengths: List[Int],
    strands: List[Int],
    reads: List[Int],
    low: Int,
    high: Int,
    mut score: Float64,
) -> Int:
    """The run length in `low ..= high` the reads' lengths make most likely under `model` (ties to the longer: reads
    lose letters in runs), its log likelihood in dB in `score`; reads with no length (-1) do not vote."""
    var best = -1
    score = -Float64.MAX
    for n in range(low, high + 1):
        var total = 0.0
        var any = False
        for k in range(len(reads)):
            if lengths[k] < 0:
                continue
            any = True
            total += model.log_probability(base, left, right, n, lengths[k], strands[reads[k]])
        if any and total >= score:
            score = total
            best = n
    return best


def phased_consensus(
    var seqs: List[List[UInt8]],
    haplotypes: List[Int],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> Clustered:
    """The two haplotypes' consensus sequences of reads already phased (`haplotypes`: 0 or 1 per read, -1 for a read
    of neither, placed where it fits best), each made as `consensus` makes one from its reads, and every read's
    alignment to its haplotype's (see `clustered_consensus`). Where the two differ only by the lengths of runs of one
    letter and a homopolymer model is loaded, each such run is decided on both haplotypes' reads together: their own
    lengths for the two only when that explains the reads `options.two_lengths_db` better than one length for both, else that
    one length in both (a run is far more often the same on the two haplotypes than not). Two clusters always, the
    second's consensus the first's when no read is on haplotype 1."""
    var n = len(seqs)
    if len(haplotypes) != n:
        raise Error(t"{len(haplotypes)} haplotypes for {n} sequences")
    if len(anchors) != 0 and len(anchors) != n:
        raise Error(t"{len(anchors)} anchors for {n} sequences")
    if len(strands) != 0 and len(strands) != n:
        raise Error(t"{len(strands)} strands for {n} sequences")
    var anchor_of = anchors.copy() if len(anchors) > 0 else List[Int](length=n, fill=BOTH)
    var strand_of = strands.copy() if len(strands) > 0 else List[Int](length=n, fill=STRAND_UNKNOWN)
    var groups = List[List[Int]](length=2, fill=List[Int]())
    var cluster_of = List[Int](length=n, fill=-1)
    for r in range(n):
        var h = haplotypes[r]
        if h != 0 and h != 1 and h != -1:
            raise Error(t"read {r}: haplotype {h} is not 0, 1 or -1")
        if h >= 0:
            cluster_of[r] = h
            if anchor_of[r] != INSIDE:
                groups[h].append(r)
    var made = List[List[UInt8]]()
    for k in range(2):
        made.append(
            _consensus_of(seqs, groups[k], anchor_of, strand_of, options, homopolymers) if len(groups[k])
            > 0 else (List[UInt8]())
        )
    if len(made[1]) == 0 and len(groups[1]) == 0:
        made[1] = made[0].copy()
    if len(made[0]) == 0 and len(groups[0]) == 0:
        made[0] = made[1].copy()
    var aligner = Aligner(warm=False)
    var runs = List[Int]()
    var single = True
    if (
        homopolymers.loaded()
        and made[0] != made[1]
        and _apart_by_runs_alone(aligner, made[0], made[1], max(options.noise_run, 2), runs, single)
    ):
        var both = groups[0].copy()
        both.extend(groups[1].copy())
        var first = List[Int]()
        var second = List[Int]()
        var lengths0 = List[Int]()
        var lengths1 = List[Int]()
        var lengths = List[Int]()
        var in_1 = _runs_in(made[0], made[1], runs)
        for k in range(len(runs)):
            var s = runs[k]
            var e = s
            while e < len(made[0]) and made[0][e] == made[0][s]:
                e += 1
            var base = Int(encode_base(made[0][s]))
            var left = Int(encode_base(made[0][s - 1])) if s > 0 else FLANK_UNKNOWN
            var right = Int(encode_base(made[0][e])) if e < len(made[0]) else FLANK_UNKNOWN
            left = left if left <= 3 else FLANK_UNKNOWN
            right = right if right <= 3 else FLANK_UNKNOWN
            # From 4 under the shorter of the two haplotypes' lengths to 4 over the longer: they may be far apart.
            var low = max(2, min(e - s, in_1[k]) - 4)
            var high = max(e - s, in_1[k]) + 4
            _read_lengths(aligner, made[0], s, e, seqs, groups[0], anchor_of, lengths0)
            _read_lengths(aligner, made[0], s, e, seqs, groups[1], anchor_of, lengths1)
            lengths = lengths0.copy()
            lengths.extend(lengths1.copy())
            var shared_score = 0.0
            var score0 = 0.0
            var score1 = 0.0
            var shared = _best_length(
                homopolymers, base, left, right, lengths, strand_of, both, low, high, shared_score
            )
            var own0 = _best_length(homopolymers, base, left, right, lengths0, strand_of, groups[0], low, high, score0)
            var own1 = _best_length(homopolymers, base, left, right, lengths1, strand_of, groups[1], low, high, score1)
            if shared < 0 or own0 < 0 or own1 < 0:
                # Too few reads to say: as the two consensus sequences have it.
                first.append(-1)
                second.append(-1)
            elif score0 + score1 - shared_score > options.two_lengths_db:
                first.append(own0)
                second.append(own1)
            else:
                first.append(shared)
                second.append(shared)
        # Both rebuilt from the first (the two differ only at these runs), unless a run had too few reads to say.
        var decided = True
        for k in range(len(runs)):
            decided = decided and first[k] >= 0
        if decided:
            var a = made[0].copy()
            made[0] = _set_runs(a, runs, first)
            made[1] = _set_runs(a, runs, second)
    return _placed(aligner, made^, seqs, cluster_of, anchor_of, alignment_costs(options.poa.scoring))
