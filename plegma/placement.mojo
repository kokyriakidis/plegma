"""What `clustered_consensus` and `phased_consensus` give: clusters of reads, each with its consensus made as
`consensus` makes one, and every read placed in its cluster's consensus (`Placement`, a CIGAR at the scoring's
costs, `alignment_costs`)."""

from std.atomic import Atomic

from karan import Aligner, Alignment, Costs, Mode, Ties
from karan.common import spread, thread_count

from .consensus import ConsensusOptions, ConsensusScratch, consensus
from .homopolymer import HomopolymerModel
from .kernel import INSIDE, LEFT, RIGHT, Scoring


@fieldwise_init
struct Placement(Copyable, Movable):
    """A read's alignment to its cluster's consensus: the `cluster`; the CIGAR (`=`, `X`, `I` a letter the read has
    and the consensus not, `D` one the consensus has and the read not), its `cost` under the scoring's costs, and
    the spans it covers, consensus letters `consensus_start ..< consensus_end` and read letters `read_start ..<
    read_end` (a read anchored at one end, or inside, covers part of the consensus); and its `margin`, how much more
    the next cheapest cluster costs it: 0 where its cluster placed it, or where two clusters tie."""

    var cluster: Int
    var cigar: String
    var cost: Int
    var consensus_start: Int
    var consensus_end: Int
    var read_start: Int
    var read_end: Int
    var margin: Int


struct Clustered(Movable):
    """`clustered_consensus`'s outcome: per cluster its `consensus`, its `members` (read indices, in input order)
    and whether it is `joined`; per read its `placements` entry."""

    var consensus: List[List[UInt8]]
    var members: List[List[Int]]
    var joined: List[Bool]
    var placements: List[Placement]

    def __init__(out self):
        self.consensus = List[List[UInt8]]()
        self.members = List[List[Int]]()
        self.joined = List[Bool]()
        self.placements = List[Placement]()


def alignment_costs(scoring: Scoring) raises -> Costs:
    """`scoring` as costs, twice its scores so an odd match converts exactly: a mismatch costs 2 (mismatch +
    match), a gap letter 2 extend + match over 2 opens. Globally the least costly alignment is the best-scoring one."""
    var per_match = Int(scoring.match_score)
    var mismatch = 2 * (Int(scoring.mismatch) + per_match)
    if scoring.convex() and 2 * Int(scoring.extend2) + per_match > 0:
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
    mut scratch: ConsensusScratch,
    seqs: List[List[UInt8]],
    reads: List[Int],
    anchors: List[Int],
    strands: List[Int],
    options: ConsensusOptions,
    homopolymers: HomopolymerModel,
    mut joined: Bool,
) raises -> List[UInt8]:
    """The consensus of `reads` alone, reusing `scratch`; of two sides not joined, the longer, and `joined` False."""
    var mine = List[List[UInt8]](capacity=len(reads))
    var mine_anchors = List[Int](capacity=len(reads))
    var mine_strands = List[Int](capacity=len(reads))
    for r in reads:
        mine.append(seqs[r].copy())
        mine_anchors.append(anchors[r])
        mine_strands.append(strands[r])
    var parts = consensus(scratch, mine^, mine_anchors, mine_strands, options, homopolymers)
    joined = len(parts) == 1
    if len(parts) > 1 and len(parts[1]) > len(parts[0]):
        return parts[1].copy()
    return parts[0].copy()


def _place_one(
    mut aligner: Aligner,
    made: List[List[UInt8]],
    seq: List[UInt8],
    given_cluster: Int,
    anchor: Int,
    costs: Costs,
) raises -> Placement:
    """Places a single read against its assigned cluster or the cheapest cluster, computing its margin."""
    var best = given_cluster
    var found: Alignment
    var margin = 0
    if best >= 0:
        found = _aligned(aligner, made[best], seq, anchor, costs)
    else:
        best = 0
        found = _aligned(aligner, made[0], seq, anchor, costs)
        var second = -1
        for k in range(1, len(made)):
            var other = _aligned(aligner, made[k], seq, anchor, costs)
            if other.cost < found.cost:
                second = found.cost
                best = k
                found = other^
            elif second < 0 or other.cost < second:
                second = other.cost
        if second >= 0:
            margin = second - found.cost
    return Placement(
        best,
        found.cigar,
        found.cost,
        found.reference_start,
        found.reference_end,
        found.query_start,
        found.query_end,
        margin,
    )


def _placed(
    mut aligner: Aligner,
    var made: List[List[UInt8]],
    var joined: List[Bool],
    seqs: List[List[UInt8]],
    cluster_of: List[Int],
    anchors: List[Int],
    costs: Costs,
    threads: Int = 1,
) raises -> Clustered:
    """The clusters' consensus sequences `made` (each `joined` or not) with every read placed, in parallel when
    `threads > 1`."""
    var out = Clustered()
    out.consensus = made^
    out.joined = joined^
    out.members = List[List[Int]](length=len(out.consensus), fill=List[Int]())
    var n = len(seqs)
    var workers = thread_count(threads, n)
    if workers > 1 and n >= 2:
        var order = List[UInt64](capacity=n)
        for r in range(n):
            order.append((UInt64(len(seqs[r])) << 32) | UInt64(r))
        sort(order)
        order.reverse()
        var placements = List[Placement](length=n, fill=Placement(0, String(), 0, 0, 0, 0, 0, 0))
        var failed = List[Bool](length=n, fill=False)
        var placements_p = placements.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)
        ref consensus_list = out.consensus

        def place_worker(
            slot: Int,
        ) {
            mut taken,
            imm order,
            imm consensus_list,
            imm seqs,
            imm cluster_of,
            imm anchors,
            imm costs,
            imm placements_p,
            imm failed_p,
        }:
            var local_aligner = Aligner(warm=False)
            while True:
                var idx = Int(taken.fetch_add(1))
                if idx >= len(order):
                    return
                var r = Int(order[idx] & 0xFFFFFFFF)
                try:
                    placements_p[unsafe_offset=r] = _place_one(
                        local_aligner, consensus_list, seqs[r], cluster_of[r], anchors[r], costs
                    )
                except:
                    failed_p[unsafe_offset=r] = True

        spread(place_worker, workers, workers)
        for r in range(n):
            if failed[r]:
                placements[r] = _place_one(aligner, out.consensus, seqs[r], cluster_of[r], anchors[r], costs)
            out.members[placements[r].cluster].append(r)
        out.placements = placements^
        return out^

    for r in range(n):
        var p = _place_one(aligner, out.consensus, seqs[r], cluster_of[r], anchors[r], costs)
        out.members[p.cluster].append(r)
        out.placements.append(p^)
    return out^
