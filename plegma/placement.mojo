"""What `clustered_consensus` and `phased_consensus` give: clusters of reads, each with its consensus made as
`consensus` makes one, and every read placed in its cluster's consensus (`Placement`, a CIGAR at the scoring's
costs, `alignment_costs`)."""

from karan import Aligner, Alignment, Costs, Mode, Ties

from .consensus import ConsensusOptions, consensus
from .homopolymer import HomopolymerModel
from .kernel import INSIDE, LEFT, RIGHT, Scoring


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
    """`clustered_consensus`'s outcome: per cluster its `consensus`, its `members` (read indices, in input order)
    and whether it is `joined`; per read its `placements` entry."""

    var consensus: List[List[UInt8]]
    var members: List[List[Int]]
    var joined: List[Bool]
    """False where the cluster's reads made two sides that could not be joined (no read spanning, `consensus`'s two
    parts): its consensus is then the longer side alone, and the other side's reads' placements in it mean nothing."""
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
    mut joined: Bool,
) raises -> List[UInt8]:
    """The consensus of `reads` alone; of two sides not joined (no read spanning), the longer, and `joined` False."""
    var mine = List[List[UInt8]](capacity=len(reads))
    var mine_anchors = List[Int](capacity=len(reads))
    var mine_strands = List[Int](capacity=len(reads))
    for r in reads:
        mine.append(seqs[r].copy())
        mine_anchors.append(anchors[r])
        mine_strands.append(strands[r])
    var parts = consensus(mine^, mine_anchors, mine_strands, options, homopolymers)
    joined = len(parts) == 1
    if len(parts) > 1 and len(parts[1]) > len(parts[0]):
        return parts[1].copy()
    return parts[0].copy()


def _placed(
    mut aligner: Aligner,
    var made: List[List[UInt8]],
    var joined: List[Bool],
    seqs: List[List[UInt8]],
    cluster_of: List[Int],
    anchors: List[Int],
    costs: Costs,
) raises -> Clustered:
    """The clusters' consensus sequences `made` (each `joined` or not) with every read placed: in its cluster (`cluster_of`), or with none
    (-1) in the one it aligns to most cheaply."""
    var out = Clustered()
    out.consensus = made^
    out.joined = joined^
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
