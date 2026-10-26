"""The consensus of a set of reads as one call: `plegma poa`'s default output (`-r 0`) as a library function, for a
caller that assembles many small read sets (a local assembler's, between two anchors) and would rather not run the
command for each.

`consensus` takes the reads, each anchored as a local assembly's reads are (`BOTH`, `LEFT`, `RIGHT`) and on a strand
(`STRAND_UNKNOWN` when not known), and returns what `plegma poa -r 0` writes: one sequence, polished; or, when no read
spans both anchors and the two sides could not be joined, the left side's and the right side's, unpolished.
"""

from .graph import Graph
from .homopolymer import STRAND_UNKNOWN, HomopolymerModel
from .kernel import BOTH, GLOBAL, LEFT, RIGHT, Topology
from .poa import PoaOptions, poa
from .polish import polish


def _check(count: Int, anchors: List[Int], strands: List[Int]) raises:
    """`consensus`'s anchors and strands for `count` sequences: none (the defaults), or one each and valid."""
    if len(anchors) != 0 and len(anchors) != count:
        raise Error(t"{len(anchors)} anchors for {count} sequences")
    if len(strands) != 0 and len(strands) != count:
        raise Error(t"{len(strands)} strands for {count} sequences")
    for r in range(len(anchors)):
        if anchors[r] != BOTH and anchors[r] != LEFT and anchors[r] != RIGHT:
            raise Error(t"read {r}: anchor {anchors[r]} is not {BOTH} (both), {LEFT} (left) or {RIGHT} (right)")
    for r in range(len(strands)):
        if strands[r] != 0 and strands[r] != 1 and strands[r] != STRAND_UNKNOWN:
            raise Error(t"read {r}: strand {strands[r]} is not 0 (forward), 1 (reverse) or {STRAND_UNKNOWN} (unknown)")


struct ConsensusOptions(Copyable, Movable):
    """`poa`'s options and the polishing rounds: `polish_rounds` scored by edit distance, then `likely_rounds` by
    likelihood (`plegma poa -P`), the defaults `plegma poa`'s."""

    var poa: PoaOptions
    var polish_rounds: Int
    var likely_rounds: Int

    def __init__(out self):
        self.poa = PoaOptions()
        self.polish_rounds = 50
        self.likely_rounds = 20


def consensus(
    var seqs: List[List[UInt8]],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> List[List[UInt8]]:
    """The consensus of `seqs` (A, C, G, T in upper case; any other letter matches none), each anchored as `anchors`
    says (`BOTH`, `LEFT` or `RIGHT`; empty: all `BOTH`) and on `strands`' strand (0, 1 or `STRAND_UNKNOWN`, for
    `homopolymers`; empty: all unknown), as `plegma poa -r 0` writes it (with `-H` when `homopolymers` is loaded).
    Raises on anchors or strands of another length than `seqs`, or values outside those."""
    _check(len(seqs), anchors, strands)
    if len(seqs) == 0:
        return [List[UInt8]()]
    if len(anchors) == 0:
        var all_spanning = List[Int](length=len(seqs), fill=BOTH)
        return consensus(seqs^, all_spanning, strands, options, homopolymers)
    var names = List[String](capacity=len(seqs))
    for k in range(len(seqs)):
        names.append(String(t"r{k}"))
    var made = poa(names, seqs, anchors, options.poa, Graph())
    var polishing = options.polish_rounds > 0 or options.likely_rounds > 0 or homopolymers.loaded()
    polishing = polishing and options.poa.mode == GLOBAL
    if len(made.consensus) > 0:
        # Reads anchored at one end: the assembly's own consensus, polished against every read but another copy's,
        # each aligned as it is anchored; two sides not joined, each as it is.
        if len(made.consensus) > 1 or not polishing:
            return made.consensus.copy()
        var members = List[Int]()
        var f = 0
        for r in range(len(seqs)):
            if f < len(made.foreign) and made.foreign[f] == r:
                f += 1
                continue
            members.append(r)
        var polished = polish(
            made.consensus[0],
            seqs,
            members,
            options.polish_rounds,
            likely_rounds=options.likely_rounds,
            anchors=anchors,
            homopolymers=homopolymers,
            strands=strands,
        )
        return [polished^]
    # Every read spanning: the heaviest bundle over every path, polished against every read.
    var letters = List[UInt8]()
    for node in made.graph.heaviest_path(Topology(made.graph).order, made.shown):
        letters.append(made.graph.bases[node])
    if not polishing:
        return [letters^]
    var members = List[Int]()
    for r in range(len(seqs)):
        members.append(r)
    var polished = polish(
        letters,
        seqs,
        members,
        options.polish_rounds,
        likely_rounds=options.likely_rounds,
        homopolymers=homopolymers,
        strands=strands,
    )
    return [polished^]
