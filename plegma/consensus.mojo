"""The consensus of a set of reads as one call: `plegma poa`'s default output (`-r 0`) as a library function, for a
caller that assembles many small read sets (a local assembler's, between two anchors) and would rather not run the
command for each.

`consensus` takes the reads, each anchored as a local assembly's reads are (`BOTH`, `LEFT`, `RIGHT`) and on a strand
(`STRAND_UNKNOWN` when not known), and returns what `plegma poa -r 0` writes: one sequence, polished; or, when no read
spans both anchors and the two sides could not be joined, the left side's and the right side's, unpolished.
"""

from .graph import Graph
from .homopolymer import HomopolymerModel
from .kernel import GLOBAL, Topology
from .poa import PoaOptions, poa
from .polish import polish


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
    mut seqs: List[List[UInt8]],
    anchors: List[Int],
    strands: List[Int],
    options: ConsensusOptions,
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> List[List[UInt8]]:
    """The consensus of `seqs` (A, C, G, T in upper case), each anchored as `anchors` says and on `strands`' strand
    (for `homopolymers`; empty: all unknown), as `plegma poa -r 0` writes it (with `-H` when `homopolymers` is
    loaded). `seqs` may change: with `options.poa.strand` a spanning read is reverse-complemented when that aligns
    better."""
    if len(seqs) == 0:
        return [List[UInt8]()]
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
