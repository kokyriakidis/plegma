"""The consensus of a set of reads as one call: `plegma poa`'s default output (`-r 0`) as a library function, for a
caller that assembles many small read sets (a local assembler's, between two anchors) and would rather not run the
command for each.

`consensus` takes the reads, each anchored as a local assembly's reads are (`BOTH`, `LEFT`, `RIGHT`) and on a strand
(`STRAND_UNKNOWN` when not known), and returns what `plegma poa -r 0` writes: one sequence, polished; or, when no read
spans both anchors and the two sides could not be joined, the left side's and the right side's, unpolished.
"""

from std.math import isfinite

from .graph import Graph
from .homopolymer import STRAND_UNKNOWN, HomopolymerModel
from .kernel import BOTH, GLOBAL, INSIDE, LEFT, RIGHT, Topology
from .poa import Poa, PoaOptions, poa
from .polish import polish


def _check(count: Int, anchors: List[Int], strands: List[Int]) raises:
    """`consensus`'s anchors and strands for `count` sequences: none (the defaults), or one each and valid."""
    if len(anchors) != 0 and len(anchors) != count:
        raise Error(t"{len(anchors)} anchors for {count} sequences")
    if len(strands) != 0 and len(strands) != count:
        raise Error(t"{len(strands)} strands for {count} sequences")
    for r in range(len(anchors)):
        var a = anchors[r]
        if a != BOTH and a != LEFT and a != RIGHT and a != INSIDE:
            raise Error(
                t"read {r}: anchor {a} is not {BOTH} (both), {LEFT} (left), {RIGHT} (right) or {INSIDE} (inside)"
            )
    for r in range(len(strands)):
        if strands[r] != 0 and strands[r] != 1 and strands[r] != STRAND_UNKNOWN:
            raise Error(t"read {r}: strand {strands[r]} is not 0 (forward), 1 (reverse) or {STRAND_UNKNOWN} (unknown)")


struct ConsensusOptions(Copyable, Movable):
    """`poa`'s options and the polishing rounds: `polish_rounds` scored by edit distance, then `likely_rounds` by
    likelihood (`plegma poa -P`), the defaults `plegma poa`'s; and for `clustered_consensus`, `noise_run`: clusters
    whose consensus sequences differ only by single letters in runs of one letter this long or longer are one (0:
    never)."""

    var poa: PoaOptions
    var polish_rounds: Int
    var likely_rounds: Int
    var noise_run: Int
    var two_lengths_db: Float64
    """`clustered_consensus` and `phased_consensus`: how many dB better two run lengths must explain the reads than
    one for the two clusters to keep different lengths, with a homopolymer model. On HG002 HiFi chr20 (the reads'
    truth haplotypes) 15 dB kept 23 of 57 unphased splits the truth has heterozygous and 4 of 421 it has homozygous."""

    def __init__(out self):
        self.poa = PoaOptions()
        self.polish_rounds = 50
        self.likely_rounds = 20
        self.noise_run = 3
        self.two_lengths_db = 15.0

    def check(self) raises:
        """Raises on options the consensus calls cannot run with (`PoaOptions.check` too), naming the value."""
        self.poa.check()
        if self.polish_rounds < 0 or self.likely_rounds < 0:
            raise Error(t"polishing rounds {self.polish_rounds},{self.likely_rounds}: each must be 0 or more")
        if self.noise_run < 0:
            raise Error(t"noise_run {self.noise_run}: must be 0 or more")
        if not isfinite(self.two_lengths_db):
            raise Error(t"two_lengths_db {self.two_lengths_db}: must be a number")


def consensus(
    var seqs: List[List[UInt8]],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> List[List[UInt8]]:
    """The consensus of `seqs` (A, C, G, T in upper case; any other letter matches none), each anchored as `anchors`
    says (`BOTH`, `LEFT` or `RIGHT`, or `INSIDE`, a read covering neither anchor, which makes no consensus; empty:
    all `BOTH`) and on `strands`' strand (0, 1 or `STRAND_UNKNOWN`, for
    `homopolymers`; empty: all unknown), as `plegma poa -r 0` writes it (with `-H` when `homopolymers` is loaded).
    Raises on anchors or strands of another length than `seqs`, or values outside those, or on `options` it cannot
    run with (`ConsensusOptions.check`)."""
    options.check()
    _check(len(seqs), anchors, strands)
    if len(seqs) == 0:
        return [List[UInt8]()]
    if len(anchors) == 0:
        var all_spanning = List[Int](length=len(seqs), fill=BOTH)
        return consensus(seqs^, all_spanning, strands, options, homopolymers)
    if INSIDE in anchors:
        # Reads covering neither anchor place no letter at either end: the others make the consensus.
        var kept = List[List[UInt8]]()
        var kept_anchors = List[Int]()
        var kept_strands = List[Int]()
        for r in range(len(seqs)):
            if anchors[r] != INSIDE:
                kept.append(seqs[r].copy())
                kept_anchors.append(anchors[r])
                if len(strands) > 0:
                    kept_strands.append(strands[r])
        return consensus(kept^, kept_anchors, kept_strands, options, homopolymers)
    var names = List[String](capacity=len(seqs))
    for k in range(len(seqs)):
        names.append(String(t"r{k}"))
    var made = poa(names, seqs, anchors, options.poa, Graph())
    return polished(made, seqs, anchors, strands, options, homopolymers)


def polished(
    made: Poa,
    seqs: List[List[UInt8]],
    anchors: List[Int],
    strands: List[Int],
    options: ConsensusOptions,
    homopolymers: HomopolymerModel,
) raises -> List[List[UInt8]]:
    """The consensus of what `poa` made of `seqs` (each anchored as `anchors` says, on `strands`' strand), polished
    (when `options` has rounds or `homopolymers` is loaded, and the alignment is global): with reads anchored at one
    end, the assembly's own consensus against every read but another copy's (`Poa.foreign`), each aligned as it is
    anchored, or two sides not joined, each as it is; with every read spanning, the heaviest bundle over every path
    (a starting MSA's rows too) against every read. What `consensus` and `plegma poa -r 0` give. A read `poa`
    reverse-complemented (`Poa.reversed`) is on the other strand of `strands` to the consensus."""
    var strand_of = strands.copy()
    for r in made.reversed:
        if r < len(strand_of) and strand_of[r] != STRAND_UNKNOWN:
            strand_of[r] = 1 - strand_of[r]
    var polishing = options.polish_rounds > 0 or options.likely_rounds > 0 or homopolymers.loaded()
    polishing = polishing and options.poa.mode == GLOBAL
    if len(made.consensus) > 0:
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
            anchors=made.anchors if len(made.anchors) == len(seqs) else anchors,
            homopolymers=homopolymers,
            strands=strand_of,
        )
        return [polished^]
    var letters = List[UInt8]()
    for node in made.graph.heaviest_path(Topology(made.graph).order, made.shown):
        letters.append(made.graph.bases[node])
    if not polishing or len(seqs) == 0:
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
        strands=strand_of,
    )
    return [polished^]
