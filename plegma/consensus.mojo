"""The consensus of a set of reads as one call: `plegma poa`'s default output (`-r 0`) as a library function, for a
caller that assembles many small read sets (a local assembler's, between two anchors) and would rather not run the
command for each.

`consensus` takes the reads, each anchored as a local assembly's reads are (`BOTH`, `LEFT`, `RIGHT`) and on a strand
(`STRAND_UNKNOWN` when not known), and returns what `plegma poa -r 0` writes: one sequence, polished; or, when no read
spans both anchors and the two sides could not be joined, the left side's and the right side's, unpolished.
"""

from std.math import isfinite

from karan import Aligner

from .bands import PolishScratch
from .graph import Graph
from .homopolymer import STRAND_UNKNOWN, HomopolymerModel
from .kernel import BOTH, GLOBAL, INSIDE, LEFT, RIGHT
from .poa import Poa, PoaOptions, PoaScratch, poa_with_scratch
from .polish import polish_with_scratch


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
    likelihood (`plegma poa -P`), the defaults `plegma poa`'s; and for `clustered_consensus`, `noise_run` and
    `two_lengths_db`."""

    var poa: PoaOptions
    var polish_rounds: Int
    var likely_rounds: Int
    var noise_run: Int
    var two_lengths_db: Float64

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


struct ConsensusScratch(Movable):
    """Reusable scratch buffers for `consensus`, `polished`, `clustered_consensus`, and `phased_consensus` across
    calls."""

    var poa: PoaScratch
    var aligner: Aligner
    var polish: PolishScratch

    def __init__(out self):
        self.poa = PoaScratch()
        self.aligner = Aligner(warm=False)
        self.polish = PolishScratch()


def consensus(
    var seqs: List[List[UInt8]],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> List[List[UInt8]]:
    """The consensus of `seqs` (A, C, G, T in upper case; any other letter matches none), each anchored as `anchors`
    says (`BOTH`, `LEFT`, `RIGHT`, or `INSIDE`; empty: all `BOTH`) and on `strands`' strand."""
    var scratch = ConsensusScratch()
    return consensus(scratch, seqs^, anchors, strands, options, homopolymers)


def consensus(
    mut scratch: ConsensusScratch,
    var seqs: List[List[UInt8]],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> List[List[UInt8]]:
    """As `consensus`, reusing `scratch` across calls."""
    options.check()
    _check(len(seqs), anchors, strands)
    if len(seqs) == 0:
        return [List[UInt8]()]
    if len(anchors) == 0:
        var all_spanning = List[Int](length=len(seqs), fill=BOTH)
        return consensus(scratch, seqs^, all_spanning, strands, options, homopolymers)
    if INSIDE in anchors:
        var kept = List[List[UInt8]]()
        var kept_anchors = List[Int]()
        var kept_strands = List[Int]()
        for r in range(len(seqs)):
            if anchors[r] != INSIDE:
                kept.append(seqs[r].copy())
                kept_anchors.append(anchors[r])
                if len(strands) > 0:
                    kept_strands.append(strands[r])
        return consensus(scratch, kept^, kept_anchors, kept_strands, options, homopolymers)
    var names = List[String](capacity=len(seqs))
    for k in range(len(seqs)):
        names.append(String(t"r{k}"))
    var made = poa_with_scratch(names, seqs, anchors, options.poa, Graph(), scratch.poa, scratch.aligner)
    return polished(scratch, made, seqs, anchors, strands, options, homopolymers)


def polished(
    made: Poa,
    seqs: List[List[UInt8]],
    anchors: List[Int],
    strands: List[Int],
    options: ConsensusOptions,
    homopolymers: HomopolymerModel,
) raises -> List[List[UInt8]]:
    """The consensus of what `poa` made of `seqs`, polished when `options` has rounds or `homopolymers` is loaded."""
    var scratch = ConsensusScratch()
    return polished(scratch, made, seqs, anchors, strands, options, homopolymers)


def polished(
    mut scratch: ConsensusScratch,
    made: Poa,
    seqs: List[List[UInt8]],
    anchors: List[Int],
    strands: List[Int],
    options: ConsensusOptions,
    homopolymers: HomopolymerModel,
) raises -> List[List[UInt8]]:
    """As `polished`, reusing `scratch` across calls."""
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
        var out = polish_with_scratch(
            scratch.aligner,
            scratch.polish,
            made.consensus[0],
            seqs,
            members,
            options.polish_rounds,
            likely_rounds=options.likely_rounds,
            anchors=made.anchors if len(made.anchors) == len(seqs) else anchors,
            homopolymers=homopolymers,
            strands=strand_of,
            threads=options.poa.threads,
        )
        return [out^]
    var letters = List[UInt8]()
    for node in made.graph.heaviest_path(scratch.poa.order_of(made.graph), made.shown):
        letters.append(made.graph.bases[node])
    if not polishing or len(seqs) == 0:
        return [letters^]
    var members = List[Int](capacity=len(seqs))
    for r in range(len(seqs)):
        members.append(r)
    var out = polish_with_scratch(
        scratch.aligner,
        scratch.polish,
        letters,
        seqs,
        members,
        options.polish_rounds,
        likely_rounds=options.likely_rounds,
        homopolymers=homopolymers,
        strands=strand_of,
        threads=options.poa.threads,
    )
    return [out^]
