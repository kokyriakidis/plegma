"""The two haplotypes' consensus sequences of reads already phased (`phased_consensus`): each made from its labelled
reads, then again with the unlabelled reads on the haplotype whose letters they carry where the two differ
(`_sides`, at the places `_differences` finds), and each homopolymer run they differ by decided on both haplotypes'
reads together (`runs`)."""

from karan import Aligner, Alignment, Costs, Mode, Ties

from .alphabet import encode_base
from .consensus import ConsensusOptions
from .homopolymer import HomopolymerModel, STRAND_UNKNOWN, reads_say_one
from .kernel import BOTH, INSIDE, LEFT
from .placement import Clustered, alignment_costs, _consensus_of, _placed
from .poa import poa
from .trace import PHASED, Trace
from .runs import _best_length, _best_shared_length, _flanks, _read_lengths, _run_pairs, _set_runs


comptime TO_BOTH = 2
"""`_sides`: a read for both haplotypes, covering no place where the two differ."""


comptime MIN_CARRIERS = 1
"""`_sides`: a place where the two haplotypes' consensus sequences differ decides reads only where each spelling is
carried by at least this many reads of the window, labelled or not, and by `MIN_CARRIER_SHARE` of those spelling it
either way. The share does the work: a spelling few of many reads carry is a read's error made the consensus's, not a
second haplotype. HG002's HLA-DPA1 thinned to 15x had haplotype 2's left half built from 1 or 2 labelled reads, whose
36 errors read as 36 places apart from haplotype 1, and every unlabelled read, carrying the truth there, went to
haplotype 1 alone. Counting instead the labelled reads behind each consensus (3 of each) set aside the true
differences low depth leaves under 3 labelled reads, and the unlabelled reads there, given to both haplotypes, made
heterozygotes homozygous (CYP4A11 among them, unflagged); a floor of 2 carriers did the same at 10x and below (NAT2,
MTHFR, HTR2A, GRIK1, XPC), a second allele there being one read in four."""


comptime MIN_CARRIER_SHARE = 0.25
"""`_sides`: the least share of the reads spelling a place either way that each spelling must have to decide reads
there, a second haplotype being about half. On HG002 thinned to 30x down to 5x (two read sets) this kept every
HLA-DPA1 consensus exact at 15x and 10x and left no wrong call unflagged in panStar; it still set aside a second allele
one read in five carried at 10x and below (NAT2, MTHFR thinned), as the reads alone cannot tell it from an error
there: at 0.15 those were kept, but so were DPA1's errors where one read in six carried them (33 edits at 15x)."""


comptime PHASED_CARRIERS = 3
"""`_sides`: a place also decides reads where pgphase's labels back it: each haplotype's spelling carried by at least
this many of its own labelled reads, and by `PHASED_AGREEMENT` of its labelled reads spelling it either way. The
share of all reads cannot tell a minority haplotype from an error (NA18941's CYP1A1: 8 reads of one haplotype, 3 of
the other, 4 of another phase set carrying the first's letter; at 3 of 15 the place was set aside, the 4 went to
both and outvoted the second's own 3), but an error does not follow the labels: one read's error made a thin
consensus's (HG002's HLA-DPA1 at 15x, a haplotype built from 1 or 2 reads) has under 3 carriers of its own."""


comptime PHASED_AGREEMENT = 0.75


comptime SITE_FLANK = 8
"""Letters of the longer haplotype either side of a place where the two differ, compared with a read's letters there
along with it: two alignments may put the same gap a letter or two apart, which compared alone would read as the
other haplotype's spelling; with the flank both are one string. Places within twice this of each other are one."""


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


def _ops(cigar: String) -> List[Tuple[Int, UInt8]]:
    """karan's run-length CIGAR as (length, letter) pairs."""
    var out = List[Tuple[Int, UInt8]]()
    var number = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            number = number * 10 + Int(byte - UInt8(ord("0")))
        else:
            out.append((number, byte))
            number = 0
    return out^


def _query_at(alignment: Alignment) -> List[Int]:
    """For each reference letter `reference_start ..= reference_end` of `alignment`, the query letter it stands
    before (a gap's letters count before the reference letter after them): the query's letters over reference
    `s ..< e` are `query_at[s - reference_start] ..< query_at[e - reference_start]`."""
    var out = List[Int](capacity=alignment.reference_end - alignment.reference_start + 1)
    var q = alignment.query_start
    var ops = _ops(alignment.cigar)
    for ref op in ops:
        var length = op[0]
        var letter = op[1]
        if letter == UInt8(ord("I")):
            q += length
            continue
        for _ in range(length):
            out.append(q)
            if letter != UInt8(ord("D")):
                q += 1
    out.append(q)
    return out^


def _differences(mut aligner: Aligner, longer: List[UInt8], shorter: List[UInt8], costs: Costs) raises -> List[_Site]:
    """Where `shorter`, placed inside `longer` (`Mode.INFIX`: a consensus cut short where its reads stop is not
    thereby another haplotype), spells it otherwise, each place widened by `SITE_FLANK` and those that then touch
    joined."""
    var out = List[_Site]()
    if len(longer) == 0 or len(shorter) == 0 or longer == shorter:
        return out^
    var placed = aligner.align(Span(longer), Span(shorter), costs, Mode.INFIX, ties=Ties.LEFT)
    var at = _query_at(placed)
    var spans = List[Tuple[Int, Int]]()
    var r = placed.reference_start
    var ops = _ops(placed.cigar)
    for ref op in ops:
        var length = op[0]
        var letter = op[1]
        if letter == UInt8(ord("=")) or letter == UInt8(ord("M")):
            r += length
            continue
        var consumed = 0 if letter == UInt8(ord("I")) else length
        var start = max(placed.reference_start, r - SITE_FLANK)
        var end = min(placed.reference_end, r + consumed + SITE_FLANK)
        if len(spans) > 0 and start <= spans[len(spans) - 1][1]:
            spans[len(spans) - 1] = (spans[len(spans) - 1][0], max(end, spans[len(spans) - 1][1]))
        else:
            spans.append((start, end))
        r += consumed
    for ref span in spans:
        var long_letters = List[UInt8]()
        for i in range(span[0], span[1]):
            long_letters.append(longer[i])
        var short_start = at[span[0] - placed.reference_start]
        var short_end = at[span[1] - placed.reference_start]
        var short_letters = List[UInt8]()
        for i in range(short_start, short_end):
            short_letters.append(shorter[i])
        out.append(_Site(span[0], span[1], short_start, short_end, long_letters^, short_letters^))
    return out^


def _haplotype_at(alignment: Alignment) -> List[Int]:
    """`_query_at` the other way round: for a haplotype placed inside a read (`alignment`'s query the haplotype, its
    reference the read; `Mode.INFIX`, so every haplotype letter is placed), for each haplotype letter `0 ..= len` the
    read letter it stands before, a gap's letters counting before the haplotype letter after them."""
    var out = List[Int](capacity=alignment.query_end - alignment.query_start + 1)
    var r = alignment.reference_start
    var ops = _ops(alignment.cigar)
    for ref op in ops:
        var length = op[0]
        var letter = op[1]
        if letter == UInt8(ord("D")):
            r += length
            continue
        for _ in range(length):
            out.append(r)
            if letter != UInt8(ord("I")):
                r += 1
    out.append(r)
    return out^


def _spelling(seq: List[UInt8], start: Int, end: Int, at: List[Int], lo: Int, hi: Int, site: _Site) -> Int:
    """Which spelling of `site` a read carries, placed over a haplotype's letters `start ..< end` (`at`: the read
    letter each haplotype letter from `start` stands before), the site lying on that haplotype's `lo ..< hi`: 0 the
    longer's spelling, 1 the shorter's, -1 neither, or a place the read does not span wholly short of its
    placement's ends, where a read running past the haplotype's end piles its overhang."""
    if lo <= start or hi >= end:
        return -1
    var letters = List[UInt8]()
    for q in range(at[lo - start], at[hi - start]):
        letters.append(seq[q])
    if letters == site.longer:
        return 0
    if letters == site.shorter:
        return 1
    return -1


def _sides(
    mut aligner: Aligner,
    made: List[List[UInt8]],
    seqs: List[List[UInt8]],
    builders: List[Int],
    reads: List[Int],
    costs: Costs,
    labels: List[Int] = List[Int](),
) raises -> List[Int]:
    """For each of `reads`, the haplotype of `made` whose letters it carries where the two differ: 0 or 1, `TO_BOTH`
    when it covers no such place (a stretch where the two read alike, or one neither reaches), -1 when it carries as
    many of one's as of the other's. Every read of `builders` (labelled or not) is placed on a haplotype
    (`Mode.INFIX`) and its letters over each place it spans wholly compared with the two spellings there -- a read
    labelled with the shorter's haplotype on the shorter, at each place's span there, every other read on the
    longer; a place
    decides only where each spelling is carried by `MIN_CARRIERS` of them and `MIN_CARRIER_SHARE` of those spelling
    it either way, or where `labels` (a haplotype per read of `seqs`, -1 none) back it (`PHASED_CARRIERS`)."""
    var out = List[Int](length=len(reads), fill=TO_BOTH)
    var long_index = 0 if len(made[0]) >= len(made[1]) else 1
    var short_index = 1 - long_index
    ref longer = made[long_index]
    var sites = _differences(aligner, longer, made[short_index], costs)
    if len(sites) == 0:
        return out^
    # Each read's spelling at each place: row per read of `builders`.
    var spelt = List[List[Int]]()
    var row_of = List[Int](length=len(seqs), fill=-1)
    for b in range(len(builders)):
        ref seq = seqs[builders[b]]
        # A read labelled with the shorter's haplotype is placed on the shorter, the places read at their spans there
        # (`_Site.short_start`): it carries that haplotype's letters, and placed on the longer it would differ from it
        # wherever the two do. NA19007's HLA-DRB1 (*09:01 and *08:03, hundreds of edits apart through the introns):
        # 77 reads placed on the other haplotype at half a second each, 39 of the window's 47 s, for 2 unlabelled
        # reads. Every other read is placed on the longer.
        var h = labels[builders[b]] if len(labels) > 0 else -1
        var on_shorter = h == short_index
        ref haplotype = made[short_index if on_shorter else long_index]
        # A read is placed inside the haplotype, unless it is the longer of the two: then the haplotype inside it.
        # A haplotype built from a few labelled reads can be a fraction of the window (HG02106's HLA-DPB1: 2,716 and
        # 3,732 letters from 7 of 71 reads), and a 15.7 kb read forced inside 3.7 kb is a 12 kb gap, which an
        # alignment's cost grows with: 1.8 s a read, 57 of the window's 59 s. The other way round the two differ by
        # what the read's errors and the haplotypes' differences make, and every haplotype letter has its read
        # letter.
        var start: Int
        var end: Int
        var at: List[Int]
        if len(seq) > len(haplotype):
            var inside = aligner.align(Span(seq), Span(haplotype), costs, Mode.INFIX, ties=Ties.LEFT)
            at = _haplotype_at(inside)
            start = 0
            end = len(haplotype)
        else:
            var placed = aligner.align(Span(haplotype), Span(seq), costs, Mode.INFIX, ties=Ties.LEFT)
            at = _query_at(placed)
            start = placed.reference_start
            end = placed.reference_end
        var row = List[Int](length=len(sites), fill=-1)
        for k in range(len(sites)):
            ref site = sites[k]
            if on_shorter:
                row[k] = _spelling(seq, start, end, at, site.short_start, site.short_end, site)
            else:
                row[k] = _spelling(seq, start, end, at, site.start, site.end, site)
        spelt.append(row^)
        row_of[builders[b]] = b
    var counted = List[Bool](length=len(sites), fill=False)
    for k in range(len(sites)):
        var carriers = [0, 0]
        # own[2 * h + c]: reads labelled h spelling the place c (0 the longer haplotype's spelling).
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
            # Labels back it: the longer haplotype's reads spell it 0, the other's 1.
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
                # Spelling 0 is the longer haplotype's.
                votes[long_index if row[k] == 0 else 1 - long_index] += 1
        if votes[0] > votes[1]:
            out[i] = 0
        elif votes[1] > votes[0]:
            out[i] = 1
        elif votes[0] > 0:
            out[i] = -1
    return out^


def _any_builds(anchor_of: List[Int]) -> Bool:
    """Whether any read is anchored (`BOTH`, `LEFT` or `RIGHT`), not inside."""
    for a in anchor_of:
        if a != INSIDE:
            return True
    return False


def phased_consensus(
    var seqs: List[List[UInt8]],
    haplotypes: List[Int],
    anchors: List[Int] = List[Int](),
    strands: List[Int] = List[Int](),
    options: ConsensusOptions = ConsensusOptions(),
    homopolymers: HomopolymerModel = HomopolymerModel(),
) raises -> Clustered:
    """The two haplotypes' consensus sequences of reads already phased (`haplotypes`: 0 or 1 per read, -1 for a read
    of neither), each made as `consensus` makes one from its reads, and every read's alignment to its haplotype's
    (see `clustered_consensus`). A read of neither builds as well: made first from the labelled reads, the two are
    made again with each such read on the haplotype whose letters it carries where they differ, on both where it
    covers no such place, on neither where it carries as many of each (`_sides`); it is then placed where it fits
    best. Where the two differ by the length of a run of one letter and a homopolymer model is loaded, each such run
    is decided on both haplotypes' reads together, whatever else the two differ by: their own lengths for the two only
    when that explains the reads `options.two_lengths_db` better than one length for both, else that one length in
    both (a run is far more often the same on the two haplotypes than not). Two clusters always: a haplotype no read
    is labelled with has the other's consensus, and with no labelled read at all both are the unlabelled reads'. With
    every read inside, all of them build as if spanning (as `clustered_consensus` has it)."""
    var n = len(seqs)
    if len(haplotypes) != n:
        raise Error(t"{len(haplotypes)} haplotypes for {n} sequences")
    if len(anchors) != 0 and len(anchors) != n:
        raise Error(t"{len(anchors)} anchors for {n} sequences")
    if len(strands) != 0 and len(strands) != n:
        raise Error(t"{len(strands)} strands for {n} sequences")
    var anchor_of = anchors.copy() if len(anchors) > 0 else List[Int](length=n, fill=BOTH)
    var strand_of = strands.copy() if len(strands) > 0 else List[Int](length=n, fill=STRAND_UNKNOWN)
    # The reads that build: all but those inside, or with every read inside all of them as if spanning.
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
        # No labelled read builds: both haplotypes the unlabelled reads' one consensus.
        groups[0] = unlabelled.copy()
        groups[1] = unlabelled.copy()
        unlabelled.clear()
    var made = List[List[UInt8]]()
    for k in range(2):
        made.append(
            _consensus_of(seqs, groups[k], build_anchor, strand_of, options, homopolymers) if len(groups[k])
            > 0 else (List[UInt8]())
        )
    if len(groups[1]) == 0:
        made[1] = made[0].copy()
    if len(groups[0]) == 0:
        made[0] = made[1].copy()
    var aligner = Aligner(warm=False)
    # Reads of neither haplotype build too: each joins the haplotype whose letters it carries where the two differ,
    # or both where it covers no such place, and the two are made again. Left out of building, they left a haplotype
    # pieced together from its few labelled reads wherever pgphase's phase set ends inside the window: HG002's
    # HLA-DPA1 thinned to 15x, its 14 reads of the window's left end all unlabelled, had haplotype 2's left half from
    # its one spanning read, 36 edits from the truth; at 10x, with no labelled read reaching that end, both
    # haplotypes stopped short, at 3,915 and 8,058 letters of 18,180. With them: 0 and 1 edits, then 0 and 0.
    if len(unlabelled) > 0:
        var costs = alignment_costs(options.poa.scoring)
        var builders = groups[0].copy()
        builders.extend(groups[1].copy())
        builders.extend(unlabelled.copy())
        var sides = _sides(aligner, made, seqs, builders, unlabelled, costs, haplotypes)
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
        var grew = [False, False]
        for i in range(len(unlabelled)):
            for k in range(2):
                if sides[i] == k or sides[i] == TO_BOTH:
                    groups[k].append(unlabelled[i])
                    grew[k] = True
        for k in range(2):
            if grew[k] and labelled[k]:
                made[k] = _consensus_of(seqs, groups[k], build_anchor, strand_of, options, homopolymers)
        # A haplotype no read is labelled with stays the other's: unlabelled reads alone (each one there because it
        # covers no place the two differ, which they do not) would make one noisy read a second haplotype.
        for k in range(2):
            if not labelled[k] and labelled[1 - k]:
                made[k] = made[1 - k].copy()
    # Each run the two differ by is decided on both haplotypes' reads together, whatever else the two differ by.
    var runs0 = List[Int]()
    var runs1 = List[Int]()
    if homopolymers.loaded() and made[0] != made[1]:
        _ = _run_pairs(aligner, made[0], made[1], max(options.noise_run, 2), runs0, runs1)
    if len(runs0) > 0:
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
            # Each haplotype's reads scored with its own consensus's flanks: the two may differ beside the run.
            var flanks0 = _flanks(made[0], s, e)
            var flanks1 = _flanks(made[1], s1, e1)
            # From 4 under the shorter of the two haplotypes' lengths to 4 over the longer: they may be far apart.
            var low = max(2, min(e - s, e1 - s1) - 4)
            var high = max(e - s, e1 - s1) + 4
            if not homopolymers.covers(base, high):
                # A run longer than the model's rows: it cannot say, so each haplotype keeps its own consensus's
                # length (searched within the rows alone, both came out at the longest row: T59 as T50).
                first.append(e - s)
                second.append(e1 - s1)
                continue
            # Each haplotype's reads measured on its own consensus: the two may differ beside the run.
            _read_lengths(aligner, made[0], s, e, seqs, groups[0], anchor_of, lengths0)
            _read_lengths(aligner, made[1], s1, e1, seqs, groups[1], anchor_of, lengths1)
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
                # Too few reads to say: as the two consensus sequences have it.
                first.append(e - s)
                second.append(e1 - s1)
                continue
            if score0 + score1 - shared_score > options.two_lengths_db:
                first.append(own0)
                second.append(own1)
            else:
                first.append(shared)
                second.append(shared)
            # A haplotype whose reads report the letter once there has no run (`reads_say_one`): the search above
            # starts at 2 (chr20:8,254,476, 42 reads T against the other haplotype's 26 TTT: written TT, the truth T).
            if reads_say_one(lengths0):
                first[k] = 1
            if reads_say_one(lengths1):
                second[k] = 1
        var a = made[0].copy()
        var b = made[1].copy()
        made[0] = _set_runs(a, runs0, first)
        made[1] = _set_runs(b, runs1, second)
    return _placed(aligner, made^, seqs, cluster_of, anchor_of, alignment_costs(options.poa.scoring))
