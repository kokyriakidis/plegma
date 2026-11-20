"""`plegma`'s command line: the usage text and the options of its two commands, parsed and checked."""

from plegma.consensus import ConsensusOptions
from plegma.homopolymer import HomopolymerModel
from plegma.kernel import GLOBAL
from plegma.poa import GUIDE_ORDER, INPUT_ORDER, LONGEST_FIRST


comptime USAGE = """Usage: plegma <command> [options] <in.fa/fq> > cons.fa/msa.fa/plegma.gfa

Commands (plegma --version prints the version):
  poa           partial order alignment of the sequences, as abPOA runs one
  ref           the MSA anchored on the first sequence (the reference): every other aligned to it alone at
                minimap2's asm5 scores (as costs: mismatch 40, a gap of k letters min(78 + 7k, 162 + 3k)), gaps
                left-normalized, its letters in the reference's columns and the runs inserted at one place aligned to
                each other on the longest; reads of one sequence past asm5's divergence fall back to pgrealign's scores

poa options:
    -m INT      alignment mode [0]: 0 global, 1 local, 2 extension
    -M INT      match score [2]
    -X INT      mismatch penalty [4]
    -O INT[,INT] gap opening penalties O1,O2 [4,24]
    -E INT[,INT] gap extension penalties E1,E2 [2,1]
                a g-long gap costs min{O1+g*E1, O2+g*E2} (convex); O2 0: O1+g*E1 (affine); O1 0: g*E1 (linear)
    -b INT      adaptive band: extra band b [10]; < 0 disables the band
    -f FLOAT    adaptive band: extra band f [0.01]: b+f*L letters either side, L the sequence's length
    -B INT      adaptive band: at most INT letters either side [64]; < 0 no cap (abPOA's band)
    -p          progressive: abPOA's guide order, most similar minimizer multisets first [the default]
    -L          the sequences longest first instead
    -n          the sequences in input order instead
    -S          seeded (global mode): each sequence pinned to the consensus path as minimap2 pins an alignment
                (asm5's seeds and chain; pins only on 19-mers unique in both), the DP only between the pins
                [the default]
    -N          not seeded: each sequence aligned to the whole graph
    -s          ambiguous strand: a sequence scoring poorly is tried reverse-complemented, the better kept
    -i FILE     start from the graph of an existing MSA (FASTA, rows with '-')
    -c INT      a one-directional assembly (reads fixed at one anchor, none spanning) stops where fewer reads
                go on [5]
                A read's anchor is a word of its FASTA header: anchor=both (the default: it spans the local assembly's
                two anchors), anchor=left or anchor=right (fixed at that anchor only, ending inside the gap). With
                reads of both sides and none spanning, each side is assembled and the two joined where they overlap.
    -H FILE     a homopolymer model (shasta2's csv, conf/HomopolymerModel-*.csv): after polishing, every run of two
                letters or more takes its most likely length given the reads' lengths there (shasta2's msa1 repair),
                each read on the strand its header says: strand=0 (as the consensus runs) or strand=1 (reverse-
                complemented to it); both weighed alike without

ref options:
    -F FILE     the reference's flanks (FASTA: left, then right): each alignment to the reference pinned as minimap2
                pins it (-x asm5 as dipcall runs it; its seeds and chain over flank + sequence + flank), exact between
                the pins, so differences are written as dipcall's truth sets write them
    -t INT      threads for the alignments to the reference [1]

Output options (both):
    -r INT      output [0]: 0 consensus FASTA, 1 MSA FASTA, 2 both, 3 GFA, 4 GFA with the consensus path
    -a INT      consensus [0]: 0 heaviest path, 1 most frequent letter per column
    -d INT      at most INT consensus sequences [1]; at 2 or more the sequences are clustered on the MSA's
                heterozygous columns and each cluster gets its own (-d and -a: reads spanning both anchors only)
    -q FLOAT    smallest share of the sequences a cluster and a column's second allele must hold [0.25]
    -P INT[,INT] rounds of consensus polishing, at most [50,20] (global alignment): each sequence aligned to the
                consensus, and every one-letter or whole-repeat-unit edit kept that lowers their summed edit distance
                to it (the first), then that raises their likelihood under a pair HMM (the second); 0,0 off
    -l          the input is a list of sequence files, one per line, each aligned on its own
    -o FILE     output to FILE [stdout]
    -v          report timings on stderr
    -h          this help
"""


struct Options(Copyable, Movable):
    """A command line, parsed: the command (`reference`: `ref`, else `poa`), the library's options for it
    (`consensus`: `poa`'s and the polishing rounds), and the command's own."""

    var reference: Bool
    var consensus: ConsensusOptions
    var start_msa: String
    var flanks: String
    var threads: Int
    var result: Int
    var consensus_rule: Int
    var max_consensus: Int
    var min_freq: Float64
    var homopolymers: HomopolymerModel
    var verbose: Bool
    var output: String
    """`-o`: the output file, empty for stdout."""
    var listed: Bool
    """`-l`: `input` lists the files to align, one per line."""
    var input: String

    def __init__(out self, reference: Bool):
        self.reference = reference
        # The guide order and the pins by default (`PoaOptions`): on the read sets the same consensus as the input
        # order unpinned, on chr20's allele sets better MSAs, on its largest VNTRs several times faster.
        self.consensus = ConsensusOptions()
        self.start_msa = String()
        self.flanks = String()
        self.threads = 1
        self.result = 0
        self.consensus_rule = 0
        self.max_consensus = 1
        self.min_freq = 0.25
        self.homopolymers = HomopolymerModel()
        self.verbose = False
        self.output = String()
        self.listed = False
        self.input = String()

    def polishing(self) -> Bool:
        """Whether a consensus is polished: some rounds (or a homopolymer model) and global alignment."""
        var rounds = self.consensus.polish_rounds > 0 or self.consensus.likely_rounds > 0
        return (rounds or self.homopolymers.loaded()) and self.consensus.poa.mode == GLOBAL


def _value(args: List[String], i: Int) raises -> String:
    if i + 1 >= len(args):
        raise Error(t"{args[i]} takes a value")
    return args[i + 1]


def _int(args: List[String], i: Int) raises -> Int:
    var text = _value(args, i)
    try:
        return Int(text)
    except:
        raise Error(t"{args[i]} takes an integer, not '{text}'")


def _float(args: List[String], i: Int) raises -> Float64:
    var text = _value(args, i)
    try:
        return Float64(text)
    except:
        raise Error(t"{args[i]} takes a number, not '{text}'")


def _pair(args: List[String], i: Int) raises -> Tuple[Int, Int]:
    """`INT` or `INT,INT`; the second -1 when absent."""
    var text = _value(args, i)
    try:
        var parts = text.split(",")
        if len(parts) == 1:
            return (Int(parts[0]), -1)
        if len(parts) == 2:
            return (Int(parts[0]), Int(parts[1]))
    except:
        pass
    raise Error(t"{args[i]} takes INT or INT,INT, not '{text}'")


def _check(options: Options) raises:
    """Values the commands cannot run with."""
    # The library's own checks (-m, -M -X -O -E, -f, -c, -P), its words naming the value.
    options.consensus.check()
    if options.result < 0 or options.result > 4:
        raise Error(t"-r {options.result}: the outputs are 0 to 4")
    if options.consensus_rule != 0 and options.consensus_rule != 1:
        raise Error(t"-a {options.consensus_rule}: the consensus rules are 0 and 1")
    if options.max_consensus < 1:
        raise Error(t"-d {options.max_consensus}: at least one consensus sequence")
    if options.min_freq < 0 or options.min_freq > 1:
        raise Error(t"-q {options.min_freq}: a share between 0 and 1")
    if options.threads < 1:
        raise Error(t"-t {options.threads}: at least one thread")


def parse(args: List[String]) raises -> Options:
    """`args` (the command first) as `Options`; raises on an unknown option, a missing or malformed value, or a value
    out of range. A missing input file is left to the caller (empty `input`)."""
    if args[0] != "poa" and args[0] != "ref":
        raise Error(t"unknown command '{args[0]}': poa or ref (plegma -h for help)")
    var reference = args[0] == "ref"
    var options = Options(reference)
    ref poa = options.consensus.poa
    var i = 1
    while i < len(args):
        var a = args[i]
        var takes_value = True
        # poa's own.
        if a == "-m" and not reference:
            poa.mode = _int(args, i)
        elif a == "-M" and not reference:
            poa.scoring.match_score = Int32(_int(args, i))
        elif a == "-X" and not reference:
            poa.scoring.mismatch = Int32(_int(args, i))
        elif a == "-O" and not reference:
            var o = _pair(args, i)
            poa.scoring.open1 = Int32(o[0])
            if o[1] >= 0:
                poa.scoring.open2 = Int32(o[1])
        elif a == "-E" and not reference:
            var e = _pair(args, i)
            poa.scoring.extend1 = Int32(e[0])
            if e[1] >= 0:
                poa.scoring.extend2 = Int32(e[1])
        elif a == "-b" and not reference:
            poa.band_b = _int(args, i)
        elif a == "-f" and not reference:
            poa.band_f = _float(args, i)
        elif a == "-B" and not reference:
            poa.band_max = _int(args, i)
        elif a == "-c" and not reference:
            poa.min_cover = _int(args, i)
        elif a == "-i" and not reference:
            options.start_msa = _value(args, i)
        elif a == "-H" and not reference:
            options.homopolymers = HomopolymerModel.load(_value(args, i))
        # ref's own.
        elif a == "-F" and reference:
            options.flanks = _value(args, i)
        elif a == "-t" and reference:
            options.threads = _int(args, i)
        # Both.
        elif a == "-o":
            options.output = _value(args, i)
        elif a == "-r":
            options.result = _int(args, i)
        elif a == "-a":
            options.consensus_rule = _int(args, i)
        elif a == "-P":
            var both = _pair(args, i)
            options.consensus.polish_rounds = both[0]
            if both[1] >= 0:
                options.consensus.likely_rounds = both[1]
        elif a == "-d":
            options.max_consensus = _int(args, i)
        elif a == "-q":
            options.min_freq = _float(args, i)
        else:
            takes_value = False
            # Flags: the last of -p, -L and -n wins.
            if a == "-p" and not reference:
                poa.order = GUIDE_ORDER
            elif a == "-L" and not reference:
                poa.order = LONGEST_FIRST
            elif a == "-n" and not reference:
                poa.order = INPUT_ORDER
            elif a == "-s" and not reference:
                poa.strand = True
            elif a == "-S" and not reference:
                poa.seeded = True
            elif a == "-N" and not reference:
                poa.seeded = False
            elif a == "-l":
                options.listed = True
            elif a == "-v":
                options.verbose = True
            elif a.startswith("-") and a.byte_length() > 1:
                raise Error(t"{args[0]}: unknown option {a} (plegma -h for help)")
            elif options.input.byte_length() > 0:
                raise Error(t"{args[0]}: one input file, not '{options.input}' and '{a}'")
            else:
                options.input = a
        i += 2 if takes_value else 1
    _check(options)
    return options^
