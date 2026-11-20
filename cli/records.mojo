"""The command's input: one FASTA or FASTQ file (gzip is fine), read once, as names, sequences and each record's
local-assembly tags from the words of its header after the name (`anchor=both|left|right`, `strand=0|1`)."""

from sporos.fasta import parse_sequences, read_file

from plegma.homopolymer import STRAND_UNKNOWN
from plegma.kernel import BOTH, LEFT, RIGHT

comptime _NEWLINE = UInt8(ord("\n"))
comptime _RETURN = UInt8(ord("\r"))
comptime _SPACE = UInt8(ord(" "))
comptime _TAB = UInt8(ord("\t"))


def upper(seq: Span[UInt8, _]) -> List[UInt8]:
    """Upper case, as abPOA reads its input: a and A are one letter to the kernel and to the graph alike."""
    var out = List[UInt8](seq)
    for k in range(len(out)):
        if out[k] >= UInt8(ord("a")) and out[k] <= UInt8(ord("z")):
            out[k] -= 32
    return out^


struct Records(Movable, Sized):
    """A file's records: `names`, `seqs` (upper case), and per record its `anchors` (`BOTH` without a tag) and
    `strands` (`STRAND_UNKNOWN` without one)."""

    var names: List[String]
    var seqs: List[List[UInt8]]
    var anchors: List[Int]
    var strands: List[Int]

    def __init__(out self, path: String) raises:
        var bytes = read_file(path)
        var records = parse_sequences(Span(bytes))
        self.names = records.names.copy()
        self.seqs = List[List[UInt8]](capacity=len(records))
        for r in range(len(records)):
            self.seqs.append(upper(records.seq(r)))
        self.anchors = List[Int](length=len(records), fill=BOTH)
        self.strands = List[Int](length=len(records), fill=STRAND_UNKNOWN)
        var headers = _headers(Span(bytes))
        for r in range(min(len(headers), len(records))):
            for word in headers[r].split():
                if word == "anchor=left":
                    self.anchors[r] = LEFT
                elif word == "anchor=right":
                    self.anchors[r] = RIGHT
                elif word == "anchor=both":
                    self.anchors[r] = BOTH
                elif word == "strand=0":
                    self.strands[r] = 0
                elif word == "strand=1":
                    self.strands[r] = 1

    def __len__(self) -> Int:
        return len(self.seqs)


def _line_end(bytes: Span[UInt8, _], start: Int) -> Int:
    var i = start
    while i < len(bytes) and bytes[i] != _NEWLINE:
        i += 1
    return i


def _headers(bytes: Span[UInt8, _]) -> List[String]:
    """Each record's header line without its leading `>` or `@` and its line break, in file order. FASTA headers are
    the lines starting with `>`; FASTQ records are four lines each (blank lines between records skipped), as sporos
    reads them."""
    var out = List[String]()
    var fastq = False
    for i in range(len(bytes)):
        var c = bytes[i]
        if c != _NEWLINE and c != _RETURN and c != _SPACE and c != _TAB:
            fastq = c == UInt8(ord("@"))
            break
    var i = 0
    while i < len(bytes):
        var c = bytes[i]
        if fastq and (c == _NEWLINE or c == _RETURN):
            i += 1
            continue
        var end = _line_end(bytes, i)
        var trimmed = end
        if trimmed > i and bytes[trimmed - 1] == _RETURN:
            trimmed -= 1
        if fastq or c == UInt8(ord(">")):
            out.append(String(StringSlice(unsafe_from_utf8=bytes[i + 1 : max(trimmed, i + 1)])))
        if fastq:
            # The sequence, '+' and quality lines.
            for _ in range(3):
                end = _line_end(bytes, end + 1)
        i = end + 1
    return out^
