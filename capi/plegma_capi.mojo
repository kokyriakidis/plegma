"""The C API of plegma (`capi/plegma.h`), built as `build/libplegma.so` (`pixi run build-capi`): the consensus of a read set
(`plegma.consensus`) for a C or C++ caller that assembles many, a homopolymer model loaded once and shared.

Handles are heap boxes the caller frees (`plegma_model_free`, `plegma_result_free`); nothing raises across the
boundary: a failure comes back as a result holding its message (`plegma_result_error`), or a null model with the
message written into the caller's buffer. A model is read-only once made, so threads may share one; every call is
re-entrant.
"""

from std.memory import Layout, alloc

from plegma.consensus import ConsensusOptions, consensus
from plegma.homopolymer import STRAND_UNKNOWN, HomopolymerModel
from plegma.kernel import BOTH, LEFT, RIGHT

comptime Bytes = Pointer[UInt8, MutUntrackedOrigin]
comptime Model = Pointer[HomopolymerModel, MutUntrackedOrigin]


struct Result(Movable):
    """A consensus call's outcome: its sequences (`parts`), or the message of what failed (`error`, NUL-ended)."""

    var parts: List[List[UInt8]]
    var error: List[UInt8]

    def __init__(out self):
        self.parts = List[List[UInt8]]()
        self.error = List[UInt8]()


def _box[T: Movable](var value: T) -> Pointer[T, MutUntrackedOrigin]:
    var allocation = alloc(Layout[T](length=1))
    allocation.unsafe_ptr().unsafe_write(value^)
    return allocation^.unsafe_leak()


def _free_model(pointer: Model):
    _ = pointer.unsafe_take_pointee()
    pointer.unsafe_free()


def _free_result(pointer: Pointer[Result, MutUntrackedOrigin]):
    _ = pointer.unsafe_take_pointee()
    pointer.unsafe_free()


def _text(pointer: Bytes, length: Int) -> String:
    return String(unsafe_from_utf8=Span(unsafe_ptr=pointer, length=length))


def _c_string_length(pointer: Bytes) -> Int:
    var n = 0
    while pointer[unsafe_offset=n] != 0:
        n += 1
    return n


def _write_error(message: String, error: OptionalPointer[UInt8, MutUntrackedOrigin], error_size: Int):
    """`message` into the caller's buffer, cut to fit, NUL-ended (nothing without a buffer)."""
    if not error or error_size <= 0:
        return
    var out = error.value()
    var bytes = message.as_bytes()
    var n = min(len(bytes), error_size - 1)
    for i in range(n):
        out[unsafe_offset=i] = bytes[i]
    out[unsafe_offset=n] = 0


@export
def plegma_model_load(
    path: Bytes, error: OptionalPointer[UInt8, MutUntrackedOrigin], error_size: Int
) abi("C") -> OptionalPointer[HomopolymerModel, MutUntrackedOrigin]:
    try:
        return _box(HomopolymerModel.load(_text(path, _c_string_length(path))))
    except e:
        _write_error(String(e), error, error_size)
        return None


@export
def plegma_model_parse(
    text: Bytes, length: Int, error: OptionalPointer[UInt8, MutUntrackedOrigin], error_size: Int
) abi("C") -> OptionalPointer[HomopolymerModel, MutUntrackedOrigin]:
    try:
        return _box(HomopolymerModel.parse(_text(text, length)))
    except e:
        _write_error(String(e), error, error_size)
        return None


@export
def plegma_model_free(model: OptionalPointer[HomopolymerModel, MutUntrackedOrigin]) abi("C"):
    if model:
        _free_model(model.value())


@export
def plegma_consensus(
    count: Int,
    sequences: Pointer[Bytes, MutUntrackedOrigin],
    lengths: Pointer[Int, MutUntrackedOrigin],
    anchors: OptionalPointer[Int32, MutUntrackedOrigin],
    strands: OptionalPointer[Int32, MutUntrackedOrigin],
    model: OptionalPointer[HomopolymerModel, MutUntrackedOrigin],
    polish_rounds: Int32,
    likely_rounds: Int32,
) abi("C") -> Pointer[Result, MutUntrackedOrigin]:
    var result = Result()
    try:
        var seqs = List[List[UInt8]](capacity=count)
        var anchor_of = List[Int](capacity=count)
        var strand_of = List[Int](capacity=count)
        for r in range(count):
            var source = sequences[unsafe_offset=r]
            var n = lengths[unsafe_offset=r]
            var seq = List[UInt8](capacity=n)
            for i in range(n):
                var c = source[unsafe_offset=i]
                seq.append(c - 32 if c >= UInt8(ord("a")) and c <= UInt8(ord("z")) else c)
            seqs.append(seq^)
            var anchor = Int(anchors.value()[unsafe_offset=r]) if anchors else BOTH
            if anchor != BOTH and anchor != LEFT and anchor != RIGHT:
                raise Error(t"read {r}: anchor {anchor} is not 0 (both), 1 (left) or 2 (right)")
            anchor_of.append(anchor)
            var strand = Int(strands.value()[unsafe_offset=r]) if strands else STRAND_UNKNOWN
            strand_of.append(strand if strand == 0 or strand == 1 else STRAND_UNKNOWN)
        var options = ConsensusOptions()
        options.polish_rounds = Int(polish_rounds)
        options.likely_rounds = Int(likely_rounds)
        if model:
            result.parts = consensus(seqs, anchor_of, strand_of, options, model.value()[])
        else:
            result.parts = consensus(seqs, anchor_of, strand_of, options)
    except e:
        result.parts.clear()
        result.error = List[UInt8](String(e).as_bytes())
        result.error.append(0)
    return _box(result^)


@export
def plegma_result_parts(result: Pointer[Result, MutUntrackedOrigin]) abi("C") -> Int:
    return len(result[].parts)


@export
def plegma_result_sequence(
    result: Pointer[Result, MutUntrackedOrigin], part: Int, length: OptionalPointer[Int, MutUntrackedOrigin]
) abi("C") -> OptionalPointer[UInt8, MutUntrackedOrigin]:
    if part < 0 or part >= len(result[].parts):
        if length:
            length.value()[] = 0
        return None
    if length:
        length.value()[] = len(result[].parts[part])
    return result[].parts[part].unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()


@export
def plegma_result_error(result: Pointer[Result, MutUntrackedOrigin]) abi("C") -> OptionalPointer[
    UInt8, MutUntrackedOrigin
]:
    if len(result[].error) == 0:
        return None
    return result[].error.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()


@export
def plegma_result_free(result: OptionalPointer[Result, MutUntrackedOrigin]) abi("C"):
    if result:
        _free_result(result.value())
