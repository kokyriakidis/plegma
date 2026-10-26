"""The C API of plegma (`capi/plegma.h`), built as `build/libplegma.so` (`pixi run build-capi`): the consensus of a read set
(`plegma.consensus`) for a C or C++ caller that assembles many, a homopolymer model loaded once and shared.

Handles are heap boxes the caller frees (`plegma_model_free`, `plegma_result_free`); nothing raises across the
boundary: a failure, or a wrong argument, comes back as a result holding its message (`plegma_result_error`), or a
null model with the message written into the caller's buffer. A model is read-only once made, so threads may share
one; every call is re-entrant.
"""

from std.math import isfinite
from std.memory import Layout, alloc

from plegma.consensus import ConsensusOptions, consensus
from plegma.homopolymer import STRAND_UNKNOWN, HomopolymerModel
from plegma.kernel import BOTH
from plegma.version import VERSION

comptime Bytes = Pointer[UInt8, MutUntrackedOrigin]
comptime Model = Pointer[HomopolymerModel, MutUntrackedOrigin]

comptime _VERSION = VERSION + "\x00"
"""`plegma_version`'s text, NUL-ended; the header's PLEGMA_VERSION_* agree (test_capi checks)."""


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


@fieldwise_init
struct _Field(TrivialRegisterPassable):
    """A plegma_options field: its byte offset (each field 4 bytes)."""

    var offset: Int


comptime _SIZE = _Field(0)
comptime _POLISH_ROUNDS = _Field(4)
comptime _LIKELY_ROUNDS = _Field(8)
comptime _PINNED = _Field(12)
comptime _MIN_COVER = _Field(16)
comptime _MATCH = _Field(20)
comptime _MISMATCH = _Field(24)
comptime _GAP_OPEN1 = _Field(28)
comptime _GAP_EXTEND1 = _Field(32)
comptime _GAP_OPEN2 = _Field(36)
comptime _GAP_EXTEND2 = _Field(40)
comptime _BAND_B = _Field(44)
comptime _BAND_F = _Field(48)
comptime _BAND_MAX = _Field(52)
comptime _OPTIONS_SIZE = 56
"""sizeof(plegma_options) in this library's plegma.h."""


struct _Options:
    """A caller's plegma_options: the fields its `size` covers; a field past it (a caller built against an older
    plegma.h) is left at its default."""

    var at: Bytes
    var size: Int

    def __init__(out self, at: Bytes):
        self.at = at
        self.size = Int(at.unsafe_bitcast[UInt32]()[])

    def has(self, field: _Field) -> Bool:
        return field.offset + 4 <= self.size

    def int(self, field: _Field, default: Int) -> Int:
        if not self.has(field):
            return default
        return Int(self.at.unsafe_offset(field.offset).unsafe_bitcast[Int32]()[])

    def float(self, field: _Field, default: Float64) -> Float64:
        if not self.has(field):
            return default
        return Float64(self.at.unsafe_offset(field.offset).unsafe_bitcast[Float32]()[])


def _set_int(at: Bytes, field: _Field, value: Int):
    at.unsafe_offset(field.offset).unsafe_bitcast[Int32]()[] = Int32(value)


def _settings(at: Bytes) raises -> ConsensusOptions:
    """The caller's plegma_options as `ConsensusOptions`, each value checked."""
    var given = _Options(at)
    if given.size < 4 or given.size > 4096:
        raise Error(t"plegma_options.size is {given.size}: set the options with plegma_options_init first")
    var options = ConsensusOptions()
    options.polish_rounds = given.int(_POLISH_ROUNDS, options.polish_rounds)
    options.likely_rounds = given.int(_LIKELY_ROUNDS, options.likely_rounds)
    options.poa.seeded = given.int(_PINNED, 1 if options.poa.seeded else 0) != 0
    options.poa.min_cover = given.int(_MIN_COVER, options.poa.min_cover)
    ref scoring = options.poa.scoring
    scoring.match_score = Int32(given.int(_MATCH, Int(scoring.match_score)))
    scoring.mismatch = Int32(given.int(_MISMATCH, Int(scoring.mismatch)))
    scoring.open1 = Int32(given.int(_GAP_OPEN1, Int(scoring.open1)))
    scoring.extend1 = Int32(given.int(_GAP_EXTEND1, Int(scoring.extend1)))
    scoring.open2 = Int32(given.int(_GAP_OPEN2, Int(scoring.open2)))
    scoring.extend2 = Int32(given.int(_GAP_EXTEND2, Int(scoring.extend2)))
    options.poa.band_b = given.int(_BAND_B, options.poa.band_b)
    options.poa.band_f = given.float(_BAND_F, options.poa.band_f)
    options.poa.band_max = given.int(_BAND_MAX, options.poa.band_max)
    if options.polish_rounds < 0 or options.likely_rounds < 0:
        raise Error(t"polishing rounds {options.polish_rounds},{options.likely_rounds}: each must be 0 or more")
    if options.poa.min_cover < 1:
        raise Error(t"min_cover {options.poa.min_cover}: must be 1 or more")
    if (
        scoring.match_score < 0
        or scoring.mismatch < 0
        or scoring.open1 < 0
        or scoring.open2 < 0
        or scoring.extend1 <= 0
        or scoring.extend2 < 0
    ):
        raise Error("scoring: match, mismatch and gap opens must be 0 or more, gap_extend1 more than 0")
    if not isfinite(options.poa.band_f) or options.poa.band_f < 0:
        raise Error(t"band_f {options.poa.band_f}: must be 0 or more")
    return options^


@export
def plegma_version() abi("C") -> Pointer[UInt8, ImmStaticOrigin]:
    return _VERSION.ptr()


@export
def plegma_options_init(settings: Bytes) abi("C"):
    var defaults = ConsensusOptions()
    settings.unsafe_bitcast[UInt32]()[] = UInt32(_OPTIONS_SIZE)
    _set_int(settings, _POLISH_ROUNDS, defaults.polish_rounds)
    _set_int(settings, _LIKELY_ROUNDS, defaults.likely_rounds)
    _set_int(settings, _PINNED, 1 if defaults.poa.seeded else 0)
    _set_int(settings, _MIN_COVER, defaults.poa.min_cover)
    ref scoring = defaults.poa.scoring
    _set_int(settings, _MATCH, Int(scoring.match_score))
    _set_int(settings, _MISMATCH, Int(scoring.mismatch))
    _set_int(settings, _GAP_OPEN1, Int(scoring.open1))
    _set_int(settings, _GAP_EXTEND1, Int(scoring.extend1))
    _set_int(settings, _GAP_OPEN2, Int(scoring.open2))
    _set_int(settings, _GAP_EXTEND2, Int(scoring.extend2))
    _set_int(settings, _BAND_B, defaults.poa.band_b)
    settings.unsafe_offset(_BAND_F.offset).unsafe_bitcast[Float32]()[] = Float32(defaults.poa.band_f)
    _set_int(settings, _BAND_MAX, defaults.poa.band_max)


@export
def plegma_model_load(
    path: OptionalPointer[UInt8, MutUntrackedOrigin], error: OptionalPointer[UInt8, MutUntrackedOrigin], error_size: Int
) abi("C") -> OptionalPointer[HomopolymerModel, MutUntrackedOrigin]:
    if not path:
        _write_error("plegma_model_load: no path", error, error_size)
        return None
    try:
        return _box(HomopolymerModel.load(_text(path.value(), _c_string_length(path.value()))))
    except e:
        _write_error(String(e), error, error_size)
        return None


@export
def plegma_model_parse(
    text: OptionalPointer[UInt8, MutUntrackedOrigin],
    length: Int,
    error: OptionalPointer[UInt8, MutUntrackedOrigin],
    error_size: Int,
) abi("C") -> OptionalPointer[HomopolymerModel, MutUntrackedOrigin]:
    if not text and length > 0:
        _write_error("plegma_model_parse: no text", error, error_size)
        return None
    try:
        return _box(HomopolymerModel.parse(_text(text.value(), length) if text else String()))
    except e:
        _write_error(String(e), error, error_size)
        return None


@export
def plegma_model_free(model: OptionalPointer[HomopolymerModel, MutUntrackedOrigin]) abi("C"):
    if model:
        _free_model(model.value())


def _letter(c: UInt8) -> UInt8:
    """A read's byte as the consensus takes it: A, C, G, T (either case) as upper case, anything else N."""
    var upper = c & 0xDF
    if upper == 0x41 or upper == 0x43 or upper == 0x47 or upper == 0x54:
        return upper
    return 0x4E


def _consensus(
    count: Int,
    sequences: OptionalPointer[OptionalPointer[UInt8, MutUntrackedOrigin], MutUntrackedOrigin],
    lengths: OptionalPointer[Int, MutUntrackedOrigin],
    anchors: OptionalPointer[Int32, MutUntrackedOrigin],
    strands: OptionalPointer[Int32, MutUntrackedOrigin],
    model: OptionalPointer[HomopolymerModel, MutUntrackedOrigin],
    settings: OptionalPointer[UInt8, MutUntrackedOrigin],
) raises -> List[List[UInt8]]:
    if count < 0:
        raise Error(t"count {count}: must be 0 or more")
    if count > 0 and (not sequences or not lengths):
        raise Error("sequences and lengths must be given for count > 0")
    var seqs = List[List[UInt8]](capacity=count)
    var anchor_of = List[Int](capacity=count)
    var strand_of = List[Int](capacity=count)
    for r in range(count):
        var n = lengths.value()[unsafe_offset=r]
        var source = sequences.value()[unsafe_offset=r]
        if n < 0 or (n > 0 and not source):
            raise Error(t"read {r}: no letters for length {n}")
        var seq = List[UInt8](capacity=n)
        for i in range(n):
            seq.append(_letter(source.value()[unsafe_offset=i]))
        seqs.append(seq^)
        anchor_of.append(Int(anchors.value()[unsafe_offset=r]) if anchors else BOTH)
        strand_of.append(Int(strands.value()[unsafe_offset=r]) if strands else STRAND_UNKNOWN)
    var options = _settings(settings.value()) if settings else ConsensusOptions()
    if model:
        return consensus(seqs^, anchor_of, strand_of, options, model.value()[])
    return consensus(seqs^, anchor_of, strand_of, options)


@export
def plegma_consensus(
    count: Int,
    sequences: OptionalPointer[OptionalPointer[UInt8, MutUntrackedOrigin], MutUntrackedOrigin],
    lengths: OptionalPointer[Int, MutUntrackedOrigin],
    anchors: OptionalPointer[Int32, MutUntrackedOrigin],
    strands: OptionalPointer[Int32, MutUntrackedOrigin],
    model: OptionalPointer[HomopolymerModel, MutUntrackedOrigin],
    settings: OptionalPointer[UInt8, MutUntrackedOrigin],
) abi("C") -> Pointer[Result, MutUntrackedOrigin]:
    var result = Result()
    try:
        result.parts = _consensus(count, sequences, lengths, anchors, strands, model, settings)
    except e:
        result.parts.clear()
        result.error = List[UInt8](String(e).as_bytes())
        result.error.append(0)
    return _box(result^)


@export
def plegma_result_parts(result: OptionalPointer[Result, MutUntrackedOrigin]) abi("C") -> Int:
    return len(result.value()[].parts) if result else 0


@export
def plegma_result_sequence(
    result: OptionalPointer[Result, MutUntrackedOrigin], part: Int, length: OptionalPointer[Int, MutUntrackedOrigin]
) abi("C") -> OptionalPointer[UInt8, MutUntrackedOrigin]:
    if not result or part < 0 or part >= len(result.value()[].parts):
        if length:
            length.value()[] = 0
        return None
    ref parts = result.value()[].parts
    if length:
        length.value()[] = len(parts[part])
    return parts[part].unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()


@export
def plegma_result_error(result: OptionalPointer[Result, MutUntrackedOrigin]) abi("C") -> OptionalPointer[
    UInt8, MutUntrackedOrigin
]:
    if not result or len(result.value()[].error) == 0:
        return None
    return result.value()[].error.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()


@export
def plegma_result_free(result: OptionalPointer[Result, MutUntrackedOrigin]) abi("C"):
    if result:
        _free_result(result.value())
