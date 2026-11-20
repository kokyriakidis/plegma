"""plegma from Python: the consensus of a read set (a local assembly's, between two anchors) through plegma's C API
(build/libplegma.so, `pixi run build-capi`), with ctypes alone.

    import plegma
    model = plegma.Model.load("HomopolymerModel-ONT-R10.4.1-sup5.0.0-2026-09-28.csv")  # optional
    parts = plegma.consensus(["ACGT...", "ACGT..."], anchors=["both", "left"], strands=[0, 1], model=model,
                             pinned=False)
    two = plegma.clustered_consensus(reads, model=model, max_clusters=2, min_freq=0.25, pinned=False)
    two.consensus, two.members, two.placements[0].cigar   # haplotypes, their reads, each read's alignment

The library is found at $PLEGMA_LIBRARY, else at ../build/libplegma.so beside this file. Calls are re-entrant and a
model is read-only once made, so threads may share one. Failures raise plegma.Error.
"""

import ctypes
import os

__all__ = [
    "Error",
    "Model",
    "Clustered",
    "Placement",
    "consensus",
    "clustered_consensus",
    "phased_consensus",
    "version",
    "BOTH",
    "LEFT",
    "RIGHT",
    "INSIDE",
    "FORWARD",
    "REVERSE",
    "UNKNOWN",
]

BOTH, LEFT, RIGHT, INSIDE = 0, 1, 2, 3
FORWARD, REVERSE, UNKNOWN = 0, 1, 2
_ANCHORS = {"both": BOTH, "left": LEFT, "right": RIGHT, "inside": INSIDE}


class Error(RuntimeError):
    """A plegma call that failed, with the library's message."""


class _Options(ctypes.Structure):
    """plegma_options (capi/plegma.h)."""

    _fields_ = [
        ("size", ctypes.c_uint32),
        ("polish_rounds", ctypes.c_int32),
        ("likely_rounds", ctypes.c_int32),
        ("pinned", ctypes.c_int32),
        ("min_cover", ctypes.c_int32),
        ("match", ctypes.c_int32),
        ("mismatch", ctypes.c_int32),
        ("gap_open1", ctypes.c_int32),
        ("gap_extend1", ctypes.c_int32),
        ("gap_open2", ctypes.c_int32),
        ("gap_extend2", ctypes.c_int32),
        ("band_b", ctypes.c_int32),
        ("band_f", ctypes.c_float),
        ("band_max", ctypes.c_int32),
        ("noise_run", ctypes.c_int32),
        ("two_lengths_db", ctypes.c_float),
    ]


_OPTION_NAMES = [name for name, _ in _Options._fields_ if name != "size"]


class _Placement(ctypes.Structure):
    """plegma_placement (capi/plegma.h)."""

    _fields_ = [
        ("cluster", ctypes.c_int64),
        ("cost", ctypes.c_int64),
        ("consensus_start", ctypes.c_int64),
        ("consensus_end", ctypes.c_int64),
        ("read_start", ctypes.c_int64),
        ("read_end", ctypes.c_int64),
    ]


def _library():
    path = os.environ.get("PLEGMA_LIBRARY") or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "build", "libplegma.so"
    )
    lib = ctypes.CDLL(path)
    p = ctypes.c_void_p
    lib.plegma_version.restype = ctypes.c_char_p
    lib.plegma_version.argtypes = []
    lib.plegma_options_init.restype = None
    lib.plegma_options_init.argtypes = [ctypes.POINTER(_Options)]
    lib.plegma_model_load.restype = p
    lib.plegma_model_load.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_size_t]
    lib.plegma_model_parse.restype = p
    lib.plegma_model_parse.argtypes = [ctypes.c_char_p, ctypes.c_size_t, ctypes.c_char_p, ctypes.c_size_t]
    lib.plegma_model_free.restype = None
    lib.plegma_model_free.argtypes = [p]
    lib.plegma_consensus.restype = p
    lib.plegma_consensus.argtypes = [
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_char_p),
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        p,
        ctypes.POINTER(_Options),
    ]
    reads_args = [
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_char_p),
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
    ]
    lib.plegma_clustered_consensus.restype = p
    lib.plegma_clustered_consensus.argtypes = reads_args + [p, ctypes.POINTER(_Options), ctypes.c_int32, ctypes.c_float]
    lib.plegma_phased_consensus.restype = p
    lib.plegma_phased_consensus.argtypes = reads_args + [ctypes.POINTER(ctypes.c_int32), p, ctypes.POINTER(_Options)]
    lib.plegma_result_members.restype = ctypes.c_size_t
    lib.plegma_result_members.argtypes = [p, ctypes.c_size_t, ctypes.POINTER(ctypes.POINTER(ctypes.c_int32))]
    lib.plegma_result_cigar.restype = ctypes.c_char_p
    lib.plegma_result_cigar.argtypes = [p, ctypes.c_size_t, ctypes.POINTER(_Placement)]
    lib.plegma_result_parts.restype = ctypes.c_size_t
    lib.plegma_result_parts.argtypes = [p]
    lib.plegma_result_sequence.restype = ctypes.POINTER(ctypes.c_char)
    lib.plegma_result_sequence.argtypes = [p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
    lib.plegma_result_error.restype = ctypes.c_char_p
    lib.plegma_result_error.argtypes = [p]
    lib.plegma_result_free.restype = None
    lib.plegma_result_free.argtypes = [p]
    return lib


_lib = _library()


def version():
    """The library's version, "major.minor.patch"."""
    return _lib.plegma_version().decode()


class Model:
    """A homopolymer model (shasta2's csv), read-only once made: `Model.load(path)` or `Model.parse(text)`."""

    def __init__(self, handle):
        self._handle = handle

    @classmethod
    def load(cls, path):
        error = ctypes.create_string_buffer(512)
        handle = _lib.plegma_model_load(os.fsencode(path), error, len(error))
        if not handle:
            raise Error(error.value.decode(errors="replace"))
        return cls(handle)

    @classmethod
    def parse(cls, text):
        data = text.encode() if isinstance(text, str) else bytes(text)
        error = ctypes.create_string_buffer(512)
        handle = _lib.plegma_model_parse(data, len(data), error, len(error))
        if not handle:
            raise Error(error.value.decode(errors="replace"))
        return cls(handle)

    def close(self):
        if self._handle:
            _lib.plegma_model_free(self._handle)
            self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __del__(self):
        self.close()


def _anchor(value):
    if isinstance(value, str):
        if value.lower() not in _ANCHORS:
            raise ValueError(f"anchor {value!r}: not 'both', 'left', 'right' or 'inside'")
        return _ANCHORS[value.lower()]
    return int(value)


def _int32(value, what):
    """`value` if a C int32_t holds it; ctypes would wrap it silently (2**32 read as 0)."""
    if not -(2**31) <= value < 2**31:
        raise ValueError(f"{what} {value}: out of a 32-bit integer's range")
    return value


def _int32s(values, what):
    return (ctypes.c_int32 * len(values))(*[_int32(int(v), what) for v in values])


class _Reads:
    """The arguments every consensus call takes, as the C API takes them (kept alive while the call runs)."""

    def __init__(self, reads, anchors, strands, model, options, caller):
        letters = [r.encode() if isinstance(r, str) else bytes(r) for r in reads]
        self.count = len(letters)
        if anchors is not None and len(anchors) != self.count:
            raise ValueError(f"{len(anchors)} anchors for {self.count} reads")
        if strands is not None and len(strands) != self.count:
            raise ValueError(f"{len(strands)} strands for {self.count} reads")
        self.settings = _Options()
        _lib.plegma_options_init(ctypes.byref(self.settings))
        for name, value in options.items():
            if name not in _OPTION_NAMES:
                raise TypeError(f"{caller}() got an unexpected option {name!r}")
            value = int(value) if name == "pinned" else value
            if isinstance(value, int):
                _int32(value, name)
            setattr(self.settings, name, value)
        self.sequences = (ctypes.c_char_p * self.count)(*letters)
        self.lengths = (ctypes.c_size_t * self.count)(*[len(s) for s in letters])
        self.anchors = None if anchors is None else _int32s([_anchor(a) for a in anchors], "anchor")
        self.strands = None if strands is None else _int32s([UNKNOWN if s is None else int(s) for s in strands], "strand")
        self.model = model._handle if model is not None else None

    def head(self):
        return (self.count, self.sequences, self.lengths, self.anchors, self.strands)


def _parts(result):
    error = _lib.plegma_result_error(result)
    if error:
        raise Error(error.decode(errors="replace"))
    parts = []
    for part in range(_lib.plegma_result_parts(result)):
        length = ctypes.c_size_t(0)
        data = _lib.plegma_result_sequence(result, part, ctypes.byref(length))
        parts.append(ctypes.string_at(data, length.value).decode() if length.value else "")
    return parts


def consensus(reads, anchors=None, strands=None, model=None, **options):
    """The consensus of `reads` (str or bytes; A, C, G, T), as `plegma poa -r 0` writes it: a list of one sequence,
    or of the left and the right side when no read spans both anchors and the two could not be joined.

    anchors: per read "both", "left" or "right" (or BOTH, LEFT, RIGHT); None: all spanning. strands: per read 0, 1
    or None (unknown), for the homopolymer model; None: all unknown. model: a Model. options: plegma_options' fields
    by name (polish_rounds, likely_rounds, pinned, min_cover, match, mismatch, gap_open1, gap_extend1, gap_open2,
    gap_extend2, band_b, band_f, band_max, noise_run, two_lengths_db; see capi/plegma.h), plegma poa's defaults
    otherwise.
    """
    given = _Reads(reads, anchors, strands, model, options, "consensus")
    result = _lib.plegma_consensus(*given.head(), given.model, ctypes.byref(given.settings))
    try:
        return _parts(result)
    finally:
        _lib.plegma_result_free(result)


class Placement:
    """A read's alignment to its cluster's consensus: the cluster, the alignment's cost (twice the options' scores,
    see capi/plegma.h), the consensus letters consensus_start ..< consensus_end and the read letters
    read_start ..< read_end it covers, and its CIGAR (=, X, I a letter the read has and the consensus not, D one the
    consensus has and the read not)."""

    __slots__ = ("cluster", "cost", "consensus_start", "consensus_end", "read_start", "read_end", "cigar")

    def __init__(self, raw, cigar):
        for name in self.__slots__[:-1]:
            setattr(self, name, getattr(raw, name))
        self.cigar = cigar

    def __repr__(self):
        return (
            f"Placement(cluster={self.cluster}, cost={self.cost}, consensus={self.consensus_start}:"
            f"{self.consensus_end}, read={self.read_start}:{self.read_end}, cigar={self.cigar!r})"
        )


class Clustered:
    """Several consensus sequences of one read set: `consensus` (one per cluster), `members` (each cluster's reads,
    input indices, ascending) and `placements` (one Placement per read, in input order)."""

    __slots__ = ("consensus", "members", "placements")

    def __init__(self, result, count):
        self.consensus = _parts(result)
        self.members = []
        for part in range(len(self.consensus)):
            ids = ctypes.POINTER(ctypes.c_int32)()
            n = _lib.plegma_result_members(result, part, ctypes.byref(ids))
            self.members.append([ids[i] for i in range(n)])
        self.placements = []
        for read in range(count):
            raw = _Placement()
            cigar = _lib.plegma_result_cigar(result, read, ctypes.byref(raw))
            self.placements.append(Placement(raw, cigar.decode() if cigar else ""))

    def __repr__(self):
        return f"Clustered({len(self.consensus)} clusters, {[len(m) for m in self.members]} reads)"


def clustered_consensus(reads, anchors=None, strands=None, model=None, max_clusters=2, min_freq=0.25, **options):
    """Several consensus sequences of `reads` (haplotypes, as a diploid caller asks abPOA's -d for): at most
    `max_clusters`, each the consensus of the reads nearest it, a cluster with under `min_freq` of the reads joining
    another. Reads, anchors (also "inside": a read on neither anchor, aligned within a cluster's consensus), strands,
    model and options as `consensus` takes them. A Clustered."""
    given = _Reads(reads, anchors, strands, model, options, "clustered_consensus")
    result = _lib.plegma_clustered_consensus(
        *given.head(), given.model, ctypes.byref(given.settings), int(max_clusters), float(min_freq)
    )
    try:
        return Clustered(result, given.count)
    finally:
        _lib.plegma_result_free(result)


def phased_consensus(reads, haplotypes, anchors=None, strands=None, model=None, **options):
    """The two haplotypes' consensus sequences of reads already phased: `haplotypes` per read 0 or 1 (-1 or None:
    neither, which builds on the haplotype whose letters it carries where the two differ, on both where it covers no
    such place, on neither where it carries as many of each, and is placed where it fits best). Where the two differ
    only by homopolymer lengths and a model is given, each such run is decided on both haplotypes' reads (their own
    lengths only when that explains the reads `two_lengths_db` better). Other arguments as `clustered_consensus`.
    A Clustered of two clusters."""
    given = _Reads(reads, anchors, strands, model, options, "phased_consensus")
    if len(haplotypes) != given.count:
        raise ValueError(f"{len(haplotypes)} haplotypes for {given.count} reads")
    haps = _int32s([-1 if h is None else int(h) for h in haplotypes], "haplotype")
    result = _lib.plegma_phased_consensus(*given.head(), haps, given.model, ctypes.byref(given.settings))
    try:
        return Clustered(result, given.count)
    finally:
        _lib.plegma_result_free(result)
