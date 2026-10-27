"""plegma from Python: the consensus of a read set (a local assembly's, between two anchors) through plegma's C API
(build/libplegma.so, `pixi run build-capi`), with ctypes alone.

    import plegma
    model = plegma.Model.load("HomopolymerModel-ONT-R10.4.1-sup5.0.0-2026-09-28.csv")  # optional
    parts = plegma.consensus(["ACGT...", "ACGT..."], anchors=["both", "left"], strands=[0, 1], model=model,
                             pinned=False)

The library is found at $PLEGMA_LIBRARY, else at ../build/libplegma.so beside this file. Calls are re-entrant and a
model is read-only once made, so threads may share one. Failures raise plegma.Error.
"""

import ctypes
import os

__all__ = ["Error", "Model", "consensus", "version", "BOTH", "LEFT", "RIGHT", "FORWARD", "REVERSE", "UNKNOWN"]

BOTH, LEFT, RIGHT = 0, 1, 2
FORWARD, REVERSE, UNKNOWN = 0, 1, 2
_ANCHORS = {"both": BOTH, "left": LEFT, "right": RIGHT}


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
            raise ValueError(f"anchor {value!r}: not 'both', 'left' or 'right'")
        return _ANCHORS[value.lower()]
    return int(value)


def consensus(reads, anchors=None, strands=None, model=None, **options):
    """The consensus of `reads` (str or bytes; A, C, G, T), as `plegma poa -r 0` writes it: a list of one sequence,
    or of the left and the right side when no read spans both anchors and the two could not be joined.

    anchors: per read "both", "left" or "right" (or BOTH, LEFT, RIGHT); None: all spanning. strands: per read 0, 1
    or None (unknown), for the homopolymer model; None: all unknown. model: a Model. options: plegma_options' fields
    by name (polish_rounds, likely_rounds, pinned, min_cover, match, mismatch, gap_open1, gap_extend1, gap_open2,
    gap_extend2, band_b, band_f, band_max, noise_run, two_lengths_db; see capi/plegma.h), plegma poa's defaults
    otherwise.
    """
    letters = [r.encode() if isinstance(r, str) else bytes(r) for r in reads]
    count = len(letters)
    if anchors is not None and len(anchors) != count:
        raise ValueError(f"{len(anchors)} anchors for {count} reads")
    if strands is not None and len(strands) != count:
        raise ValueError(f"{len(strands)} strands for {count} reads")
    settings = _Options()
    _lib.plegma_options_init(ctypes.byref(settings))
    for name, value in options.items():
        if name not in _OPTION_NAMES:
            raise TypeError(f"consensus() got an unexpected option {name!r}")
        setattr(settings, name, int(value) if name == "pinned" else value)
    sequences = (ctypes.c_char_p * count)(*letters)
    lengths = (ctypes.c_size_t * count)(*[len(s) for s in letters])
    anchor_array = None if anchors is None else (ctypes.c_int32 * count)(*[_anchor(a) for a in anchors])
    strand_array = (
        None if strands is None else (ctypes.c_int32 * count)(*[UNKNOWN if s is None else int(s) for s in strands])
    )
    result = _lib.plegma_consensus(
        count,
        sequences,
        lengths,
        anchor_array,
        strand_array,
        model._handle if model is not None else None,
        ctypes.byref(settings),
    )
    try:
        error = _lib.plegma_result_error(result)
        if error:
            raise Error(error.decode(errors="replace"))
        parts = []
        for part in range(_lib.plegma_result_parts(result)):
            length = ctypes.c_size_t(0)
            data = _lib.plegma_result_sequence(result, part, ctypes.byref(length))
            parts.append(ctypes.string_at(data, length.value).decode() if length.value else "")
        return parts
    finally:
        _lib.plegma_result_free(result)
