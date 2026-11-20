"""Tracing: the decisions plegma takes, on stderr, for whoever needs to see why a consensus or an MSA came out as it
did.

`PLEGMA_TRACE` names the topics, comma-separated, or `all`:

- `poa`: the reads by anchor and their order, a read reverse-complemented (`-s`), a one-sided read found to be another
  copy's, the two sides of an assembly with no spanning read and how (and where) they were joined, a one-sided
  consensus cut.
- `pins`: the pins each pinned alignment keeps (chain length, pins, windows).
- `polish`: every round's edits and cost, the reads left out as outliers, the run lengths a homopolymer model changed.
- `clustered`: the clusters `clustered_consensus` starts from, each merge and why, the reads moved per refinement.
- `phased`: the places the two haplotypes differ and where each unlabelled read went.
- `reference`: per sequence of `on_reference`, how it was paired with the reference, and the fallback.

Each line is written whole (`plegma[topic] ...`), so lines from threads (`plegma ref -t`, a caller's threads) do not
interleave within a line. The variable is read once per top-level call (`Trace()`), and tracing stays out of inner
loops: off, it costs one branch per traced decision.
"""

from std.os import getenv
from std.sys import stderr

comptime POA = 1
comptime PINS = 2
comptime POLISH = 4
comptime CLUSTERED = 8
comptime PHASED = 16
comptime REFERENCE = 32
comptime _ALL = 63


def _name(topic: Int) -> StaticString:
    if topic == POA:
        return "poa"
    if topic == PINS:
        return "pins"
    if topic == POLISH:
        return "polish"
    if topic == CLUSTERED:
        return "clustered"
    if topic == PHASED:
        return "phased"
    return "reference"


struct Trace(ImplicitlyCopyable):
    """The topics `PLEGMA_TRACE` turns on, read when made."""

    var topics: Int

    def __init__(out self):
        self.topics = 0
        var wanted = getenv("PLEGMA_TRACE")
        if wanted.byte_length() == 0:
            return
        for word in wanted.split(","):
            var name = String(word.strip())
            if name == "all" or name == "1":
                self.topics = _ALL
            for topic in [POA, PINS, POLISH, CLUSTERED, PHASED, REFERENCE]:
                if name == _name(topic):
                    self.topics |= topic

    def __init__(out self, topics: Int):
        """Tracing exactly `topics` (a test's, say), whatever the environment says."""
        self.topics = topics

    @always_inline
    def on(self, topic: Int) -> Bool:
        return self.topics & topic != 0

    def line(self, topic: Int, message: String):
        """`message` as one line under `topic`, when it is on."""
        if self.on(topic):
            print(String("plegma[", _name(topic), "] ", message, "\n"), end="", file=stderr)
