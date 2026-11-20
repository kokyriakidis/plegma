"""plegma: partial order alignment and multiple sequence alignment of DNA in Mojo.

The public API is what this package exports; everything else in its modules
is implementation and may change.

- Sequences: `encode_base`, `encode_into`, `CodeBuffer` (ASCII to 2-bit
  codes, `N_CODE` for anything else).
- The graph: `Graph` (`add`, `fuse_nodes` along a graph alignment's nodes,
  `from_msa`, `msa`, `heaviest_path`), `cigar_runs` -> `Runs`.
- Whole-graph alignment (abPOA's modes, scoring and adaptive band, memory
  following the band): `align_to_graph` -> `GraphAlignment` under `Scoring`
  (`GLOBAL`, `LOCAL`, `EXTENSION`; for reads fixed at one end of a local
  assembly `PREFIX` and `SUFFIX`), over a `Topology`, through a reusable
  `Workspace` for runs of alignments, between two nodes (`align_between`;
  `SOURCE` and `SINK` for the graph's ends).
- Pinned alignment to the graph: `align_pinned` (minimap2's pins between a
  sequence and the consensus path, each on a k-mer unique in both, the
  kernel only between them), over `window_nodes` with a `WindowScratch`.
- Consensus polishing: `polish` (one-letter and whole-repeat-unit edits kept
  while they lower the sequences' summed edit distance to the consensus, then
  while they raise its likelihood under a pair HMM), optionally with a
  `HomopolymerModel` (shasta2's empirical run-length error model, `load`ed
  from its csv): every run then at its most likely length given the reads'.
- The POA of a read set: `poa` -> `Poa` under `PoaOptions`, each read
  anchored `BOTH`, `LEFT` or `RIGHT` (a local assembly's reads: spanning it,
  or fixed at one anchor; with none spanning, each side assembled and the two
  joined, or one side alone, one-directional).
- One call from reads to consensus: `consensus` under `ConsensusOptions`
  (what `plegma poa -r 0` writes, with a `HomopolymerModel` as `-H`), the
  function behind the C API (`capi/plegma.h`, `capi/plegma.hpp`,
  `build/libplegma.so`) a local assembler calls instead of the command.
- Several consensus sequences of one read set (haplotypes, as a diploid caller asks abPOA's `-d` for):
  `clustered_consensus` -> `Clustered` (per cluster its consensus and members; per read its `Placement`, the
  CIGAR to its cluster's consensus at `alignment_costs` of the scoring), a read covering neither anchor `INSIDE`;
  `phased_consensus` the same for reads already phased into two haplotypes.
- Order and clusters: `guide_order` (abPOA's progressive order by minimizer
  multiset Jaccard), `cluster_rows` (several consensus sequences, by MSA
  columns where the sequences disagree).
- The MSA anchored on a reference: `on_reference` (each other sequence aligned
  to the first exactly, or pinned as minimap2 pins it given the reference's
  flanks, its letters in the reference letters' columns), for calling variants
  against that reference.
- `VERSION`, as `plegma --version` and the C API's `plegma_version` give it.
- minimap2's pinned alignment (2.31, -x asm5 as dipcall runs it): `pins` ->
  `Pins` (the points its alignment of a sequence to the reference passes
  through), `pins_in_context` (the same with flanks as the seeds' context),
  `pinned_pairing` (the alignment exact between them).
"""

from .alphabet import N_CODE, CodeBuffer, encode_base, encode_into
from .clustered import clustered_consensus
from .clusters import cluster_rows
from .consensus import ConsensusOptions, consensus
from .graph import Graph, Runs, cigar_runs
from .homopolymer import STRAND_UNKNOWN, HomopolymerModel
from .kernel import (
    BOTH,
    EXTENSION,
    GLOBAL,
    INSIDE,
    LEFT,
    LOCAL,
    PREFIX,
    RIGHT,
    SINK,
    SOURCE,
    SUFFIX,
    GraphAlignment,
    Scoring,
    Topology,
    Workspace,
    align_between,
    align_to_graph,
)
from .order import guide_order
from .phased import phased_consensus
from .pinned import WindowScratch, align_pinned, window_nodes
from .pins import Pins, pinned_pairing, pins, pins_in_context
from .placement import Clustered, Placement, alignment_costs
from .poa import Poa, PoaOptions, poa
from .polish import polish
from .reference import on_reference
from .version import VERSION
