"""plegma: partial order alignment and multiple sequence alignment of DNA in Mojo.

The public API is what this package exports; everything else in its modules
is implementation and may change.

- Sequences: `encode_base`, `encode_into`, `CodeBuffer` (ASCII to 2-bit
  codes, `N_CODE` for anything else).
- The graph: `Graph` (`add`, `fuse_nodes` along a graph alignment's nodes,
  `from_msa`, `msa`, `heaviest_path`), `cigar_runs` -> `Runs`.
- Whole-graph alignment (abPOA's modes, scoring and adaptive band, memory
  following the band): `align_to_graph` -> `GraphAlignment` under `Scoring`
  (`GLOBAL`, `LOCAL`, `EXTENSION`), over a `Topology`, through a reusable
  `Workspace` for runs of alignments, between two nodes (`align_between`;
  `SOURCE` and `SINK` for the graph's ends).
- Consensus polishing: `polish` (one-letter and whole-repeat-unit edits kept
  while they lower the sequences' summed edit distance to the consensus, then
  while they raise its likelihood under a pair HMM).
- Order and clusters: `guide_order` (abPOA's progressive order by minimizer
  multiset Jaccard), `cluster_rows` (several consensus sequences, by MSA
  columns where the sequences disagree).
- The MSA anchored on a reference: `on_reference` (each other sequence aligned
  to the first exactly, or pinned as minimap2 pins it given the reference's
  flanks, its letters in the reference letters' columns), for calling variants
  against that reference.
- minimap2's pinned alignment (2.31, -x asm5 as dipcall runs it): `pins` ->
  `Pins` (the points its alignment of a sequence to the reference passes
  through), `pins_in_context` (the same with flanks as the seeds' context),
  `pinned_pairing` (the alignment exact between them).
"""

from .alphabet import N_CODE, CodeBuffer, encode_base, encode_into
from .clusters import cluster_rows
from .graph import Graph, Runs, cigar_runs
from .kernel import (
    EXTENSION,
    GLOBAL,
    LOCAL,
    SINK,
    SOURCE,
    GraphAlignment,
    Scoring,
    Topology,
    Workspace,
    align_between,
    align_to_graph,
)
from .order import guide_order
from .pins import Pins, pinned_pairing, pins, pins_in_context
from .polish import polish
from .reference import on_reference
