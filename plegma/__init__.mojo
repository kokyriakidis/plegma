"""plegma: fast partial order alignment of DNA in Mojo.

The public API is what this package exports; everything else in its modules
is implementation and may change.

- Sequences: `encode_base`, `encode_into`, `CodeBuffer` (ASCII to 2-bit
  codes, `N_CODE` for anything else).
- Anchors between two sequences: `find_anchors` -> `Anchor` (every MEM with
  its copies in each), and its views `mums`, `rare`, `smems`,
  `rare_or_maximal`.
- Chains: `chain` -> `Chain` under `ChainCosts` (MUMmer's `heaviest`,
  `chainx`, concave and period-aware indels, rarity-weighted anchors),
  `linked_cost`, `without_overlaps`.
- Backbones: `backbone` over `Level`s (`MUMS`, `REFERENCE_UNIQUE`, `RARE`,
  `SMEMS`, `ALL`, `RARE_OR_MAXIMAL`), seeding again inside every gap over a
  budget; `windows` -> `Window`, the gaps an exact aligner fills; `split_anchors`, anchors placed where they are
  sure (the flanks' MUMs, then the longest match unique in each window that leaves its length difference to one
  side and holds enough letters for it, recursively).
- Pairwise alignment through a backbone: `anchored_alignment` -> `Anchored`.
- The graph: `Graph` (`add`, `fuse` along a path's CIGAR, `fuse_nodes`,
  `msa`), `cigar_runs` -> `Runs`; windows aligned to it between two nodes,
  `align_window` -> `WindowAlignment` under `GapAffine`, and a whole sequence
  through a backbone, `graph_alignment`.
- Whole-graph alignment (abPOA's modes, scoring and adaptive band, memory
  following the band): `align_to_graph` -> `GraphAlignment` under `Scoring`
  (`GLOBAL`, `LOCAL`, `EXTENSION`), over a `Topology`, through a reusable
  `Workspace` for runs of alignments; consensus by
  `Graph.heaviest_path` (or `heaviest_path_of` a cluster's paths).
- Consensus polishing: `polish` (one-letter and whole-repeat-unit edits kept
  while they lower the sequences' summed edit distance to the consensus, then
  while they raise its likelihood under a pair HMM).
- MSA refinement: `refine` (each distinct sequence realigned exactly to the
  others' columns on their sum of pairs with gap opening, in a band centred on
  the sequences' optimal pairwise alignments; with `group_band`, every split of
  a UPGMA tree over the rows realigned as two profiles; mismatches, gap letters
  and gap openings priced apart), `objective` (what it lowers), `sum_of_pairs`,
  `gap_opens`.
- MSA by anchored progressive placement: `progressive` (each distinct sequence,
  longest first, placed exactly into the MSA so far through anchors every
  placed sequence votes for, MUMs then SMEMs then optimal alignments; refined;
  for a tandem array, also a start aligned copy by copy, the cheaper kept);
  its votes' costs `UNIT_VOTES`, `PRICED_VOTES`, `AFFINE_VOTES`, `TWO_PIECE_VOTES`; `on_reference`, the MSA
  anchored on its first sequence (each other aligned to it exactly), for calling variants against that reference.
- Tandem arrays: `tandem_repeat` -> `Repeat` (a sequence's array and period),
  `tandem_layout` -> `Layout` (a region's arrays cut into copies phased by one
  unit, the copies aligned as symbols on their sum of pairs, then as letters).
- Order and clusters: `guide_order` (abPOA's progressive order by minimizer
  multiset Jaccard), `cluster_rows` (several consensus sequences, by MSA
  columns where the sequences disagree).
"""

from .alphabet import N_CODE, CodeBuffer, encode_base, encode_into
from .align import Anchored, anchored_alignment, graph_alignment
from .anchors import Anchor, find_anchors, mums, rare, rare_or_maximal, smems
from .backbone import (
    ALL,
    Level,
    MUMS,
    RARE,
    RARE_OR_MAXIMAL,
    REFERENCE_UNIQUE,
    SMEMS,
    Window,
    backbone,
    split_anchors,
    windows,
)
from .chain import Chain, ChainCosts, chain, linked_cost, without_overlaps
from .clusters import cluster_rows
from .graph import Graph, Runs, cigar_runs
from .kernel import (
    EXTENSION,
    GLOBAL,
    LOCAL,
    GraphAlignment,
    Scoring,
    Topology,
    Workspace,
    align_between,
    align_to_graph,
)
from .order import guide_order
from .poa import GapAffine, SINK, SOURCE, WindowAlignment, align_window
from .polish import polish
from .progressive import AFFINE_VOTES, PRICED_VOTES, TWO_PIECE_VOTES, UNIT_VOTES, on_reference, progressive
from .refine import gap_opens, objective, refine, sum_of_pairs
from .units import Layout, Repeat, tandem_layout, tandem_repeat
