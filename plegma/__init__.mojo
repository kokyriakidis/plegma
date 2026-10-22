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
  budget; `windows` -> `Window`, the gaps an exact aligner fills.
- Pairwise alignment through a backbone: `anchored_alignment` -> `Anchored`.
- The graph: `Graph` (`add`, `fuse` along a path's CIGAR, `fuse_nodes`,
  `msa`), `cigar_runs` -> `Runs`; windows aligned to it between two nodes,
  `align_window` -> `WindowAlignment` under `GapAffine`, and a whole sequence
  through a backbone, `graph_alignment`.
- Whole-graph alignment (abPOA's modes, scoring and adaptive band, memory
  following the band): `align_to_graph` -> `GraphAlignment` under `Scoring`
  (`GLOBAL`, `LOCAL`, `EXTENSION`), over a `Topology`; consensus by
  `Graph.heaviest_path`.
"""

from .alphabet import N_CODE, CodeBuffer, encode_base, encode_into
from .align import Anchored, anchored_alignment, graph_alignment
from .anchors import Anchor, find_anchors, mums, rare, rare_or_maximal, smems
from .backbone import ALL, Level, MUMS, RARE, RARE_OR_MAXIMAL, REFERENCE_UNIQUE, SMEMS, Window, backbone, windows
from .chain import Chain, ChainCosts, chain, linked_cost, without_overlaps
from .graph import Graph, Runs, cigar_runs
from .kernel import EXTENSION, GLOBAL, LOCAL, GraphAlignment, Scoring, Topology, align_to_graph
from .poa import GapAffine, SINK, SOURCE, WindowAlignment, align_window
