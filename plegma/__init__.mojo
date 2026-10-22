"""plegma: fast partial order alignment of DNA in Mojo.

The public API is what this package exports; everything else in its modules
is implementation and may change.

- Sequences: `encode_base`, `encode_into`, `CodeBuffer` (ASCII to 2-bit
  codes, `N_CODE` for anything else).
- Anchors between two sequences: `find_anchors` -> `Anchor` (every MEM with
  its copies in each), and its views `mums`, `rare`, `smems`.
- Chains: `chain` -> `Chain` under `ChainCosts` (MUMmer's `heaviest`,
  `chainx`, concave and period-aware indels, rarity-weighted anchors),
  `linked_cost`, `without_overlaps`.
- Backbones: `backbone` over `Level`s (`MUMS`, `REFERENCE_UNIQUE`, `RARE`,
  `SMEMS`, `ALL`), recursing into gaps over a budget; `windows` -> `Window`,
  the gaps an exact aligner fills.
"""

from .alphabet import N_CODE, CodeBuffer, encode_base, encode_into
from .anchors import Anchor, find_anchors, mums, rare, smems
from .backbone import ALL, Level, MUMS, RARE, REFERENCE_UNIQUE, SMEMS, Window, backbone, windows
from .chain import Chain, ChainCosts, chain, linked_cost, without_overlaps
