"""plegma: fast partial order alignment of DNA in Mojo.

The public API is what this package exports; everything else in its modules
is implementation and may change.

- Sequences: `encode_base`, `encode_into`, `CodeBuffer` (ASCII to 2-bit
  codes, `N_CODE` for anything else).
"""

from .alphabet import N_CODE, CodeBuffer, encode_base, encode_into
