# Seeding strategies, run 2: anchors outside tandem repeats (2026-10-22)

The 3,623 real vntr-eval haplotype pairs of run 1, against each pair's optimum (F1 of aligned letter pairs /
cost over optimum / share at the optimum).

| strategy | all pairs | hotspot VNTR | hotspot other | matched controls | anchors / pair |
| :-- | :-- | :-- | :-- | :-- | --: |
| mum51;stay (run 1's best) | 0.8322 / 1.0215 / 0.54 | 0.6917 / 1.0250 / 0.26 | 0.8645 / 1.0162 / 0.51 | 0.9678 / 1.0080 / 0.87 | 8.7 |
| mum51;stay;V (verified) | 0.9302 / 1.0037 / 0.83 | 0.8705 / 1.0044 / 0.67 | 0.9576 / 1.0023 / 0.90 | 0.9769 / 1.0018 / 0.97 | 7.2 |
| mum51;stay;R=12:500 (filtered after chaining) | 0.9917 / 1.0004 / 0.97 | 0.9830 / 1.0005 / 0.95 | 0.9960 / 1.0003 / 0.98 | 0.9996 / 1.0003 / 0.99 | 1.5 |
| **mum51;stay;M=12:500 (masked before chaining)** | **0.9997 / 1.0000 / 1.00** | **0.9995 / 1.0000 / 1.00** | 1.0000 / 1.0000 / 1.00 | 1.0000 / 1.0000 / 1.00 | 1.2 |
| mum51,mum31,mum21;stay;M=12:500 (now -S) | 0.9996 / 1.0001 / 1.00 | 0.9995 / 1.0000 / 1.00 | 0.9999 / 1.0006 / 1.00 | 0.9996 / 1.0005 / 1.00 | 1.6 |
| ms51,ms31,ms21;stay;M=12:500 (MUMs and SMEMs) | 0.9959 / 1.0033 / 0.99 | 0.9958 / 1.0002 / 0.99 | 0.9973 / 1.0102 / 0.98 | 0.9933 / 1.0180 / 0.97 | 1.6 |
| ms51,ms31,ms21;stay;M=12:500;V | 0.9993 / 1.0000 / 1.00 | 0.9991 / 1.0000 / 1.00 | 1.0000 / 1.0000 / 1.00 | 0.9996 / 1.0005 / 1.00 | 1.6 |

The mask: a k-mer start is masked when the same 12-mer starts again within 500 letters of the partner; a seed
any of whose 12-mers is masked is dropped before chaining, at every level. Inside an array every exact match,
unique or not, may pair different copies; masked, the anchors fall in flanks and unique junctions and the DP
aligns each array whole (median largest window 3.5 kb). SMEMs still cost a little (up to 1.8% at the matched
controls) unless every anchor is verified.
