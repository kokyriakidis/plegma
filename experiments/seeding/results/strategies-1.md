# Seeding strategies, run 1 (2026-10-22)

Two-piece costs (4, 6, 2, 24, 1), ties left, anchors trimmed 10 letters at each end unless `T=0`, unanchored
letters charged 0.1 in the chain unless `s=1`, chain look-back 200.

## 3,623 real haplotype pairs (vntr-eval hap32), against each pair's optimum

| strategy | F1 vs optimum | cost / optimum | at optimum | hotspot VNTR: F1 / cost / at opt |
| :-- | --: | --: | --: | --: |
| none (karan end to end) | 1.0000 | 1.0000 | 1.00 | 1.0000 / 1.0000 / 1.00 |
| mum51;stay | 0.8322 | 1.0215 | 0.54 | 0.6917 / 1.0250 / 0.26 |
| mum51,mum31,mum21,mum15;stay | 0.8077 | 1.0327 | 0.49 | 0.6640 / 1.0366 / 0.21 |
| smem51..15;stay | 0.7631 | 1.0597 | 0.38 | 0.6060 / 1.0643 / 0.13 |
| ms51..15;stay (MUMs and SMEMs) | 0.7631 | 1.0597 | 0.38 | 0.6060 / 1.0643 / 0.13 |
| mum51,smem51,...,mum15,smem15;stay (MUMs first) | 0.7616 | 1.0561 | 0.39 | 0.6020 / 1.0603 / 0.14 |
| combo8_51..15;stay (the -S seeds) | 0.7485 | 1.0685 | 0.36 | 0.5871 / 1.0722 / 0.10 |
| ms51..15;stay;T=0 | 0.7540 | 1.0703 | 0.37 | 0.5980 / 1.0738 / 0.13 |
| ms51..15;stay;s=1 | 0.7461 | 1.1458 | 0.36 | 0.5747 / 1.1695 / 0.10 |

Non-repeat SV controls: every strategy at the optimum. In VNTRs every exact match shorter than the array,
unique or not, can be false homology between copies: the fewest, longest, unique anchors (MUMs of 51) lose
least; SMEMs (the longest match, whatever copy it lands on) and re-seeding to 15 lose more.

## 400 simulated VNTR pairs, against the truth

The optimum itself has F1 0.856 against the truth (0.68 with identical copies, where every placement of a
copy-number gap ties and the truth picks one arbitrarily): F1 against the truth is informative where copies
vary. MUMs of 51 come closest (0.863-0.870 for periods over 18), SMEMs lower it (0.834-0.845); trimming
raises the share of anchors on the truth from 0.61 to 0.79.
