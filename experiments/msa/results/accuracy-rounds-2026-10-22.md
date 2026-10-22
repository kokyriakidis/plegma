# Consensus accuracy with refinement rounds (2026-10-22)

Simulated sets of the benchmark (uniform and ONT-like errors, 500 bp to 100 kb, 10 to 200 reads), consensus
scored by edit distance to the truth; `-` where a tool failed (abPOA: memory or time). On the 63 sets every
configuration finished: abPOA 1,284 errors, abPOA -p 1,275, plegma -S -p 1,263, with one round of refinement
(-R 1) 1,140 (better than abPOA on 23 sets, worse on 8), -R 2 1,128.


| config | total errors | exact sets | failed | worst set | time s |
| :-- | --: | --: | --: | :-- | --: |
| abpoa | 1284 | 31/74 | 10 | sim_uni_L30000_N10_r15.fa (366) | 772.7 |
| abpoa_p | 1280 | 34/74 | 10 | sim_uni_L30000_N10_r15.fa (368) | 913.5 |
| plegma_Sp | 1611 | 32/74 | 0 | sim_uni_L30000_N10_r15.fa (345) | 135.3 |
| plegma_SpR1 | 1403 | 34/74 | 0 | sim_uni_L30000_N10_r15.fa (309) | 243.8 |
| plegma_SpR2 | 1396 | 32/74 | 0 | sim_uni_L30000_N10_r15.fa (303) | 345.3 |
| plegma_pR1 | 1415 | 35/74 | 0 | sim_uni_L30000_N10_r15.fa (310) | 944.5 |

| set | abpoa | abpoa_p | plegma_Sp | plegma_SpR1 | plegma_SpR2 | plegma_pR1 |
| :-- | --: | --: | --: | --: | --: | --: |
| sim_ont_L10000_N10_r05 | 4 | 3 | 4 | 4 | 4 | 4 |
| sim_ont_L10000_N10_r10 | 114 | 106 | 116 | 114 | 110 | 115 |
| sim_ont_L10000_N200_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_ont_L10000_N200_r10 | 21 | 27 | 24 | 27 | 24 | 26 |
| sim_ont_L10000_N50_r05 | 0 | 0 | 1 | 1 | 1 | 1 |
| sim_ont_L10000_N50_r10 | 37 | 31 | 37 | 32 | 32 | 33 |
| sim_ont_L2000_N10_r05 | 6 | 4 | 3 | 3 | 3 | 3 |
| sim_ont_L2000_N10_r10 | 29 | 27 | 24 | 23 | 23 | 23 |
| sim_ont_L2000_N200_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_ont_L2000_N200_r10 | 10 | 12 | 10 | 9 | 9 | 8 |
| sim_ont_L2000_N50_r05 | 1 | 1 | 1 | 1 | 1 | 1 |
| sim_ont_L2000_N50_r10 | 9 | 10 | 12 | 13 | 13 | 14 |
| sim_ont_L30000_N10_r05 | 32 | 35 | 32 | 30 | 31 | 30 |
| sim_ont_L30000_N10_r10 | 324 | 309 | 310 | 298 | 296 | 298 |
| sim_ont_L30000_N200_r05 | - | - | 3 | 4 | 4 | 4 |
| sim_ont_L30000_N200_r10 | - | - | 103 | 89 | 91 | 92 |
| sim_ont_L30000_N50_r05 | 7 | 8 | 7 | 7 | 7 | 7 |
| sim_ont_L30000_N50_r10 | - | - | 116 | 110 | 109 | 110 |
| sim_ont_L500_N10_r05 | 1 | 0 | 1 | 1 | 1 | 1 |
| sim_ont_L500_N10_r10 | 7 | 8 | 5 | 4 | 4 | 4 |
| sim_ont_L500_N200_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_ont_L500_N200_r10 | 1 | 1 | 1 | 0 | 0 | 0 |
| sim_ont_L500_N50_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_ont_L500_N50_r10 | 3 | 1 | 2 | 2 | 2 | 2 |
| sim_uni_L100000_N10_r01 | - | - | 0 | 0 | 0 | 0 |
| sim_uni_L100000_N10_r05 | - | - | 13 | 11 | 12 | 11 |
| sim_uni_L10000_N10_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L10000_N10_r05 | 1 | 1 | 1 | 1 | 1 | 1 |
| sim_uni_L10000_N10_r10 | 19 | 22 | 16 | 14 | 13 | 14 |
| sim_uni_L10000_N10_r15 | 97 | 133 | 139 | 122 | 118 | 122 |
| sim_uni_L10000_N200_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L10000_N200_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L10000_N200_r10 | 3 | 0 | 1 | 1 | 1 | 0 |
| sim_uni_L10000_N200_r15 | 16 | 15 | 12 | 10 | 10 | 9 |
| sim_uni_L10000_N50_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L10000_N50_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L10000_N50_r10 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L10000_N50_r15 | 27 | 15 | 21 | 12 | 13 | 12 |
| sim_uni_L2000_N10_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L2000_N10_r05 | 1 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L2000_N10_r10 | 6 | 7 | 4 | 3 | 3 | 3 |
| sim_uni_L2000_N10_r15 | 23 | 36 | 34 | 27 | 26 | 27 |
| sim_uni_L2000_N200_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L2000_N200_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L2000_N200_r10 | 1 | 0 | 2 | 0 | 2 | 0 |
| sim_uni_L2000_N200_r15 | 7 | 6 | 4 | 2 | 3 | 2 |
| sim_uni_L2000_N50_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L2000_N50_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L2000_N50_r10 | 0 | 1 | 1 | 1 | 1 | 1 |
| sim_uni_L2000_N50_r15 | 5 | 3 | 5 | 1 | 1 | 1 |
| sim_uni_L30000_N10_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L30000_N10_r05 | 4 | 5 | 4 | 3 | 3 | 3 |
| sim_uni_L30000_N10_r10 | 93 | 76 | 77 | 58 | 62 | 58 |
| sim_uni_L30000_N10_r15 | 366 | 368 | 345 | 309 | 303 | 310 |
| sim_uni_L30000_N200_r01 | 0 | - | 0 | 0 | 0 | 0 |
| sim_uni_L30000_N200_r05 | - | - | 0 | 0 | 0 | 0 |
| sim_uni_L30000_N200_r10 | - | - | 4 | 0 | 2 | 0 |
| sim_uni_L30000_N200_r15 | - | - | 37 | 12 | 14 | 16 |
| sim_uni_L30000_N50_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L30000_N50_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L30000_N50_r10 | - | 5 | 3 | 1 | 2 | 3 |
| sim_uni_L30000_N50_r15 | - | - | 69 | 36 | 34 | 39 |
| sim_uni_L500_N10_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L500_N10_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L500_N10_r10 | 3 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L500_N10_r15 | 6 | 4 | 5 | 5 | 5 | 5 |
| sim_uni_L500_N200_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L500_N200_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L500_N200_r10 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L500_N200_r15 | 0 | 0 | 0 | 1 | 1 | 1 |
| sim_uni_L500_N50_r01 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L500_N50_r05 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L500_N50_r10 | 0 | 0 | 0 | 0 | 0 | 0 |
| sim_uni_L500_N50_r15 | 0 | 0 | 2 | 1 | 1 | 1 |
