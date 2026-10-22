# MSA benchmarks: plegma against abPOA, SPOA, spoars, POASTA and Theseus

`simulate.py L N rate seed` writes a random reference of `L` letters (`ref_L.fa`) and `N` full-length reads of it
with errors at `rate` (substitutions, deletions and insertions alike; `reads_L_N_rate.fa`). `long.sh` runs plegma
`-S`, Theseus and abPOA on three such sets under a 64 GB address-space cap and a 1 h limit, and scores each
consensus by its edit distance to the reference with karan. Results, one file per run, are in `results/`.
