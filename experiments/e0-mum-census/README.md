# E0: pairwise MUM census on the vntr-eval regions

The first measurement behind `docs/design.md`, made while surveying vntr-eval (`docs/research/vntr-eval.md`,
section 5.4): how many pairwise MUMs fall in VNTR arrays and flanks, how many lie on some optimal alignment, and
what chaining them costs. It is plain C and Python, standard library only, and predates plegma's own code.

```sh
cc -O3 -o build/mumcheck experiments/e0-mum-census/mumcheck.c
python3 -I experiments/e0-mum-census/test_mumcheck.py build/mumcheck           # against brute force
python3 -I experiments/e0-mum-census/mum_driver.py <vntr-eval>/regions build/mumcheck build/mum20.tsv 7 20
python3 -I experiments/e0-mum-census/region_stats.py <vntr-eval>/regions build/region_stats.tsv
python3 -I experiments/e0-mum-census/mum_summary.py build/mum20.tsv build/region_stats.tsv
```

`region_stats.tsv` is the per-region size table of the 149 packages (vntr-eval da5c634). Check each script's usage
line before running it; they were written for one survey and may need small fixes.
