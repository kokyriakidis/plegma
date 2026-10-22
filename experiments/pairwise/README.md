# E1-E5: pairwise anchors, chains and filled backbones

`pairwise.mojo` runs on one vntr-eval region: CHM13 is the partner and every other distinct hap32 sequence a
query.
- **E1, the census.** For each anchor kind (MUM, mumreference, rare at caps 2/8/50, SMEM, all MEMs) and minimum
  length (15/21/31/51), it records how many anchors fall in the repeat array. It then tests up to 200 of each
  against the optimum: an anchor lies on some optimal alignment when the unit costs before and after it add up to
  the pair's edit distance.
- **E2-E5, the backbones.** For each configuration (chainer, anchor kind, length, recursion, window budget), it
  builds the backbone and fills its windows exactly with karan, then compares the result with the optimum.

```sh
pixi run build-pairwise
python3 -I experiments/pairwise/run.py <vntr-eval>/regions <out dir> 8         # every region, 8 jobs
python3 -I experiments/pairwise/summarise.py <out dir> --md results.md          # per stratum
```

Results go in `results/`, one file per run, named for the plegma commit it ran at.
