# Seeding strategies on pairs: are the anchors where the true or the optimal alignment puts them?

Every sequence joins the graph through one pairwise alignment, so a seeding strategy is judged pair by pair:
its backbone (anchors chained, trimmed, gaps between them filled exactly by karan) against the true alignment
where one is known and against the pair's optimum always.

- `simulate_vntr.py` makes VNTR pairs with their true alignment (flanks, copies of a unit with per-copy
  variants, whole-copy gains and losses, point mutations; every letter keeps its ancestral origin).
- `real_pairs.py` takes the vntr-eval packages' hap32 haplotypes: CHM13 against every other distinct one and
  random pairs, tagged with region, stratum, period and copies. No truth: scored against the optimum.
- `seeding.mojo` (`pixi run build-seeding`) runs strategies written as `levels[;options]` (see its header:
  `mum51,smem31,...`, `stay`, `s=`, `T=`, `n=`, `H=`) and writes one row per pair and strategy: cost against
  the optimum, F1 of aligned letter pairs against the truth and the optimum, anchors and how many lie on the
  truth, the largest window, time.
- `run.py` runs it in parallel chunks and summarises per strategy, and by period, copy variation and
  copy-number difference (simulated) or stratum, period and pair kind (real, `--meta`).

```sh
python3 -I experiments/seeding/simulate_vntr.py build/sim 400 11
python3 -I experiments/seeding/real_pairs.py <vntr-eval>/regions build/real 10 1
python3 -I experiments/seeding/run.py build/sim.fa build/sim.truth.tsv build/out_sim 9 "none" "mum51;stay" ...
python3 -I experiments/seeding/run.py build/real.fa - build/out_real 9 "none" "mum51;stay" ... --meta build/real.meta.tsv
```

Results are in `results/`.
