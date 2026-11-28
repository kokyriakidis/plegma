# plegma

Partial order alignment (POA) of DNA in Mojo: multiple sequence alignment and consensus. DNA only. The goal is to be
the fastest, the leanest in memory and the most accurate POA aligner; rivals are abPOA, SPOA, POASTA and Theseus
(Marco-Sola-Lab/theseus-lib). The first job is re-aligning the VNTRs the HPRC Minigraph-Cactus graph aligns wrongly
(test regions and scoring: github.com/benedictpaten/vntr-eval; the realignment pipeline it feeds, and plegma's closest
rival there: github.com/benedictpaten/pgrealign).

Two commands over one library: `plegma poa`, the general POA (each sequence aligned to the whole graph, abPOA's
modes and adaptive band, `plegma/kernel.mojo`; consensus polished, `plegma/polish.mojo`), and `plegma ref`, the MSA
anchored on a reference for pgrealign (`plegma/reference.mojo`: each sequence aligned to the reference alone at
minimap2's asm5 costs, pinned as minimap2 pins it with the reference's flanks, `plegma/pins.mojo`). `docs/design.md`
maps the code; `docs/research/` holds the rival tools' studies.

## Dependencies

`deps/karan` (exact pairwise alignment) and `deps/sporos` (FASTA/FASTQ reading) are git submodules of the user's
own libraries; every build puts both on the include path. Change them where plegma needs it: commit and push in
the library first (its own tests passing), then move the submodule here. All three pin the same Mojo nightly.

## Mojo skills

Modular's Mojo skills live in `.claude/skills` and correct pretrained Mojo, which is out of date. Hooks in
`.claude/settings.json` inject them: the full `mojo-syntax` skill and `references/idiomatic-mojo.md` at session
start (also after `/clear` and compaction) and into every subagent, and a reminder on every prompt
(`PLEGMA_SKILLS_EVERY_PROMPT=full` injects the full text on every prompt instead). Follow them over pretrained
habits:

- `mojo-syntax` before writing, editing or reviewing any Mojo, with `references/idiomatic-mojo.md`.
- `mojo-gpu-fundamentals` before touching GPU code; `mojo-python-interop` before Mojo-Python interop.
- Include the skills rule in every subagent prompt.

`pixi run sync-skills` refreshes them from upstream; `.claude/skills/UPSTREAM` records the commit.

## Working here

- The public API is what `plegma/__init__.mojo` exports; keep its docstring and the README in step when adding or
  renaming anything public.
- Build and test on xaros: `rsync -az --delete --exclude .pixi --exclude build --exclude .git ~/plegma/
  xaros.tail2f508.ts.net:plegma/` (this carries the submodules' files), then `~/.pixi/bin/pixi run test` there.
  Format there with `pixi run format` (120 columns) and rsync the `.mojo` files back.
- Toolchain: one Mojo nightly, pinned exactly in `pixi.toml` with its `max` (1.3.0.dev2026112705, the same as sporos
  and karan). Move it on purpose: change the pin, `pixi lock`, build and test on xaros.
- `plegma/alphabet.mojo` started as a copy of sporos's; the two may drift.
- Benchmark with min-of-N timings; other jobs share the machine.
- Finished work merges straight into `main` and is pushed; no pull requests.
