# plegma

Fast partial order alignment (POA) of DNA in Mojo: multiple sequence alignment and consensus. DNA only. The goal
is to be the fastest POA aligner by a wide margin; rivals are abPOA, SPOA, POASTA and Theseus (Marco-Sola-Lab/theseus-lib).

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
- Build and test on xaros: `rsync -az --exclude .pixi --exclude build --exclude .git ~/plegma/
  xaros.tail2f508.ts.net:plegma/`, then `~/.pixi/bin/pixi run test` there. Format there with `pixi run format`
  (120 columns) and rsync the `.mojo` files back.
- Toolchain: one Mojo nightly, pinned exactly in `pixi.toml` (1.2.0.dev2026102005, the same as sporos). Move it on
  purpose: change the pin, `pixi lock`, build and test on xaros.
- `plegma/alphabet.mojo` started as a copy of sporos's; the two may drift.
- Benchmark with min-of-N timings; other jobs share the machine.
- Finished work merges straight into `main` and is pushed; no pull requests.
