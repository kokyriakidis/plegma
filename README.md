# plegma

Partial order alignment (POA) and multiple sequence alignment of DNA in Mojo. Rivals are
[abPOA](https://github.com/yangao07/abPOA), [SPOA](https://github.com/rvaser/spoa),
[POASTA](https://github.com/broadinstitute/poasta) and [Theseus](https://github.com/Marco-Sola-Lab/theseus-lib).

*πλέγμα* (plegma) is Greek for a woven mesh or lattice, and is what a POA graph is: sequences threaded into one
another where they agree.

## Commands

`build/plegma` (`pixi run build-cli`) has two commands over one FASTA/FASTQ file, each writing a consensus, the MSA or
the graph as GFA (`-r`):

- **`plegma poa`**: the general POA, with abPOA's flags, defaults and output modes (global, local and extension
  alignment under its adaptive band; memory follows the band rather than the graph times the sequence). In global
  mode the consensus is then polished (`-P`): one-letter and whole-repeat-unit edits are kept while they lower the
  sequences' summed edit distance to it, then while they raise the sequences' likelihood under a pair HMM, which mends
  most of the errors a POA consensus makes in tandem repeats and short runs (13 errors on 180 real VNTR read sets,
  abPOA 357; 55 on simulated sets, abPOA 458). With `-S` each sequence is pinned to the consensus so far as
  minimap2 pins an alignment (its seeds and chain; pins only on k-mers unique in both) and the DP runs only between
  the pins: the same consensus with `-p` (13 and 58 errors), better allele MSAs (pgrealign qc pair F1 0.878 against
  0.864 on 377 chr20 regions), 1.4 to 7 times faster on the largest chr20 VNTRs.
- **`plegma ref`**: the MSA anchored on the first sequence, the reference (CHM13's allele in a pangenome region): every
  other aligned to it alone at minimap2's asm5 scores, gaps left-normalized, its letters in the reference's columns
  and the runs inserted at one place aligned to each other on the longest, so each sequence's alignment with the
  reference is its own and a copy-number change stays one event. With `-F` (the reference's flanks) each alignment is
  pinned as minimap2 pins it (`plegma/pins.mojo`, a port of its seeds, chain and pins at dipcall's settings, exact
  between the pins), so differences are written as dipcall-made truth sets write them; `-t` spreads the alignments
  over threads. Spliced into chr20 of the HPRC graph by pgrealign (`experiments/chr20`), `ref -F`'s graphs genotype
  HG002 better than pgrealign's star: structural F1 0.737 against 0.614, 0.703 against 0.556 in repeats, 0.452
  against 0.347 on vntr-eval's hard set (truvari refine, which scores haplotypes however their records are split:
  0.758 against 0.667, the hard set 0.476 against 0.385). On the largest chr20 VNTRs it is 17-19 times faster than
  abPOA in a tenth to a hundredth of its memory.

The library (`plegma/__init__.mojo` lists its API) holds both: the graph and its whole-graph aligner
(`plegma/kernel.mojo`), consensus polishing, and `on_reference` with minimap2's pinned alignment.

## Building

plegma builds with [pixi](https://pixi.sh), which fetches the pinned Mojo nightly. Its pairwise alignment comes from
[karan](https://github.com/kokyriakidis/karan) and its FASTA/FASTQ reading from
[sporos](https://github.com/kokyriakidis/sporos), both git submodules:

```sh
git clone --recurse-submodules https://github.com/kokyriakidis/plegma
pixi run test       # build and run the tests
pixi run build-cli  # the command line, build/plegma
pixi run fuzz       # random inputs through both commands, every output checked
pixi run format     # mojo format, 120 columns
```

## License

MIT
