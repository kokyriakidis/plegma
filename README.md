# plegma

Fast partial order alignment (POA) of DNA in Mojo: multiple sequence alignment and consensus from a graph that
each sequence is woven into. The goal is to be the fastest POA aligner by a wide margin; rivals are
[abPOA](https://github.com/yangao07/abPOA), [SPOA](https://github.com/rvaser/spoa),
[POASTA](https://github.com/broadinstitute/poasta) and [Theseus](https://github.com/Marco-Sola-Lab/theseus-lib).

*πλέγμα* (plegma) is Greek for a woven mesh or lattice, and is what a POA graph is: sequences threaded into one
another where they agree.

Status: early. `build/plegma` (`pixi run build-cli`) runs a POA of one FASTA/FASTQ file with abPOA's flags and
outputs (consensus, MSA, GFA), in abPOA's global, local and extension modes, under its adaptive band; memory
follows the band rather than the graph times the sequence. In global mode the consensus is then polished
(`-P`): one-letter and whole-repeat-unit edits are kept while they lower the sequences' summed edit distance to
it, then while they raise the sequences' likelihood under a pair HMM, which mends most of the errors a POA
consensus makes in tandem repeats and short runs. `-I` refines the MSA (and the graph built from it) on its
sum of pairs: each sequence realigned in turn, exactly, to the others' columns, anchored by the sequences'
optimal pairwise alignments; on vntr-eval's 149 regions it beats abPOA in every stratum on both its unit and
affine scores. The package also holds the pairwise layer (sporos
anchors, optimal chaining, recursive backbones; `docs/design.md`, `experiments/`).

## Building

plegma builds with [pixi](https://pixi.sh), which fetches the pinned Mojo nightly. Its seeds come from
[sporos](https://github.com/kokyriakidis/sporos) and its pairwise alignment from
[karan](https://github.com/kokyriakidis/karan), both git submodules:

```sh
git clone --recurse-submodules https://github.com/kokyriakidis/plegma
pixi run test     # build and run the tests
pixi run build    # precompile the package to build/plegma.mojoc
pixi run format   # mojo format, 120 columns
```

## License

MIT
