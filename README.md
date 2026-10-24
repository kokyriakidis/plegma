# plegma

Fast partial order alignment (POA) of DNA in Mojo: multiple sequence alignment and consensus from a graph that
each sequence is woven into. The goal is to be the fastest POA aligner by a wide margin; rivals are
[abPOA](https://github.com/yangao07/abPOA), [SPOA](https://github.com/rvaser/spoa),
[POASTA](https://github.com/broadinstitute/poasta) and [Theseus](https://github.com/Marco-Sola-Lab/theseus-lib).

*πλέγμα* (plegma) is Greek for a woven mesh or lattice, and is what a POA graph is: sequences threaded into one
another where they agree.

Status: early. `build/plegma` (`pixi run build-cli`) runs a POA of one FASTA/FASTQ file with abPOA's flags and outputs
(consensus, MSA, GFA), in abPOA's global, local and extension modes, under its adaptive band; memory follows the band
rather than the graph times the sequence. In global mode the consensus is then polished (`-P`): one-letter and
whole-repeat-unit edits are kept while they lower the sequences' summed edit distance to it, then while they raise the
sequences' likelihood under a pair HMM, which mends most of the errors a POA consensus makes in tandem repeats and short
runs. `-I` refines the MSA (and the graph built from it) on its sum of pairs: each sequence realigned in turn, exactly,
to the others' columns, anchored by the sequences' optimal pairwise alignments; on vntr-eval's 149 regions it beats
abPOA in every stratum on both its unit and affine scores. `-A` builds the MSA without the graph, anchored on the first
sequence (the reference, CHM13's allele in a pangenome region): every other aligned to it exactly at pgrealign's scores,
gaps left-normalized, and the runs different sequences insert at one place aligned to each other, so each sequence's
alignment with the reference is its own optimal one and a copy-number change stays one event; spliced into chr20 of the
HPRC graph by pgrealign, its graphs genotype HG002 better than pgrealign's own (structural F1 0.658 against 0.614, 0.606
against 0.556 in repeats; small 0.9704 against 0.9691; `experiments/chr20`). `-A -N` instead places each distinct
sequence, longest first, exactly into the MSA so far in a band centred on the columns the placed sequences' optimal
alignments to it (at vntr-eval's affine costs) vote for (an anchor audit, `experiments/anchors`, showed seeds inside
repeats pair letters no optimal alignment pairs), and, where every sequence holds one clear tandem array, a second start
whose arrays are cut into copies phased by one unit and aligned copy by copy; each is refined (`-I`, with `-G`'s group
moves, which realign every split of a tree over the rows as two profiles), knowing the repeat's period (bands reaching a
copy, group moves splitting the sequences by copy number), on the sum of pairs over every pair of input sequences with a
gap opening priced at half a mismatch (`-W`), and the lowest kept; on vntr-eval's 149 regions its medians beat the
unit-aware realigner's on both its unit and affine costs in every stratum (hotspot VNTRs 1.032 / 1.080 times the
pairwise optima, against 1.048 / 1.093; pgrealign's star 1.225 / 1.204). The package also holds the pairwise layer
(sporos anchors, optimal chaining, recursive backbones; `docs/design.md`, `experiments/`).

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
