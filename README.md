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
  abPOA 357; 58 on simulated sets, abPOA 458). By default the sequences join in abPOA's guide order (`-p`) and
  each is pinned to the consensus so far as minimap2 pins an alignment (`-S`: its seeds and chain, pins only on
  k-mers unique in both), the DP running only between the pins: the same consensus as aligning each in input order
  to the whole graph (`-n -N`), better allele MSAs (pgrealign qc pair F1 0.872 against 0.864 on 377 chr20 regions),
  up to 7 times faster and in a sixth of the memory on the largest chr20 VNTRs.
  Reads of a local assembly may be fixed at one anchor only (`anchor=left` or `anchor=right` in their FASTA
  headers, `plegma/poa.mojo`): aligned from the graph's start to anywhere (or from anywhere to its end), whatever
  passes the graph's end growing it; with no read spanning both anchors each side is assembled and the two joined
  where they overlap, or one side alone stops where fewer than `-c` reads (5) go on (a telomere). The consensus is
  polished against every read as it is anchored. On simulated local assemblies (`experiments/msa/anchored.py`) it
  reproduces 5 kb truths exactly in 0.2 s and 20 MB, where pericles (theseus-lib's consensus branch), which needs a
  spanning read, takes 17 to 63 s and its majority vote leaves 4 to 29 errors. A one-sided read much further from the
  spanning reads' consensus than they are (another copy's, sharing the one anchor) gets no vote.
  With `-H` (a homopolymer model in shasta2's csv, `plegma/homopolymer.mojo`) every run then takes its most likely
  length given the reads' lengths there, each read on the strand its header gives (`strand=0|1`), as shasta2's msa1
  repair chooses it. In shasta2's local assemblies of HG002 chr21 (its Final graph reassembled with plegma, scored by
  lassen against HG002 v1.1) it leaves 6749 errors, shasta2's msa1 repair 7353 and its abPOA default 11302.
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

## Library API (C and C++)

`pixi run build-capi` builds `build/libplegma.so`, the consensus of a read set as a call: `capi/plegma.h` (C) and
`capi/plegma.hpp` (C++, header-only over it). A homopolymer model is loaded once (from a file, or from its text in
memory) and shared by every thread; each call takes the reads with their anchors and strands and returns what
`plegma poa -r 0 [-H model]` writes, byte for byte (checked on 2,000 of shasta2's local assemblies). Errors come back
as messages (C) or `std::runtime_error` (C++), never an abort.

```cpp
#include "plegma.hpp"
plegma::Model model = plegma::Model::load("HomopolymerModel-ONT-R10.4.1-sup5.0.0-2026-09-28.csv");
std::vector<plegma::Read> reads = {{spanning, plegma::Anchor::Both, plegma::Strand::Forward},
                                   {fixedLeft, plegma::Anchor::Left, plegma::Strand::Reverse}};
std::vector<std::string> parts = plegma::consensus(reads, &model, plegma::Options().pinned(false));
```

`plegma_options` (C) or `plegma::Options` (C++) carries `plegma poa`'s options: polishing rounds (`-P`), pinning
(`-S`; off is `-N`) and `-c`. A local assembler of one place's reads turns pinning off: there the pins only cost
time (shasta2's chr21: the same consensus in 2.4 times less).

Link with `-Lbuild -lplegma -Wl,-rpath,<plegma>/build`; the library finds the Mojo runtime it needs (in plegma's
pixi environment) by its own run path. `pixi run test-capi` tests it (threads sharing a model included).
shasta2's Plegma local assembly method uses it when configured with `-DPLEGMA_DIR=<plegma checkout>`
(github.com/kokyriakidis/shasta2, branch la7-dump). On its HG002 chr21 assembly (the Final graph reassembled, lassen
against HG002 v1.1):

| consensus | errors | wall | CPU | peak memory |
|---|---|---|---|---|
| abPOA (shasta2's default) | 11,302 | 0:34 | 249 s | 3.7-4.0 GB |
| plegma (library, `-H`, unpinned) | 6,750 | 0:31 | 223 s | 3.2 GB |

(shasta2's msa1 repair: 7,353 errors.) Both timed twice in the same state of the machine, 10 threads.

Per local assembly, in one process, on 2,093 of them: abPOA (spanning reads, as shasta2 runs it) 30.9 s and a 6.0 GB
peak; plegma on the same reads 5.9 s and 131 MB (under 1 kb, level: 1.96 against 1.97 s); plegma on every read
with the model 11.2 s and 158 MB.

The library (`plegma/__init__.mojo` lists its API) holds both: the graph and its whole-graph aligner
(`plegma/kernel.mojo`), consensus polishing (optionally with a `HomopolymerModel`), and `on_reference` with minimap2's
pinned alignment.

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
