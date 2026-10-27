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

## Using plegma from other tools

`pixi run build-capi` builds `build/libplegma.so`: the consensus of a read set as one call, from C, C++, Python or
Mojo. Each call takes the reads with where they are anchored (`both`, `left`, `right`, a local assembly's reads
between two anchors) and their strands, and returns what `plegma poa -r 0 [-H model]` writes, byte for byte: one
sequence, or the left and the right side when no read spans both anchors and the two could not be joined. A
homopolymer model (shasta2's csv) is loaded once and shared by every thread; every call is re-entrant. Wrong
arguments and failures come back as messages (C), `std::runtime_error` (C++) or `plegma.Error` (Python), never an
abort. Letters other than A, C, G, T are read as N; an empty read is allowed (adjacent anchors); reads far from the
consensus (over 30% edits) shape the graph but are left out of polishing.

C++ (`capi/plegma.hpp`, header-only over the C API in `capi/plegma.h`):

```cpp
#include "plegma.hpp"
plegma::Model model = plegma::Model::load("HomopolymerModel-ONT-R10.4.1-sup5.0.0-2026-09-28.csv");
std::vector<plegma::Read> reads = {{spanning, plegma::Anchor::Both, plegma::Strand::Forward},
                                   {fixedLeft, plegma::Anchor::Left, plegma::Strand::Reverse}};
std::vector<std::string> parts = plegma::consensus(reads, &model, plegma::Options().pinned(false));
```

With CMake, `find_package(plegma CONFIG REQUIRED PATHS <plegma>/capi)` and `target_link_libraries(tool PRIVATE
plegma::plegma)`; by hand, `-I<plegma>/capi -L<plegma>/build -lplegma -Wl,-rpath,<plegma>/build`. `build-capi` copies
the Mojo runtime libraries the library needs beside it, where it finds them by its own run path (`$ORIGIN`): build/
alone is what a tool links against and ships, with no pixi environment at run time (glibc 2.35 and libstdc++ of GCC
12 or newer).

Python (`python/plegma.py`, ctypes over the same library, no build step of its own):

```python
import plegma  # python/ on the path; the library at $PLEGMA_LIBRARY or ../build/libplegma.so
model = plegma.Model.load("HomopolymerModel-ONT-R10.4.1-sup5.0.0-2026-09-28.csv")
parts = plegma.consensus(reads, anchors=["both", "left"], strands=[0, 1], model=model, pinned=False)
two = plegma.clustered_consensus(reads, model=model, max_clusters=2, min_freq=0.25, pinned=False)
two.consensus, two.members, two.placements[0].cigar  # also plegma.phased_consensus(reads, haplotypes, ...)
```

Mojo: `plegma.consensus(seqs^, anchors, strands, ConsensusOptions(), model)`, the function behind the others.

The options (`plegma_options` in C, set by `plegma_options_init` first; `plegma::Options` in C++; keyword arguments in
Python) are `plegma poa`'s: polishing rounds (`-P`), pinning (`-S`; off is `-N`), `-c`, the scoring (`-M -X -O -E`)
and the band (`-b -f -B`). A local assembler of one place's reads turns pinning off: there the pins only cost time
(shasta2's chr21: the same consensus in 2.4 times less). `plegma_options` starts with its own size, so a tool built
against an older `plegma.h` keeps working with a newer library. `plegma_version()` (C), `plegma::version()` (C++),
`plegma.version()` (Python) and `plegma --version` give the version. `pixi run test-capi` and `pixi run test-python`
test the two (threads sharing a model, wrong arguments, empty and odd reads included).

Haplotypes, for a diploid caller (abPOA's `-d`): `plegma_clustered_consensus` (C++ `plegma::clusteredConsensus`)
gives up to `max_clusters` consensus sequences of one read set, each from the reads nearest it, a cluster under
`min_freq` of the reads joining another; `plegma_phased_consensus` (`plegma::phasedConsensus`) gives the two of reads
already phased. Either returns each cluster's reads and each read's CIGAR to its cluster's consensus; a read on
neither anchor (`PLEGMA_ANCHOR_INSIDE`) builds nothing and is aligned within one. Clusters apart only by homopolymer
lengths are one unless, with a model, two lengths explain the reads `two_lengths_db` better than one (15 dB;
without a model, unless one differs by more than one letter): reads' length errors in runs split one haplotype's
reads far more often than a true difference does. pgphase (github.com/kokyriakidis/pgphase, branch `plegma-msa`,
`PGPHASE_MSA=plegma`) uses them as its default MSA (`--msa abpoa` for abPOA); at pgphase 081cfcc, on HG002 HiFi chr12, chr18 and chr20 against the v1.1 truth (its
noisy regions' het calls; a HiFi homopolymer model learnt on chr20, `two_lengths_db` 5):

| MSA | true het | false het | wrong het | missed het | phased correct | discordant | unphased |
|---|---|---|---|---|---|---|---|
| abPOA | 28,133 | 1,320 | 1,831 | 1,136 | 1,036,220 | 5,173 | 912 |
| plegma | 29,146 | 233 | 1,475 | 712 | 1,036,664 | 4,509 | 1,132 |

pgphase writes plegma's hets unphased where the reads are too noisy to phase them (a long homopolymer HiFi reads one
letter short as often as not); abPOA's need the reads, since written so they add 420 false hets for 150 true.

shasta2's Plegma local assembly method uses the library when configured with `-DPLEGMA_DIR=<plegma checkout>`
(github.com/kokyriakidis/shasta2, branch la7-dump). On its HG002 chr21 assembly (the Final graph reassembled, lassen
against HG002 v1.1):

| consensus | errors | wall | CPU | peak memory |
|---|---|---|---|---|
| abPOA (shasta2's default) | 11,302 | 0:34 | 249 s | 3.7-4.0 GB |
| plegma (library, `-H`, unpinned) | 6,750 | 0:20 | 129 s | 2.3 GB |

(shasta2's msa1 repair: 7,353 errors.) Each timed twice, 10 threads, on an otherwise idle machine.

Per local assembly, in one process, on 2,093 of them (every 17th that reaches a consensus method), routed as shasta2's
default routes them (abPOA with 6 or more spanning reads up to 5 kb, theseus otherwise), plegma on every read with
the model:

| shasta2's default sends it to | local assemblies | shasta2's default | plegma |
|---|---|---|---|
| abPOA | 1,638 | 8.10 s | 3.53 s |
| theseus, spanning reads (over 5 kb) | 22 | 2.42 s | 2.19 s |
| theseus, every read (under 6 spanning) | 433 | 0.161 s | 0.135 s |
| total | 2,093 | 10.69 s | 5.86 s |

The 1,852 under 1 kb (some 16 reads of 200 letters) take 1.80 s by shasta2's default and 1.74 s by plegma; peak
memory over all 2,093 is 288 MB against 152 MB. On abPOA's own input for all 2,093 (spanning reads alone), abPOA
takes 30.7 s and a 6.0 GB peak.

`plegma poa`'s band is abPOA's adaptive band, `b + f L` letters either side of the row's best cells and of the
graph's remaining path, capped at 64 letters (`-B`; a negative value leaves it uncapped, as abPOA). The cap leaves
shasta2's chr21 assembly byte for byte the same and its long local assemblies take 2.6 times fewer cells.

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
