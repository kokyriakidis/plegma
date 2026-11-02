# plegma

Fast, accurate partial order alignment (POA) of DNA: consensus sequences, multiple sequence alignments and POA
graphs, as a command line tool and as a library for C, C++, Python and Mojo.

*πλέγμα* (plegma) is Greek for a woven mesh or lattice, which is what a POA graph is: sequences threaded into one
another where they agree.

plegma is built for the consensus problems inside genome tools:

- **Read consensus** (`plegma poa`): abPOA's flags and output, followed by polishing that fixes most of the errors a
  POA consensus makes in tandem repeats and homopolymer runs, and optionally a homopolymer run-length model.
- **Local assembly**: reads fixed at one anchor or both, as an assembler collects them between two anchors, with
  each read's strand fed to the run-length model.
- **Haplotypes**: several consensus sequences of one read set (clustered, or two from reads already phased), with
  each read's CIGAR to its haplotype, for a diploid variant caller.
- **Reference-anchored MSA** (`plegma ref`): a pangenome region's alleles aligned to its reference allele the way
  minimap2 and dipcall write differences, so the MSA's variants match a truth set's representation.

## Results

| Task | Data, truth | plegma | Alternatives |
|---|---|---|---|
| shasta2 local assembly consensus | HG002 ONT R10.4.1 chr21 assembly, errors against HG002 v1.1 (lassen) | **6,755** | msa1 7,353; abPOA (default) 11,302 |
| pgphase het calls in noisy regions | HG002 HiFi chr12, 18, 20 against v1.1: false / wrong / missed hets | **190 / 679 / 527** | abPOA 1,345 / 2,134 / 1,130 |
| Read consensus | 180 read sets simulated from real VNTR haplotypes, errors against truth | **8** | abPOA 357 |
| Pangenome region MSA (`ref -F`) | HPRC graph chr20 realigned by pgrealign, HG002 structural F1 | **0.737** | pgrealign's star 0.614 |

Details, versions and how each number was made are under [Benchmarks](#benchmarks).

## Building

plegma builds with [pixi](https://pixi.sh), which fetches the pinned Mojo nightly. It is developed and tested on
Linux x86-64. Its pairwise alignment comes from [karan](https://github.com/kokyriakidis/karan) and its FASTA/FASTQ
reading from [sporos](https://github.com/kokyriakidis/sporos), both git submodules.

```sh
git clone --recurse-submodules https://github.com/kokyriakidis/plegma
cd plegma
pixi run build-cli    # the command line tool: build/plegma
pixi run build-capi   # the shared library: build/libplegma.so (with the Mojo runtime beside it)
```

## Command line

```
plegma <poa|ref> [options] <in.fa/fq> > out
```

Input is one FASTA or FASTQ file (gzip is fine). `-r` picks the output: `0` consensus FASTA (default), `1` MSA FASTA,
`2` both, `3` GFA, `4` GFA with the consensus path. `plegma -h` lists every option.

### `plegma poa`

```sh
plegma poa reads.fa > consensus.fa                      # polished consensus
plegma poa -r 1 reads.fa > msa.fa                       # the MSA
plegma poa -d 2 reads.fa > haplotypes.fa                # up to two consensus sequences
plegma poa -N -H HomopolymerModel.csv la.fa > cons.fa   # a local assembly, with a run-length model
```

The flags, defaults and alignment modes (`-m` 0 global, 1 local, 2 extension) are abPOA's.

| Option | Default | Meaning |
|---|---|---|
| `-M -X -O -E` | `2 4 4,24 2,1` | match, mismatch, gap open and extend (two pieces: a gap of g costs min(O1 + g E1, O2 + g E2)) |
| `-b -f -B` | `10 0.01 64` | adaptive band: b + f L letters either side of the best cells, capped at B (`-B -1`: no cap, as abPOA) |
| `-p` / `-L` / `-n` | `-p` | order: abPOA's guide order (minimizer similarity), longest first, or input order |
| `-S` / `-N` | `-S` | each sequence pinned to the consensus so far at minimap2's chained seeds (DP only between pins), or aligned to the whole graph |
| `-P` | `50,20` | polishing rounds: by edit distance, then by pair-HMM likelihood (`0,0` off) |
| `-H FILE` | | homopolymer model (shasta2's csv): every run takes its most likely length given the reads' |
| `-d` `-q` | `1 0.25` | at most d consensus sequences; the smallest share of reads a cluster needs |
| `-c` | `5` | a one-sided assembly stops where fewer reads go on |
| `-s` | | try reverse complements of poorly scoring sequences |
| `-i FILE` | | start from an existing MSA |
| `-l` | | the input is a list of files, each aligned on its own |

Pinning (`-S`) pays on long, divergent sequences such as a pangenome region's alleles; on reads of one place it only
costs time, so local assembly uses `-N`.

**Local assembly reads.** Words in a read's FASTA header say where it is fixed and on which strand it was read:

```
>read1 anchor=both strand=0     spans both anchors (the default)
>read2 anchor=left strand=1     starts at the left anchor, ends anywhere; reverse-complemented to the consensus
>read3 anchor=right             ends at the right anchor
```

Reads fixed at one anchor are aligned from that end to anywhere and may grow the graph past the other end. When no
read spans both anchors, each side is assembled and the two joined where they overlap; a side alone stops where
fewer than `-c` reads go on (a telomere). A one-sided read much further from the spanning reads' consensus than they
are (another repeat copy sharing the anchor) gets no vote. The strand matters only with `-H`.

### `plegma ref`

```sh
plegma ref -F flanks.fa -t 8 -r 1 alleles.fa > msa.fa
```

The first sequence is the reference. Every other is aligned to it alone at minimap2's asm5 costs, gaps
left-normalized, its letters placed in the reference's columns, and insertions at the same place aligned to each
other, so a copy-number change stays one event. With `-F` (a FASTA of the reference's left and right flanks) each
alignment is pinned where minimap2 would pin it at dipcall's settings and is exact between the pins, so differences
are written as dipcall-made truth sets write them. `-t` sets the number of threads.

## Library

`pixi run build-capi` builds `build/libplegma.so`. One call takes a read set (each read with its anchor and strand)
and returns what `plegma poa -r 0 [-H model]` writes, byte for byte: one sequence, or the left and right sides when
no read spans both anchors and they could not be joined.

- **Threads**: every call is re-entrant; a homopolymer model is loaded once and shared read-only by every thread.
- **Errors**: wrong arguments and failures come back as messages (C), `std::runtime_error` (C++) or `plegma.Error`
  (Python), never an abort.
- **Input**: letters other than A, C, G, T are read as N; empty reads are allowed (adjacent anchors); reads with
  more than 30% edits to the consensus shape the graph but are left out of polishing.
- **Compatibility**: `plegma_options` starts with its own size, so a tool built against an older `plegma.h` keeps
  working with a newer library. `plegma_version()`, `plegma::version()`, `plegma.version()` and `plegma --version`
  give the version.

### C++

`capi/plegma.hpp` is header-only over the C API in `capi/plegma.h`.

```cpp
#include "plegma.hpp"

plegma::Model model = plegma::Model::load("HomopolymerModel-ONT-R10.4.1-sup5.0.0-2026-09-28.csv");
std::vector<plegma::Read> reads = {{spanning, plegma::Anchor::Both, plegma::Strand::Forward},
                                   {fixedLeft, plegma::Anchor::Left, plegma::Strand::Reverse}};
std::vector<std::string> parts = plegma::consensus(reads, &model, plegma::Options().pinned(false));
```

With CMake:

```cmake
find_package(plegma CONFIG REQUIRED PATHS <plegma>/capi)
target_link_libraries(mytool PRIVATE plegma::plegma)
```

By hand: `-I<plegma>/capi -L<plegma>/build -lplegma -Wl,-rpath,<plegma>/build`. `build-capi` copies the Mojo runtime
libraries beside `libplegma.so`, where it finds them by its run path (`$ORIGIN`), so `build/` alone is what a tool
links against and ships: no pixi environment at run time, glibc 2.35 and libstdc++ of GCC 12 or newer.

### Python

`python/plegma.py` is ctypes over the same library; it has no build step of its own.

```python
import plegma  # python/ on the path; the library at $PLEGMA_LIBRARY or ../build/libplegma.so

model = plegma.Model.load("HomopolymerModel-ONT-R10.4.1-sup5.0.0-2026-09-28.csv")
parts = plegma.consensus(reads, anchors=["both", "left"], strands=[0, 1], model=model, pinned=False)

two = plegma.clustered_consensus(reads, model=model, max_clusters=2, min_freq=0.25, pinned=False)
two.consensus, two.members, two.placements[0].cigar
```

### Mojo

`plegma.consensus(seqs^, anchors, strands, ConsensusOptions(), model)` is the function behind the others.
`plegma/__init__.mojo` lists the package's public API: the graph and its aligner, pinned alignment, polishing, the
homopolymer model, clustering and `on_reference`.

### Options

`plegma_options` (C, filled with defaults by `plegma_options_init`), `plegma::Options` (C++) and keyword arguments
(Python) carry `plegma poa`'s options:

| Field | CLI | Default |
|---|---|---|
| `polish_rounds`, `likely_rounds` | `-P` | 50, 20 |
| `pinned` | `-S` / `-N` | 1 (a local assembler sets 0) |
| `min_cover` | `-c` | 5 |
| `match`, `mismatch`, `gap_open1`, `gap_extend1`, `gap_open2`, `gap_extend2` | `-M -X -O -E` | 2, 4, 4, 2, 24, 1 |
| `band_b`, `band_f`, `band_max` | `-b -f -B` | 10, 0.01, 64 |
| `noise_run` | | 3 (clustered: clusters apart only by single letters in runs this long are one) |
| `two_lengths_db` | | 15 (clustered and phased, with a model: how much better two run lengths must explain the reads) |

### Haplotypes

For a diploid caller (abPOA's `-d`):

- `plegma_clustered_consensus` (C++ `plegma::clusteredConsensus`, Python `plegma.clustered_consensus`) gives up to
  `max_clusters` consensus sequences of one read set, each from the reads nearest it; a cluster under `min_freq` of
  the reads joins another.
- `plegma_phased_consensus` (`plegma::phasedConsensus`, `plegma.phased_consensus`) gives the two consensus sequences
  of reads already phased.

Both return each cluster's reads and each read's CIGAR to its cluster's consensus. A read on neither anchor
(`PLEGMA_ANCHOR_INSIDE`) builds nothing and is only aligned within a cluster. Clusters apart only by homopolymer
lengths are merged unless, with a model, two lengths explain the reads `two_lengths_db` better than one (without a
model, unless one differs by more than one letter): read length errors in runs split one haplotype's reads far more
often than a true difference does. Likewise, clusters apart by one short tandem repeat unit need 30% of the reads.

## How it works

- **Alignment**: each sequence is aligned to the graph under abPOA's scoring and adaptive band, capped at 64
  letters either side; memory follows the band, not the graph times the sequence. By default sequences join in
  abPOA's guide order and are pinned to the consensus so far as minimap2 pins an alignment (its seeds and chain, pins
  only on k-mers unique in both), with the DP only between pins.
- **Polishing**: one-letter and whole-repeat-unit edits to the consensus are kept while they lower the reads' summed
  edit distance to it, then while they raise the reads' likelihood under a pair HMM.
- **Homopolymer model**: with shasta2's empirical run-length model, every run of two or more letters takes its most
  likely length given the reads' lengths there, each read on its own strand, as shasta2's msa1 repair does. A read
  votes on a run only when it has the letters on both sides of it.
- **Reference MSA**: `on_reference` aligns each allele exactly to the reference (or between minimap2's pins, from a
  port of its asm5 seeding and chaining, `plegma/pins.mojo`).

## Benchmarks

### shasta2

shasta2's `Plegma` local assembly method calls the library in-process when configured with
`-DPLEGMA_DIR=<plegma checkout>` (github.com/kokyriakidis/shasta2, branch `la7-dump`). The Final assembly graph of
an HG002 ONT R10.4.1 chr21 run was reassembled with each consensus method, then scored by
[lassen](https://github.com/paoloshasta/lassen) against HG002 v1.1 (65 mappings, 83.3 Mb aligned). plegma at
0737184, with the R10.4.1 homopolymer model, unpinned:

| Consensus | Errors | Mismatches | Insertions | Deletions | vs default |
|---|---|---|---|---|---|
| abPOA (shasta2's default) | 11,302 | 76 | 1,891 | 9,335 | |
| shasta2's msa1 repair | 7,353 | 52 | 3,488 | 3,813 | −34.9% |
| **plegma** | **6,755** | **49** | **3,213** | **3,493** | **−40.2%** |

Over the 25 longest segments: abPOA 11,046, msa1 7,194, plegma 6,588 (better than the default in 21 segments, worse
in 1). Most of what remains is runs of 10 or more letters at the reads' information limit.

Speed (measured at 005e225, each twice, 10 threads, an otherwise idle machine): the whole reassembly takes 0:20
wall, 129 s CPU and 2.3 GB with plegma, 0:34, 249 s and 3.7–4.0 GB with abPOA. Per local assembly, in one process,
on 2,093 of them (every 17th that reaches a consensus method), each routed as shasta2's default routes it, with
plegma given every read and the model:

| shasta2's default sends it to | Local assemblies | shasta2's default | plegma |
|---|---|---|---|
| abPOA (6+ spanning reads, up to 5 kb) | 1,638 | 8.10 s | 3.53 s |
| theseus, spanning reads (over 5 kb) | 22 | 2.42 s | 2.19 s |
| theseus, every read (under 6 spanning) | 433 | 0.161 s | 0.135 s |
| **Total** | **2,093** | **10.69 s** | **5.86 s** |

Peak memory over all 2,093 is 288 MB for the default against 152 MB for plegma. On abPOA's own input for all of them
(spanning reads alone), abPOA takes 30.7 s and a 6.0 GB peak. The band cap (`-B 64`) leaves this assembly byte for
byte the same and cuts the cells of its long local assemblies 2.6-fold.

### pgphase

[pgphase](https://github.com/kokyriakidis/pgphase) uses the haplotype calls as its default MSA (`--msa abpoa` for
abPOA). At pgphase 6a03ef2 (plegma be505b6), on HG002 HiFi chr12, chr18 and chr20 against the v1.1 truth, both MSAs
scored on the het calls of the same noisy regions, with a HiFi homopolymer model learnt on chr20 (`two_lengths_db` 5
on reads already phased, 8 when clustering them all):

| MSA | True het | False het | Wrong het | Missed het | Phased correct | Discordant | Unphased |
|---|---|---|---|---|---|---|---|
| abPOA | 33,026 | 1,345 | 2,134 | 1,130 | 1,036,198 | 5,194 | 913 |
| **plegma** | **35,038** | **190** | **679** | **527** | **1,036,841** | **4,416** | 1,048 |

pgphase writes plegma's hets unphased where the reads are too noisy to phase them (HiFi reads a long homopolymer one
letter short as often as not); abPOA's would need the reads, since writing them adds 420 false hets for 150 true.

### Read consensus

`experiments/msa/accuracy.py` at 0737184, default options: each consensus scored by its edit distance to the
truth. The VNTR sets are 10 or 30 full-length reads of a real VNTR haplotype from vntr-eval
(`experiments/msa/vntr_reads.py`) at 1%, 5% or 10% errors, substitutions, insertions and deletions alike. The
simulated sets are 10, 50 or 200 reads of random 0.5, 2 or 10 kb sequences, with uniform errors (1–15%) or
homopolymer-biased indels (5% or 10%).

| Read sets | plegma errors | plegma exact | abPOA errors | abPOA exact |
|---|---|---|---|---|
| 180 VNTR read sets | **8** | **174** | 357 | 122 |
| 54 simulated sets | **60** | **43** | 458 | 27 |

On simulated local assemblies (`experiments/msa/anchored.py`) plegma reproduces 5 kb truths exactly in 0.2 s and
20 MB, where pericles (theseus-lib's consensus branch), which needs a spanning read, takes 17 to 63 s and leaves 4 to
29 errors.

### Pangenome region MSA

Spliced into chr20 of the HPRC graph by [pgrealign](https://github.com/benedictpaten/pgrealign)
(`experiments/chr20`), `ref -F`'s graphs genotype HG002 better than pgrealign's own star MSA:

| Structural F1 | plegma `ref -F` | pgrealign star |
|---|---|---|
| chr20 | **0.737** | 0.614 |
| chr20 repeats | **0.703** | 0.556 |
| vntr-eval hard set | **0.452** | 0.347 |
| chr20, truvari refine | **0.758** | 0.667 |
| hard set, truvari refine | **0.476** | 0.385 |

truvari refine scores haplotypes however their records are split. On the largest chr20 VNTRs `ref` is 17–19 times
faster than abPOA in a tenth to a hundredth of its memory. In `poa`, the default pinning gives the same consensus as
aligning in input order to the whole graph (`-n -N`), better allele MSAs (pgrealign qc pair F1 0.872 against 0.864
on 377 chr20 regions), and is up to 7 times faster in a sixth of the memory on the largest chr20 VNTRs.

## Development

```sh
pixi run test         # unit tests
pixi run test-capi    # C/C++ API: threads sharing a model, wrong arguments, empty and odd reads
pixi run test-python  # Python bindings
pixi run fuzz         # random inputs through every option combination, every output checked
pixi run format       # mojo format, 120 columns
```

| Path | Contents |
|---|---|
| `plegma/` | the library (Mojo); `__init__.mojo` lists the public API |
| `cli/` | the `plegma` command |
| `capi/` | C API (`plegma.h`), C++ wrapper (`plegma.hpp`), CMake package, the library's tests |
| `python/` | Python bindings and their tests |
| `tests/` | unit tests |
| `experiments/` | benchmarks and evaluation harnesses (`msa`, `chr20`, `validate`) |
| `docs/` | design notes and research reports |
| `deps/` | karan and sporos (submodules) |

## License

MIT
