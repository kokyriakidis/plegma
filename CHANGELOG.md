# Changelog

## Unreleased

### Fixed

- `plegma poa` on a gzipped FASTA crashed: the header words (`anchor=`, `strand=`) were read again from the raw file
  as text. The input is now read once, gzip or not, FASTA or FASTQ, CRLF or LF.
- Out-of-range option values (`-r 7`, `-d 0`, `-q 2`, `-P -1`, `-t 0`, ...) were accepted; `-r 7` wrote nothing and
  exited 0. Malformed numbers, unknown options and missing files ended in an unhandled exception. Errors are now one
  `plegma: ...` line on stderr and exit status 1.
- `-i` with MSA rows of different lengths read past a row's end (garbage letters, exit 0) or cut a row short; it now
  fails.
- The consensus between spanning reads' ends weighed edges in 32 bits: past some 46,000 reads it overflowed and the
  consensus collapsed to a letter.
- Under `-S` (the default) a right-anchored read could end on letters a left read carried past the right anchor. With
  no spanning read and the two sides not joined, right reads were fused onto the left side's path where the two
  overlapped, splitting the MSA.
- The banded 16-bit kernel ignored the scores' size: at every score a hundredfold the rows wrapped and the MSA
  changed. Default scores were never affected.
- `-s` with `-H`: a reverse-complemented read kept its header's strand for the homopolymer model. The library's
  `consensus` with `PoaOptions.strand` had the same fault.
- Polishing: nine or more reads whose band could not reach both ends of a repeat unit wrapped a unit deletion's
  gain to the top of 32 bits; the round was reverted and polishing stopped early.
- Homopolymer model csv: run lengths past 65,535 were read into another row, integers past 63 bits were written
  outside the table, and a trailing comma or an exponent without digits was accepted. These now fail.
- `phased_consensus`: with no labelled read both haplotypes came out empty (every read placed as an insertion); with
  every read inside, both were empty.
- C API: a `count` of `SIZE_MAX` (or a length of 2^62) aborted the caller's process in the allocator; it now fails
  with a message. `plegma_options_init(NULL)` crashed.
- Python: integers a C `int32_t` cannot hold wrapped silently (`anchors=[2**32]` ran as both).
- `poa` with no reads and no starting graph read past an empty list.
- 16-bit kernel: when two predecessors' rows had bases more than 32,767 apart (a source chain joining a long path far
  down, as a right read past the left anchor makes), the difference wrapped and the lower row's scores came in far
  too high; a 20 kb read then scored 67,510 where 40,000 was the most possible. A mismatch penalty past some 16,000
  also wrapped the 16-bit lanes. Both now take full width or 32-bit lanes.
- With no spanning read, the two sides were never joined when they overlapped wholly, or when one right read
  reached past the left anchor: the join refused every pin with the right side starting before the left or ending
  where it does, and two unpolished sequences came out. HG00512's HLA-DRB1 now joins at a pin with 7 reads either
  side, not 2.
- `ref -F` wrote an indel where the window between two pins put it, not where minimap2 writes it: minimap2 shifts
  each gap left through the matched run before it last (mm_fix_cigar), across pins. 132 of 141 real chr20 allele
  pairs now have minimap2's CIGAR exactly (102 before); the other 9 are equal-cost ties. Past 2^25 letters `-F`
  silently left the pins out; it now fails.
- Polishing tried whole repeat units of up to 24 letters, but a read's band (16 letters either side) reaches a unit
  away only up to 16: units of 17 to 24 letters were scored wrongly and, with noisy reads, rarely taken. The band is
  24 now (simulated arrays of such units a unit or three off: 142 of 192 exact before, 191 now; 0 to 14% more time).
  The likelihood stage let one read whose band could not reach both ends of a unit veto deleting it (log 0); such a
  read now counts the unit's letters as edits.
- `phased_consensus`: with the two haplotypes' first consensus sequences of one length, swapping the labels 0 and 1
  changed the result, not only its order.
- `plegma poa -d` and `-a` were ignored when any read was anchored at one end; they now fail with a message.
- Joins inside tandem arrays: one 19-mer pin, unique in each side, could be another copy's, and the join dropped or
  doubled copies. Each way to place the two sides is now scored by how well they agree over the whole overlap; where
  two placements on different diagonals agree nearly as well, the sides are left apart. Inside a tandem array
  (its period found from where its k-mers recur, so copies with indels count) every copy number is a candidate:
  the join's meeting point moved by whole copies to the same place in another copy, aligned out from there with
  the far ends free, and compared with the others on the letters both align. Simulated local assemblies, 504
  held-out cases: wrong joins 38 to 0 (main grid 8 to 0). Where no read crosses an array and its copies differ by
  fewer letters than the consensus's own errors the copy number cannot be told, and the sides stay apart: 29 of
  the joins right before were right on two or three letters of evidence and are refused now.
- A one-sided assembly's cut (`-c`) and the join's read support took a read's length for its reach along the
  consensus; they now take the consensus columns its letters are in. `-c 1` skipped the cut.
- A one-sided read whose alignment lost its path in the band (at 10% errors) was never aligned again, as an end is
  always reached; under a third of a perfect score it is now aligned once more with the band doubled.
- `clustered_consensus`: a partial read that tells the clusters nothing went to the first input row's cluster and
  counted there for `min_freq` and the 30% one-unit rule, so input order decided whether two haplotypes stayed
  apart. Such reads now count for no cluster. A read over the letters where the two consensus sequences differ at
  one cost against both lies between them and counts for its own: counted for none, an A21 haplotype whose
  consensus a homopolymer model made A22 fell under `min_freq` against an A20 one (pgphaseMojo, HG002 chr20: 13 het
  run-length indels read homozygous or lost).
- A read given as spanning whose end disagrees with the reads anchored at that end (anchored by counting letters to
  the anchor where the haplotype's length is not the reference's) set where the consensus ended: HG03009's
  HLA-DQA1 haplotype 1 ended 5.3 kb before its right reads' anchor, HG00512's and NA19003's HLA-DRB1 began 2 kb
  before their left reads'. Each spanning read's first and last 100 letters are now checked against those of the
  reads anchored there; at an end no other spanning read shares, where most disagree, it is anchored at its other
  end alone, or, disagreeing at both, builds nothing (another copy's). Only at an end one-sided reads are anchored
  at: reads all spanning are as given. (Counted by majority alone, a haplotype without a gene the other has lost its
  12 spanning reads to 20 right reads of the other: HG00512's and NA20509's second DRB slot read the gene present.)
- A one-sided assembly began where the most reads began node by node: reads whose anchored ends differ by a few
  letters spread their votes, and 8 short reads sharing another start outvoted 11 (HG002's HLA-DRB1 under GRCh38,
  151 edits in the first 300 letters). The start, and a right side's end, are now voted over 16 nodes.
- One-sided reads running past their own anchor (anchored by counting letters to it, not by alignment): a later
  read's pinned alignment onto the first such read's letters there was thrown away for an unpinned one, whose band
  lost the read's path, and the consensus collapsed (panStar, NA19007's HLA-DRB1 haplotype 2: 2,152 edits from
  its allele; HG03009's HLA-DQA1 haplotype 1 4.6 kb short). Its letters past the anchor are now inserted there and
  the rest keep their pinned places: NA19007's consensus is the reads' sequence between the anchors, edit for edit;
  HG03009's is again what 0.4.0 made, right for 14 kb and then off where its one read anchored at both ends stops
  5.3 kb short of the right anchor the right reads reach.
- The consensus depended on input order (abPOA's guide order breaks ties by position, and takes the earlier read
  of the first pair first): ties now go by the sequences, and the same reads in any order give one consensus.
- The likelihood stage judged its run-slip gate on as little as one inserted letter, and on clean reads never ran.
- `align_pinned`'s score left out inserted letters between anchors and path nodes skipped between them (`-s`
  compares these scores).
- `ref -F`: the runs inserted at one place were aligned to each other with the fallback costs the pinned mode
  otherwise refuses.

### Changed

- Options are checked in one place (`PoaOptions.check`, `ConsensusOptions.check`) by the command, the C API and
  `consensus`. The command now refuses `-c 0` and `-E 0`, as the C API did.
- `plegma poa -r 0` and the library's `consensus` share one code path (`consensus.polished`).
- `poa` records the reads it reverse-complemented (`Poa.reversed`).

### Added

- `PLEGMA_TRACE=poa,pins,polish,clustered,phased,reference` (or `all`): the decisions behind a result, one line each
  on stderr. It replaces `PLEGMA_PINS_TRACE`.
- `clustered_consensus` with no spanning read and one side only (a telomere, or a window whose other end no read
  reaches) splits that side's reads into haplotypes by their own MSA; it gave one consensus before. A column past a
  read's far end counts as no letter (`cluster_rows`' `MISSING`), and clusters are weighed over the letters both
  reach. HG002's chr18 q-end, 12 left reads: one consensus 62 and 69 edits from the two assembled haplotypes'
  ends (hg002v1.1), now 4 and 8 (pgphaseMojo's 11-read call: 0 and 6). Simulated reads (8-30 of one or two
  haplotypes, letters apart near the anchor): at 1% errors 37 of 300 one-haplotype sets split (the spanning path: 50)
  and 292 of 300 two-haplotype sets split purely (279); at 3%, 49 (25) and 264 (190).
- `Clustered.joined` (`clustered_consensus`, `phased_consensus`): per cluster, False where its reads made two sides
  that could not be joined (no read spanning). The consensus is still the longer side, as before, but a caller can
  now tell it from one reaching both anchors. Mojo only; the C API and Python do not expose it yet.
- `Placement.margin`: how much more the next cheapest cluster costs a read placed by cost (a `phased_consensus` read
  with no haplotype, an inside read of `clustered_consensus`); 0 for reads placed by their cluster, and where
  two clusters tie (the read then goes to the first). Mojo only, as `joined`.
- `pixi run test-all`: unit tests, C API (shared and static), Python and fuzzing. `experiments/validate/golden.py`
  checks a refactor byte for byte.

### Internal

- `clustered.mojo` is split into `placement.mojo`, `runs.mojo`, `phased.mojo` and `clustered.mojo`, and
  `polish.mojo` into `bands.mojo`, `likelihood.mojo`, `run_lengths.mojo` and `polish.mojo`. The command line is split
  into `cli/options.mojo`, `cli/records.mojo`, `cli/gfa.mojo` and `cli/plegma_cli.mojo`. `plegma.clustered` and
  `plegma.polish` still export the names callers import from them.

## 0.4.0 (2026-10-27)

See the git history before this file.
