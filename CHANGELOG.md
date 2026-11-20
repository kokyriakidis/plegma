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
  every copy number is a candidate (the best placement shifted by whole periods). Simulated assemblies (held-out
  seeds): wrong joins 38 to 10; where no read crosses the array and its copies differ by fewer letters than the
  consensus's own errors, the sides now stay apart (a few joins right only by chance before are refused too).
- A one-sided assembly's cut (`-c`) and the join's read support took a read's length for its reach along the
  consensus; they now take the consensus columns its letters are in. `-c 1` skipped the cut.
- A one-sided read whose alignment lost its path in the band (at 10% errors) was never aligned again, as an end is
  always reached; under a third of a perfect score it is now aligned once more with the band doubled.
- `clustered_consensus`: a partial read that tells the clusters nothing went to the first input row's cluster and
  counted there for `min_freq` and the 30% one-unit rule, so input order decided whether two haplotypes stayed
  apart. Such reads now count for no cluster.
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
- `pixi run test-all`: unit tests, C API (shared and static), Python and fuzzing. `experiments/validate/golden.py`
  checks a refactor byte for byte.

### Internal

- `clustered.mojo` is split into `placement.mojo`, `runs.mojo`, `phased.mojo` and `clustered.mojo`, and
  `polish.mojo` into `bands.mojo`, `likelihood.mojo`, `run_lengths.mojo` and `polish.mojo`. The command line is split
  into `cli/options.mojo`, `cli/records.mojo`, `cli/gfa.mojo` and `cli/plegma_cli.mojo`. `plegma.clustered` and
  `plegma.polish` still export the names callers import from them.

## 0.4.0 (2026-10-27)

See the git history before this file.
