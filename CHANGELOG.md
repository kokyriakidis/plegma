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
