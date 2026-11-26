# plegma's design

Two commands share one library. This file maps the code to them, for whoever changes it next.

## `plegma poa`: the general POA

The command line is `cli/`: `options.mojo` (the usage text, parsing, every value checked, the library's own checks
through `ConsensusOptions.check`), `records.mojo` (the input read once, gzip or not, FASTA or FASTQ, with each record's
`anchor=` and `strand=` header words), `gfa.mojo` (`-r 3`, `-r 4`) and `plegma_cli.mojo` (each command, the outputs;
an error is one `plegma: ...` line and exit status 1).

`plegma.poa` (`plegma/poa.mojo`) adds the sequences one at a time, in abPOA's guide order with ties (and the first
pair's order) by the sequences, so any input order gives one consensus (`-p`, the default,
`plegma/order.mojo`), input order (`-n`) or longest first (`-L`).

- **Alignment to the graph** (`plegma/kernel.mojo`, `plegma/kernel_simd.mojo`, `plegma/topology.mojo`): banded DP over
  the nodes in topological order (`Topology`, reusable in-place via `compute` and `compute_allowed`), abPOA's scoring
  and adaptive band, with SIMD lane helpers separated into `kernel_simd.mojo`. Memory follows the band: a row lives
  until its last successor is done, and the traceback keeps one byte per cell (three more at nodes with several
  predecessors). `Workspace.align` aligns between any two nodes (`SOURCE` and `SINK` for the ends), so a sequence can
  also be aligned window by window between anchors.
- **Seeded** (`-S`, the default; `-N` aligns to the whole graph; `plegma/pinned.mojo` `align_pinned`): minimap2's
  pins (`pins.pins`) between the sequence and the consensus path, kept only where the 19-mer centred on a pin is
  unique in both (in near-identical tandem copies the chain may pair the wrong copy), each pin a one-letter anchor on
  its path node; the kernel then aligns each stretch between anchors to the subgraph between their nodes
  (`window_nodes`), in parallel across windows when `threads > 1`, in the band the whole sequence would get. Where no
  pin survives, the one window is the whole graph.
- **The graph & CIGAR** (`plegma/graph.mojo`, `plegma/cigar.mojo`): flat arrays of nodes, columns and edges. A
  sequence joins along its aligned nodes (`fuse_nodes`) or copies an identical read directly (`add_twin_path`);
  `node_order_into` reuses Kahn's scratch buffers; `cigar.mojo` provides zero-allocation byte-level CIGAR walkers
  (`pairing_from_cigar`, `query_at`, `haplotype_at`) alongside `cigar_runs`.
- **Anchored QC & two-side joining** (`plegma/anchored.mojo`, `plegma/join.mojo`): separated from `poa.mojo` so
  spanning-read end verification (`ends_checked`), paralog/foreign-read filtering (`find_foreign`), one-directional
  reach/trimming (`read_reach`, `trim_consensus`), repeat period detection (`repeat_period`), and two-side overlap /
  tandem copy-shift joining (`join_sides`) are independently testable and parallelized across reads and candidates.
- **Consensus & reusable scratch** (`plegma/consensus.mojo` `ConsensusScratch`, `PoaScratch`, `PolishScratch`):
  reusable workspaces eliminate repeated `Aligner`, `Topology`, `Workspace`, and DP buffer allocations across rounds,
  clusters, and repeated calls.
- **Polishing** (`plegma/polish.mojo`, `plegma/bands.mojo`, `plegma/likelihood.mojo`, `plegma/run_lengths.mojo`):
  edit-distance, pair-HMM likelihood, and homopolymer stages share `PolishScratch` and direct CIGAR parsing in
  `_Paths.fill`, with multi-threaded parallel read alignment, segment-level parallelism, and deterministic bit-exact
  read-level parallelism when a single segment dominates.

### Reads fixed at one anchor (`plegma/poa.mojo`)

A local assembly's reads span both anchors (`BOTH`) or are fixed at one (`LEFT`, `RIGHT`; `anchor=` in the FASTA
header). The spanning reads build the graph first; a left read is then aligned in the kernel's `PREFIX` mode (from
a source and its first letter, every letter, ending at any node; letters past a sink free, so they join as new nodes
and grow the graph), a right read in `SUFFIX` (from any node, letters before a source free, ending at a sink). Both
modes run in the adaptive band (`PREFIX` without the diagonal toward an end it does not have), and under `-S` only
the first or last pinned window takes the read's mode. A one-sided read is aligned only to the part of the graph
between the anchors (`_between`: what reaches the spanning reads' ends, or is reached from their starts): pinned,
the pins' windows ending at any sink, and its letters on nodes outside that part inserted (`_past_anchor`: the read's
letters past its own anchor, one run at its anchored end); unpinned, aligned to that part alone.
The consensus runs between the spanning reads' ends (`Graph.heaviest_path_spanning`), every read weighing the
edges. Before any of it, each spanning read's first and last `END_LETTERS` are checked against the reads anchored at
that end (`_ends_checked`), where one-sided reads are: one whose end most of them disagree with is anchored at its
other end alone, or, at neither, builds nothing and joins the MSA after. A one-sided side's anchored end is voted on
over `ANCHOR_SPREAD` nodes (`Graph._spread_vote`), as reads' anchored ends differ by a few letters.

With no spanning read each side is assembled alone, the longest read first: the left consensus followed forward
while most reads go on (`Graph.heaviest_path`), the right one backward (`Graph.heaviest_path_backward`). Both sides:
joined where they agree (`_join`): each way to place the right side over the left's end by 200 letters or more
(minimap2's pins between them on a 19-mer unique in both, grouped by diagonal; karan overlaps of their ends) is
aligned from its meeting point over the whole overlap, the fewest edits a letter wins, at its best-supported pin;
where another placement off its path agrees nearly as well the sides stay apart. Inside a tandem array (`_period`:
the distance at which the overlap's 8-mers recur) the candidates are the copy numbers: the meeting point moved by
whole copies to the same place in another copy (`_same_phase`), each aligned out from there, compared on the left
letters both align (`_copy_join`); the copy number is taken only where every other agrees clearly worse. The graph is rebuilt on the joined sequence. One side: its consensus cut
where fewer than `-c` reads (5) reach, each read's reach the consensus columns its letters are in (`_reach`; a
one-directional assembly). Polishing aligns each read as it is anchored (`polish`'s `anchors`: karan's `PREFIX` or
`SUFFIX`, or the consensus inside a longer read's prefix or suffix), counts disagreement against the reads covering
each position, and lets a read with a free end vote on the letters it covers: in its forward and backward matrices
the consensus before its start and after its end costs nothing, and an edit there costs it nothing either. The
likelihood stage, which does not model free ends, leaves a segment some read covers only part of to the edit
distance.

## Several consensus sequences of one read set (the library, the C API)

- `plegma/clustered.mojo` `clustered_consensus`: the spanning reads first split as `cluster_rows` splits an MSA's rows
  (with none spanning and one side only, that side's reads, `MISSING` past each one's far end, each cluster then
  weighed over the letters all reach, `_cut`),
  then k-means-like rounds (each read to the consensus it differs from least, each consensus made again), then
  clusters merged that differ only by run lengths the reads fit as one, by one short tandem unit with too few reads,
  or not at all.
- `plegma/phased.mojo` `phased_consensus`: reads already phased; each haplotype from its labelled reads, then again with
  the unlabelled reads on the haplotype whose letters they carry where the two differ (`_sides`).
- Shared: `plegma/placement.mojo` (`Clustered`, `Placement`: each cluster's consensus, made as `consensus` makes
  one, the longer side where its two sides could not be joined (`Clustered.joined` False), and every read's CIGAR to
  its cluster's) and `plegma/runs.mojo` (the homopolymer runs two consensus sequences
  differ by, and the lengths a model finds most likely for their reads). `clustered.mojo` re-exports the shared
  names and `phased_consensus` for callers that import them from it.

## `plegma ref`: the MSA anchored on a reference

`plegma/reference.mojo` `on_reference`, for pgrealign: every distinct sequence's letters paired with the reference's
by its own alignment to the reference, then placed in the reference letters' columns.

- **The pairings** (`_pairings`): independent, spread over `-t` threads, an aligner each, the longest first. Exact
  (karan, asm5's scores made costs: mismatch 40, a gap of k letters min(78 + 7k, 162 + 3k), ties left), or with the
  reference's flanks pinned as minimap2 pins it (`plegma/seeds.mojo` `ReferenceIndex` built once and shared read-only
  across threads, `plegma/chain.mojo` `rmq_chains`, `plegma/pins.mojo` `pinned_pairing_with_index`: minimap2 2.31's
  minimizers, seed filters, RMQ chaining with its float arithmetic, re-chaining, bad-end and bad-seed filters, pins at
  the minimizers' middles at dipcall's settings; karan exact between consecutive pins; then, as minimap2 last does,
  each gap shifted left through the matched run before it, `_left_aligned`). The port is checked against minimap2's own
  chains (`--print-chains`) and against a Python port, `experiments/chr20/mm2like.py`.
- **Fallback** (`_fall_back`): when most sequences replace stretches of the reference and those stretches are much of
  their letters (reads of one sequence, not alleles), those are aligned again at pgrealign's costs in parallel over
  `-t` threads.
- **Columns** (`_on_reference_columns`): each letter on its reference letter's column; the runs different sequences
  insert at one place aligned to each other (recursively, on the longest) in parallel across independent insertion
  slots over `-t` threads, so their copies share columns.

`rmq_chains` in `chain.mojo` is the hot loop of `ref -F`: minimap2's RMQ chaining with the inner window's anchors in a
linked list over a count tree, one 32-byte record per anchor, a tabulated gap penalty and an outer segment tree of
(priority, leaf) pairs.

## Tracing

`PLEGMA_TRACE` (comma-separated topics, or `all`) writes the decisions behind a result on stderr, one whole line each
(`plegma/trace.mojo`): `poa` (reads by anchor, order, another copy's reads, reverse-complemented reads, the two
sides' join or cut), `pins` (chains, pins kept, windows), `polish` (outliers, each round's edit distance and edits or
revert, each stage's length), `clustered` (first split, moves, each merge and why), `phased` (where unlabelled reads
went), `reference` (mode, fallback). The variable is read once per top-level call; off, a trace point is one branch.

```sh
PLEGMA_TRACE=poa,polish plegma poa -N reads.fa > /dev/null
```

## Testing and evaluation

- `pixi run test`: unit tests per module (`tests/`); `pixi run fuzz`: random inputs through both commands, every
  output's invariants checked, and the command line's own behaviour (gzip and CRLF input, bad options) first
  (`experiments/validate/fuzz.py`); `pixi run test-all`: those, the C API (shared and static) and Python.
- Stage oracles kept as scripts: `experiments/validate/anchored_grid.py` (simulated local assemblies against their
  truth: joins, cuts, distance) and `experiments/validate/minimap2_cigars.py` (`ref -F`'s pair CIGARs against
  minimap2's own).
- Refactors: `experiments/validate/golden.py` writes every command's and library call's outputs over fixed inputs;
  run it before and after a change meant to keep behaviour and `diff -r` the two directories.
- Consensus accuracy: `experiments/msa/accuracy.py` (simulated reads and real VNTR read sets against abPOA and
  others).
- Genotyping (`ref`): `experiments/chr20` runs pgrealign's realignment with plegma on chr20 of the HPRC graph and calls
  HG002 with vg (`plegma_arm.sh` for the whole chromosome; `evaluate.py` for 387 sampled regions in minutes;
  `qc_loop.py` for pgrealign qc's pair F1 in 90 s, co-optimal with `build/coopt`).
