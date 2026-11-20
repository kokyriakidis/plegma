# plegma's design

Two commands share one library. This file maps the code to them, for whoever changes it next.

## `plegma poa`: the general POA

The command line is `cli/`: `options.mojo` (the usage text, parsing, every value checked, the library's own checks
through `ConsensusOptions.check`), `records.mojo` (the input read once, gzip or not, FASTA or FASTQ, with each record's
`anchor=` and `strand=` header words), `gfa.mojo` (`-r 3`, `-r 4`) and `plegma_cli.mojo` (each command, the outputs;
an error is one `plegma: ...` line and exit status 1).

`plegma.poa` (`plegma/poa.mojo`) adds the sequences one at a time, in abPOA's guide order (`-p`, the default,
`plegma/order.mojo`), input order (`-n`) or longest first (`-L`).

- **Alignment to the graph** (`plegma/kernel.mojo`): banded DP over the nodes in topological order, abPOA's scoring
  (match score, mismatch, convex two-piece gaps) and adaptive band (`b + f*L` either side of the predecessors' best
  cells and of the diagonal the node's heaviest remaining path leaves). Memory follows the band: a row lives until its
  last successor is done, and the traceback keeps one byte per cell (three more at nodes with several predecessors).
  `Workspace.align` aligns between any two nodes (`SOURCE` and `SINK` for the ends), so a sequence can also be aligned
  window by window between anchors.
- **Seeded** (`-S`, the default; `-N` aligns to the whole graph; `plegma/pinned.mojo` `align_pinned`): minimap2's
  pins (`pins.pins`) between the sequence and the consensus path, kept only where the 19-mer centred on a pin is
  unique in both (in near-identical tandem copies the chain may pair the wrong copy), each pin a one-letter anchor on
  its path node; the kernel then aligns each stretch between anchors to the subgraph between their nodes
  (`window_nodes`), in the band the whole sequence would get. Where no pin survives, the one window is the whole graph.
- **The graph** (`plegma/graph.mojo`): flat arrays of nodes, columns and edges. A sequence joins along its aligned
  nodes (`fuse_nodes`); an MSA's rows make a graph directly (`from_msa`); `msa` writes the rows back.
- **Consensus** (`plegma/consensus.mojo` `polished`, the one path `plegma.consensus`, the C API and `plegma poa -r 0`
  share): abPOA's heaviest bundle (`Graph.heaviest_path` over the member paths), then polished in global mode.
  `-d` clusters the sequences on the MSA's heterozygous columns (`plegma/clusters.mojo` `cluster_rows`), one
  consensus per cluster; `-a 1` takes each column's most frequent letter.
- **Polishing** (`plegma/polish.mojo` `polish`, the edit distance stage, driving the others): one-letter and
  whole-repeat-unit edits kept while they lower the summed edit distance of the sequences to it; then
  `plegma/likelihood.mojo` (`_likely`), while they raise their likelihood under a pair HMM; then, with a homopolymer
  model, `plegma/run_lengths.mojo` (`_runs`), every run at its most likely length. All three work on what
  `plegma/bands.mojo` holds: each read aligned to the consensus (`_Paths`, a one-sided read as it is anchored,
  `past_the_end`) and its forward and backward matrices in a band around that alignment (`_Bands`).

### Reads fixed at one anchor (`plegma/poa.mojo`)

A local assembly's reads span both anchors (`BOTH`) or are fixed at one (`LEFT`, `RIGHT`; `anchor=` in the FASTA
header). The spanning reads build the graph first; a left read is then aligned in the kernel's `PREFIX` mode (from
a source and its first letter, every letter, ending at any node; letters past a sink free, so they join as new nodes
and grow the graph), a right read in `SUFFIX` (from any node, letters before a source free, ending at a sink). Both
modes run in the adaptive band (`PREFIX` without the diagonal toward an end it does not have), and under `-S` only
the first or last pinned window takes the read's mode. A one-sided read is aligned only to the part of the graph
between the anchors (`_between`: what reaches the spanning reads' ends, or is reached from their starts), unpinned
once a read has grown the graph past the far anchor, as the pins' windows end at any sink. The consensus runs between the spanning reads' ends
(`Graph.heaviest_path_spanning`), every read weighing the edges.

With no spanning read each side is assembled alone, the longest read first: the left consensus followed forward
while most reads go on (`Graph.heaviest_path`), the right one backward (`Graph.heaviest_path_backward`). Both sides:
joined at the middle pin of minimap2's chain between them (`pins.pins`), or a short karan overlap of their ends; the
graph is rebuilt on the joined sequence. One side: its consensus cut where fewer than `-c` reads (5) reach (a
one-directional assembly). Polishing aligns each read as it is anchored (`polish`'s `anchors`: karan's `PREFIX` or
`SUFFIX`, or the consensus inside a longer read's prefix or suffix), counts disagreement against the reads covering
each position, and lets a read with a free end vote on the letters it covers: in its forward and backward matrices
the consensus before its start and after its end costs nothing, and an edit there costs it nothing either. The
likelihood stage, which does not model free ends, leaves a segment some read covers only part of to the edit
distance.

## Several consensus sequences of one read set (the library, the C API)

- `plegma/clustered.mojo` `clustered_consensus`: the spanning reads first split as `cluster_rows` splits an MSA's rows,
  then k-means-like rounds (each read to the consensus it differs from least, each consensus made again), then
  clusters merged that differ only by run lengths the reads fit as one, by one short tandem unit with too few reads,
  or not at all.
- `plegma/phased.mojo` `phased_consensus`: reads already phased; each haplotype from its labelled reads, then again with
  the unlabelled reads on the haplotype whose letters they carry where the two differ (`_sides`).
- Shared: `plegma/placement.mojo` (`Clustered`, `Placement`: each cluster's consensus, made as `consensus` makes
  one, and every read's CIGAR to its cluster's) and `plegma/runs.mojo` (the homopolymer runs two consensus sequences
  differ by, and the lengths a model finds most likely for their reads). `clustered.mojo` re-exports the shared
  names and `phased_consensus` for callers that import them from it.

## `plegma ref`: the MSA anchored on a reference

`plegma/reference.mojo` `on_reference`, for pgrealign: every distinct sequence's letters paired with the reference's
by its own alignment to the reference, then placed in the reference letters' columns.

- **The pairings** (`_pairings`): independent, spread over `-t` threads, an aligner each, the longest first. Exact
  (karan, asm5's scores made costs: mismatch 40, a gap of k letters min(78 + 7k, 162 + 3k), ties left), or with the
  reference's flanks pinned as minimap2 pins it (`plegma/pins.mojo` `pinned_pairing`: minimap2 2.31's minimizers,
  seed filters, RMQ chaining with its float arithmetic, re-chaining, bad-end and bad-seed filters, pins at the
  minimizers' middles at dipcall's settings; karan exact between consecutive pins). The port is checked against
  minimap2's own chains (`--print-chains`) and against a Python port, `experiments/chr20/mm2like.py`.
- **Fallback** (`_fall_back`): when most sequences replace stretches of the reference and those stretches are much of
  their letters (reads of one sequence, not alleles), those are aligned again at pgrealign's costs.
- **Columns** (`_on_reference_columns`): each letter on its reference letter's column; the runs different sequences
  insert at one place aligned to each other (recursively, on the longest), so their copies share columns.

`_chains` in `pins.mojo` is the hot loop of `ref -F` (about half its time): minimap2's RMQ chaining with the inner
window's anchors in a linked list over a count tree, one 32-byte record per anchor, a tabulated gap penalty and an
outer segment tree of (priority, leaf) pairs.

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
- Refactors: `experiments/validate/golden.py` writes every command's and library call's outputs over fixed inputs;
  run it before and after a change meant to keep behaviour and `diff -r` the two directories.
- Consensus accuracy: `experiments/msa/accuracy.py` (simulated reads and real VNTR read sets against abPOA and
  others).
- Genotyping (`ref`): `experiments/chr20` runs pgrealign's realignment with plegma on chr20 of the HPRC graph and calls
  HG002 with vg (`plegma_arm.sh` for the whole chromosome; `evaluate.py` for 387 sampled regions in minutes;
  `qc_loop.py` for pgrealign qc's pair F1 in 90 s, co-optimal with `build/coopt`).
