# plegma's design

Two commands share one library. This file maps the code to them, for whoever changes it next.

## `plegma poa`: the general POA

`cli/plegma_cli.mojo` `_poa` adds the sequences one at a time, in abPOA's guide order (`-p`, the default,
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
- **Consensus**: abPOA's heaviest bundle (`Graph.heaviest_path` over the member paths), then polished in global mode
  (`plegma/polish.mojo`): one-letter and whole-repeat-unit edits kept while they lower the summed edit distance of the
  sequences to it, then while they raise their likelihood under a pair HMM. `-d` clusters the sequences on the MSA's
  heterozygous columns (`plegma/clusters.mojo`), one consensus per cluster.

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

## Testing and evaluation

- `pixi run test`: unit tests per module (`tests/`); `pixi run fuzz`: random inputs through both commands, every
  output's invariants checked (`experiments/validate/fuzz.py`).
- Consensus accuracy: `experiments/msa/accuracy.py` (simulated reads and real VNTR read sets against abPOA and
  others).
- Genotyping (`ref`): `experiments/chr20` runs pgrealign's realignment with plegma on chr20 of the HPRC graph and calls
  HG002 with vg (`plegma_arm.sh` for the whole chromosome; `evaluate.py` for 387 sampled regions in minutes;
  `qc_loop.py` for pgrealign qc's pair F1 in 90 s, co-optimal with `build/coopt`).
