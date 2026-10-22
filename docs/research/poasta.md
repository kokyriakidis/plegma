# POASTA review for plegma

Target: https://github.com/broadinstitute/poasta. I read the default branch `stable` at commit `8186a1f` ("release: Version 1.0", 2026-10-18). This **v1.0.0** is a rewrite released three days before this review. The paper describes the older engine (tag `paper`, v0.1.0). I also read the `banded-poa`, `po-poa` and `main` branches, the GitHub issues, and the paper (van Dijk et al., *Bioinformatics* 2025, PMC11755094). I did not build or run anything.

Paths below are relative to `/private/tmp/claude-501/-Users-kokyriakidis/8500d74c-03af-4a3f-8699-f1899147a1f9/scratchpad/src/poasta` (stable). Other branches are cloned next to it as `poasta-banded-poa`, `poasta-po-poa` and `poasta-main`. **[inferred]** marks conclusions I reasoned out but did not read directly in code or the paper.

---

## 0. One-paragraph summary

POASTA aligns one sequence at a time to a POA DAG. It treats alignment as a shortest-path problem over the product graph (node × query offset × gap state) and solves it with A*. The A* uses a bucket queue and a greedy, depth-first, zero-cost match extension, which is the graph analogue of WFA's `extend`. It opens gaps only where a mismatch occurs, and it prunes states using superbubble exits. The default heuristic, "minimum gap cost", charges only for the length difference between the rest of the query and the range of remaining path lengths [lo, hi] to the graph end. v1.0 adds an optional exact gap-chaining seed heuristic (after A*PA) that uses k=8 seeds, plus linear, affine and two-piece affine costs and semi-global mode.

Results are exact per pairwise step. Exactness is enforced by a randomized differential test against a full Gotoh DP, which includes repeat-alphabet graphs. There is no seeding or banding in the default path, no guide tree, no anchoring and no cap on work. Sequences are added in input order.

For VNTRs, the default heuristic degenerates toward Dijkstra wherever haplotypes differ in copy number, because the [lo, hi] path-length interval is wide **[inferred]**. The chain heuristic must index every seed occurrence, which grows roughly quadratically in a tandem array **[inferred]**. Tie-breaking between equally optimal gap placements is an arbitrary but deterministic order.

---

## 1. Graph representation and memory layout

**Store and view split** (`src/graph/poa.rs`):
- **Store.** It only grows (`POAGraph`, poa.rs:156-187).
  - `symbols: Vec<u8>`: one base per node, so nodes are single characters.
  - `aligned_nodes: Vec<Vec<NodeId>>`: the MSA-column cluster, as in SPOA.
  - `edges: Vec<POAEdge>`, with an `edge_ix: FxHashMap<(src,tgt), usize>` for dedup.
  - Each edge carries `POAEdgeData { weight: usize, sequence_ids: Vec<usize> }` (poa.rs:67-71). That is one heap `Vec` per edge, with 8 bytes per traversing sequence.
  - `add_edge` appends the sequence id and adds weight (poa.rs:349-377).
- **View.** It is rebuilt from scratch by `post_process()` after **every** added sequence (poa.rs:417-541).
  - `succ`/`pred` are CSR `Adjacency { start, items, edges }` indexed by **topological rank** (poa.rs:107-139, 485-516).
  - It also holds `rank_symbols`, plus `depth_min`/`depth_max` from the start (poa.rs:518-540).
  - The rank is a reverse post-order from an iterative DFS (poa.rs:556-586).
  - Sentinels `#` (start) and `$` (end) bracket the graph. Their edges exist only in the view (poa.rs:422-475).
- **Edge order matters.** "Newest edge first" (poa.rs:476-481). This decides the topological order and which of several equally optimal alignments is returned (CLAUDE.md "Invariants"; architecture.md:114-116).
- **Index types.** Node index type `Ix = u32` by default. The visited storage requires `node.index() == rank` (`assert_index_is_rank`, visited.rs:426-433).
- **Folding an alignment in** (`add_alignment`, poa.rs:683-818), SPOA-style:
  - A match reuses the node.
  - A mismatch reuses a node of the same symbol in the `aligned_nodes` cluster, or creates one and links it into the cluster (poa.rs:750-777).
  - An insertion creates an unlinked node (poa.rs:778-782).
  - Unaligned leading and trailing bases become fresh chains (poa.rs:723-736).
- **Superbubble index.** `src/graph/bubbles/`:
  - Gärtner et al. 2018 finder on the DFS topological order (finder.rs:1-60).
  - `BubbleIndex` holds per-node `(min,max)` distance to the end in one reverse-rank pass (index.rs:84-98).
  - It holds a CSR node→bubble membership with distance to each exit. This is computed as `dist_to_end(node) − dist_to_end(exit)` (index.rs:100-122), using the property that superbubbles are rank intervals and are nested or disjoint (index.rs:180-219).
  - **It is rebuilt for every alignment in the POA flow** (astar/mod.rs:296-299; index.rs:5-8). lasagna builds it once.
- **A* state storage** (`src/align/astar/visited.rs`):
  - `BlockStorage<N,O,S,B=8>`: one `BlockRow` per group of 8 ranks (visited.rs:292-302).
  - Each row is a **contiguous window** of 8×8 blocks over query offsets. Each block cell is `[u32; S]`, with S = 1, 3 or 5 states (visited.rs:91-97).
  - The window starts at 2 blocks and doubles toward whichever side is touched (visited.rs:137-154).
  - A lookup is a shift, a subtract, a bounds check and a load (visited.rs:109-118).
  - Unreached cells read as `INF = 1<<28` (score.rs:18).
  - Separately, `ReachedBits` keeps one bit per query offset per rank for reached `Match` states. Pruning uses it (visited.rs:172-286).
  - `reset` clears only rows touched by the previous alignment, and the engine keeps its storage between alignments (visited.rs:335-353; astar/mod.rs:183-191).
  - Memory consequence **[inferred]**: memory is per 8-rank group × (max offset touched − min offset touched). A rank group touched at two distant offsets pays for the whole span between them. Off-diagonal exploration in repeats therefore inflates memory toward a band of the full span, not just the cells actually visited.
- **Queue.** `LayeredQueue` is a bucket queue (a deque of layers indexed by `f = g + h`) with O(1) push and pop. Each layer has one vector per state (queue.rs:28-104). Within a bucket, pops are LIFO with Match first, then deletions, then insertions (architecture.md:223-231).
- **Integer widths.**
  - Costs are `u8` (cost_models/mod.rs:13-20; aln_graph.rs:422-428).
  - Scores are `u32`, capped by `INF = 2^28`.
  - Query offsets use a generic `O: OffsetType` (`u16`/`u32`/`u64`; the CLI uses `u32`, astar/mod.rs:51-54).

## 2. The algorithm

### 2.1 State graph and main loop
- **States.** A state is `(node rank, query offset)` plus an `AlignState` (Match, Deletion, Insertion, Deletion2, Insertion2) (aln_graph.rs:16-43; align/mod.rs:94-102).
  - Linear costs use 1 state, affine 3, two-piece 5 (astar/mod.rs:91-161).
  - The two pieces are independent tracks: `D2 = min(M+go2+ge2, D2+ge2)`, with no transitions between pieces (aln_graph.rs:430-438).
- **Main loop** (`run_astar_search`, astar/mod.rs:340-540):
  1. Seed the initial states (372-376).
  2. Pop a state. Drop it if a better score is already stored (395-401).
  3. If the heuristic is dynamic, recompute h and requeue the state if its priority rose (405-419).
  4. Stop at the end state (421-425).
  5. Drop the state if bubble-exit pruning says it cannot improve (427-431).
  6. Mark it reached and expand it.
- **Match states are never expanded one step.** They are handed to the greedy DFA (444-515). Gap states go through `expand_all`, which closes the gap into Match for free and extends it by one (aln_graph.rs:494-534).

### 2.2 Depth-first greedy match extension (WFA `extend` on a DAG)
- `DepthFirstGreedyAlignment` (dfa.rs:99-269) keeps a stack of `(state, successor iterator)`. It walks **every** matching successor depth first at the same score and marks each one visited.
- It returns to the caller only on a mismatch, at the graph end, at the query end, or at an ends-free end state (dfa.rs:158-218, 220-268).
- A successor whose stored score is already ≤ the current score is skipped silently (dfa.rs:253-261).
- Matched states are also pruning-checked inside the DFA (dfa.rs:181-185).
- **Gap lemma.** Indel states are opened **only at mismatch points** (`expand_mismatch`, aln_graph.rs:591-642), at the query end (563-589) or at the graph end (536-561).
  - The paper's justification: match edges cost 0, so where an unvisited match successor exists, the insertion and deletion edges from the current state can be ignored.
  - Because the DFA keeps iterating over the remaining successors after it returns a Mismatch, gaps toward mismatching siblings are still opened (astar/mod.rs:451-498).

### 2.3 Heuristic: minimum remaining gap cost
`MinimumGapCost::bound` (heuristic/min_gap.rs:127-183):
- Inputs: `lo = min_dist_to_end(v) − 1`, `hi = max_dist_to_end(v) − 1`, query offset `i`, query length `m`.
  - If `i + lo > m`, deletions are unavoidable, so `h = gap_cost(state, i+lo−m)`.
  - If `i + hi < m`, insertions are unavoidable, so `h = gap_cost(state, m−i−hi)`.
  - **Otherwise h = 0.**
- No gap-open is charged if the state is already in a gap of the matching kind (min_gap.rs:159-176). Under ends-free modes the bound is 0 inside free regions, which keeps it admissible (min_gap.rs:148-171).
- The heuristic knows nothing about sequence content: mismatches are never charged.
- On a graph whose remaining paths span a wide range of lengths, h is 0 for every state whose offset lies between those bounds **[inferred from the formula]**.

### 2.4 Gap-chaining seed heuristic (new in v1.0, opt-in with `--heuristic chain`)
- `SmoothHeuristic` (heuristic/chain/mod.rs:186-396) is A*PA's smooth gap-chaining seed heuristic adapted to DAGs and affine costs.
- **Seeds.** The query is cut into seeds of k bases (default 8, CLI range 4-32, r ∈ {0,1} substitutions; cli/mod.rs:100-117).
- **Seed matching.** A sweep in rank order with a seed trie finds **every** k-node path that spells a seed. The doc says the index must be complete for admissibility (chain/finder.rs:1-12, 103+). The same doc notes up to 1,487 k-paths at one node on a dense graph.
- **Points and chains.** Matches become points `T = (e·(lo+i) − P, −e·(hi+i) − P)`. h is the seed potential minus the best chain credit, found by a 2D dominance query over "contours" (chain/mod.rs:6-18, 162-173; contours.rs:1-11). The result is maxed with min-gap-cost (chain/mod.rs:333-352).
- **Match pruning.** Matches behind the search front are deleted, which makes h dynamic (`on_expand`, chain/mod.rs:355-391). Without it the heuristic "is quadratic on long reads" (cli/mod.rs:122-124).
- **Charge per seed.** min(miss, mismatch, k·e, e + ⌊o/2⌋) for consistency (chain/costs.rs:44-55).
- **Proofs.** The module points to `docs/chain-heuristic-proofs.md` and `docs/chain-heuristic-explained.md`. **Neither file is in the repo**: `docs/` holds only `architecture.md`.

### 2.5 Superbubble-exit pruning ("superfluous states")
- **Idea.** Every path from a node inside a superbubble passes through the bubble's exit. If the exit has already been reached in Match at offsets `R` with known scores, a candidate `(v,i)` can only matter if it can reach the exit at an offset and score that no reached entry dominates.
- **Implementation** (`ReachedBubbleExitsMatch::can_improve_bubble`, pruning.rs:63-239):
  - The candidate's exit window is `[i+min_dist_to_exit, i+max_dist_to_exit]`.
  - The check scans the exit's `ReachedBits` below, inside and above that window (pruning.rs:105-222).
  - It compares the candidate's score with the reached exit score plus an implicitly opened gap (`can_improve_at_offset`, pruning.rs:264-319; `implicit_gap_cost`, 247-261).
- **Gap-state cases** have extra rules.
  - A candidate already in a gap pays no reopen, so the exit must pay the open of that piece.
  - A deletion candidate can delete straight through the exit (pruning.rs:75-99, 130-158).
- `prune` takes the conjunction over all bubbles containing the node (visited.rs:517-541).
- Both the DFA and the queue pops apply pruning.
- Per the paper, this pruning was the fix for issue #20, "quadratic (or even exponential) behaviour with increasing alignment score".

### 2.6 Scoring, modes, exactness
- **Cost models.** Linear, affine and two-piece affine, all with `u8` costs. CLI defaults are mismatch 4, open 6, extend 2, open2 24, extend2 1 (cli/poasta.rs:87-110). The match cost is fixed at 0. There is **no substitution matrix** (issue #24, open, requested by the Cactus developers).
- **Modes** (align/mod.rs:79-91; aln_graph.rs:159-277):
  - `Global`.
  - `EndsFree` with four bounded budgets.
  - The CLI's `semi-global` means a free graph prefix and suffix with the query aligned end to end (cli/mod.rs `AlignmentSpan::mode`).
  - Global mode charges end gaps, which caused the "strange results" in issue #25.
- **Exactness claim.** The result is the optimal sequence-to-graph alignment under the cost model, provided match = 0 and all other costs are non-negative (paper).
  - Pruning has no formal lemma in the paper ("argued in the text and the supplement").
  - In v1.0, exactness rests on `DpOracle`, a full-matrix Gotoh DP for all models in global and semi-global modes (astar/dp_oracle.rs), and on a randomized differential test.
  - The differential test runs every heuristic, seed setting and pruning setting, checks admissibility and consistency on every edge, and draws half its cases from the **AC alphabet to create repeats** (differential_tests.rs:79-90).
  - The suite carries regression tests for real exactness bugs in pruning: two-piece pricing pruned an optimal first state, and deletion candidates were pruned without a reopen charge (astar/mod.rs ~1301-1350).
  - So exactness has been violated before in development and is enforced by tests, not by proof **[the two-piece cases post-date the paper; whether the v0.1 affine engine had analogous bugs, I did not check]**.
- **Backtrace** uses scores only, with no traceback pointers (visited.rs `gap_step` 718-809, `linear_step` 587-642). At a Match cell it prefers the diagonal predecessor, then the free close from Deletion, then from Insertion (visited.rs ~741-773). Predecessors are tried in newest-edge-first order. This plus LIFO queue order defines which co-optimal alignment comes out (CLAUDE.md invariants).

### 2.7 Memory use and high divergence
- **Memory** scales with the visited band, as described in §1. Paper Table 1 (M. tuberculosis, 342 genomes, ANI 99.3%):

  | Length | Time | Peak memory |
  |---|---|---|
  | 250 kbp | 5.3 h | 63.8 GB |
  | 500 kbp | 24 h | 120 GB |
  | 1 Mbp | 69 h | 231 GB |

  SPOA and abPOA could not run these; the authors estimate 0.95, 3.5 and 13 TB.
- **High divergence.**
  - The paper states POASTA loses to SPOA when sequence diversity is high: the only loss of 13 sets was at ANI 82.6%. There are fewer exact matches for the greedy extension, and the A* is "less predictable and less cache-friendly" than SPOA's SIMD full matrix.
  - With min-gap-cost, h ignores mismatches, so the explored region grows with roughly the optimal cost divided by the gap extend cost, around every diagonal that admissible paths allow **[inferred]**.
  - Nothing caps the work: no state budget, no timeout, no fallback (grep finds none).

## 3. Seeding and banding

- **Default path.** No seeding and no banding. Exactness comes from A* plus pruning.
- **`--heuristic chain`.** Uses seeds (exact or 1-substitution k-mers), but only to **compute h**, never to restrict the search. Every seed occurrence must be found (chain/finder.rs:3-12), so there is no frequency masking. That is the opposite of what repeat-heavy input wants.
- **`banded-poa` branch** (2026-04, not merged; `poasta-banded-poa/src/align/engine/band_doubling.rs`):
  - A scalar banded Gotoh POA with **Ukkonen-style band doubling** and per-node multi-band intervals propagated from predecessors (band_doubling.rs:1-24).
  - Optimality certificate: accept the band-k result if `k ≥ forced + score/min_gap_extend + 1`, where `forced` is the length difference that indels must absorb. Otherwise double k, and at saturation run full width (band_doubling.rs:137-184).
  - So the authors are themselves exploring band-doubling DP as an alternative to A*.
- **No anchors and no MUMs.** Each sequence is aligned to the whole graph.

## 4. Progressive strategy, order dependence, outputs

- **Order.** Sequences are aligned strictly in input order, one at a time, single-threaded. The first sequence seeds the graph (align/mod.rs:295-333; bin/poasta.rs:220-258). There is no guide tree, no re-alignment and no iterative refinement, so the result depends on order.
- **Seed graph.** An existing FASTA MSA can be loaded as a starting graph (`-I/--graph`; io/fasta.rs `load_graph_from_fasta_msa`).
- **Weights.** FASTQ quality scores become per-base weights, floored at 1 (io/seq.rs:59-118).
- **PO-PO (graph-to-graph) alignment.** WIP prototype on the `po-poa` branch: a dense three-matrix DP (`poasta-po-poa/src/align/engine/dp/popoa/mod.rs:24-42`). PR #35 was closed with no completed items. It is not available.
- **Outputs.**
  - **FASTA MSA.** Columns are clusters from `aligned_nodes`, assigned in rank order (io/fasta.rs ~140-162). Issue #34, ragged MSA rows on a TC microsatellite, is fixed in v1.0, and `tests/microsat.fa` is the regression fixture.
  - **GFA 1.0.** **One S segment per base**, with no unchopping. It has L lines per edge and P lines walked through `sequence_ids` binary search (io/gfa.rs:48-98). Issue #29 asks for nodes split so paths are preserved exactly.
  - **Consensus.** Lee heaviest bundle (graph/consensus.rs:1-60). Issue #37 (open) reports it running out along single-read overhangs because sentinel edges weigh 0. The STRdust author (tandem-repeat genotyping) had to trim terminal nodes below 0.35 support.
  - The README still documents a `.poasta` binary format, DOT output and `view`. In v1.0, `view` and `stats` are "not yet implemented" (bin/poasta.rs:54-61), and the README is stale.
- **lasagna** aligns reads to a fixed GFA DAG with multi-threading and writes GAF (bin/lasagna.rs). Cyclic GFA is rejected (io/gfa.rs tests).

## 5. Tandem repeats

Nothing in the paper or the code targets tandem repeats. The paper does not evaluate repeat-rich data. Everything below is **[inferred]** from the code unless cited.

1. **Min-gap-cost collapses to Dijkstra in VNTRs.**
   - h is non-zero only when the remaining query lies outside the [lo, hi] interval of remaining path lengths (min_gap.rs:148-180).
   - A VNTR POA graph built from haplotypes with different copy numbers has hi − lo roughly equal to the copy-number range times the unit length, so h = 0 for almost every state in the array.
   - The search then explores every state with g < OPT, inside the reachable diagonals: about N_nodes × (OPT/ge + number of distinct path-length shifts).
   - Bubble pruning engages only at superbubble exits. Overlapping copy-number branches tend to form one large complex region rather than many small superbubbles, so exits are rare inside the array.
2. **Greedy extension does not go into the "wrong copy"; it goes into all copies.**
   - The DFA follows every matching successor (dfa.rs:228-267), not one.
   - In a periodic region, the zero-cost closure from a single popped state covers every diagonal reachable through branch choices that keep phase. These are many states, each written to visited storage and checked for pruning.
   - This is exact, but the work per pop is large, and block-row windows widen, which costs memory.
3. **Many equal-cost paths.**
   - A copy-number gap of length L·p has many co-optimal placements along the array, and the graph also allows zero-cost path switching between haplotypes.
   - The bucket queue does not blow up on ties: equal f values share a bucket, and LIFO within the bucket stays depth-first.
   - But which co-optimal placement is returned depends on LIFO order, newest-edge-first adjacency and the diagonal-first backtrace (CLAUDE.md; visited.rs ~741-773).
   - There is no left- or right-normalization of gaps and no unit-boundary awareness. Successive haplotypes can place the "same" copy-number event at different phases, which fragments the graph and the MSA columns.
4. **The chain heuristic is quadratic-ish in arrays.**
   - Every seed inside a VNTR unit matches every copy, so the number of matches is about (#seeds) × (#copies).
   - Example: a 30 kb array with a 30 bp unit and k = 8 gives about 3,750 seeds × 1,000 copies ≈ 3.75 M matches. All are required for admissibility (finder.rs:3-12).
   - Match pruning and contours handle them exactly, but setup cost and memory grow with the square of the array length.
   - Off-diagonal chains differ by a multiple of the unit, and wide [lo, hi] makes the shift penalty cheap, so the credit stays high and h stays weak.
5. **Is exact gap-affine optimality what we want?**
   - It optimizes the metric the evaluation uses, so as a **kernel** it is the right target: per sequence-to-graph step, POASTA returns the true minimum.
   - It fails in three places:
     - **Progressive.** The total graph cost depends on input order and earlier tie-breaks. Exactness per step does not give an optimal MSA.
     - **Ties.** Co-optimal placements are chosen by implementation order, not by homology (unit phase).
     - **Feasibility.** The optimum may be cheaper with gaps that split units, because mismatch-laden sliding onto a neighbouring copy can beat one long gap under 4/6/2. Two-piece affine reduces this, and POASTA supports it.
   - Field evidence: issue #34 (microsatellite MSA output, fixed). Issue #37: STRdust found poasta with linear costs matched rust-bio accuracy on simulated repeat reads while being **3-20× faster**, after the consensus fix.
6. **No guard rails.** Nothing caps or falls back when work explodes (§2.7).

## 6. Performance numbers (paper, v0.1 engine)

- **Data.** 13 sets of 140-2,385 sequences from five bacterial housekeeping genes (dnaG, nusA, pgk, pyrG, rpoB), mean length 1-4 kbp, ANI 82-97%.
- **vs SPOA (exact).**
  - POASTA was faster on 12 of 13 sets: **4.1× average, 9.8× maximum** speedup.
  - The one loss was the most divergent set, ANI 82.6%.
  - SPOA used **2.6× more memory** on average.
- **vs abPOA (adaptive band, heuristic).**
  - abPOA was **about 3.5× faster** than POASTA but used more memory on nearly all sets, because it allocates the full matrix.
  - abPOA was optimal in **99.8%** of alignments. Its misses "started or ended at unexpected nodes" and reported costs *below* the global optimum, which should be impossible. The authors attribute this to bands that miss indels larger than the band.
- **Scaling with length** (M. tuberculosis, Table 1): from 250 kbp to 1 Mbp, time rose about 13× and memory about 3.6×. That is clearly superlinear in time.
- **Scaling with divergence.** The paper says POASTA is best at ANI 90-100% with sequences over 1.5 kbp, and degrades as matches become scarce. The main text gives no numbers.
- **v1.0.** The CHANGELOG claims "twice as fast as the previous version". I found no published benchmark for v1.0. The `poa-bench` repo is listed as related.
- **Cost parameters** used in the benchmarks are not stated in the main text. Figures 2-3 use linear Δx = 4, Δg = 2. The CLI defaults are 4/6/2.

## 7. Newer development

Branches from `git ls-remote`:
- **`stable` = v1.0.0** (HEAD, default).
  - The rewrite: two-piece and linear costs, the chain heuristic, semi-global mode, consensus.
  - Block storage, reached bits and an O(V+E) bubble index.
  - The differential-test harness.
- **`main`** (2025-11): the older line, merged with `poa_seqgraph` (a sequence-labelled POA graph, PR #23).
- **`new-version`** (2026-09): the pre-release of v1.0.
- **`banded-poa`** (2026-04): band-doubling banded Gotoh with an Ukkonen optimality certificate (§3). Not merged. It signals that the authors see a deterministic banded DP as competitive with A*.
- **`po-poa`** (2025-11, PR #35 closed): an unfinished skeleton for graph-to-graph alignment.
- **`graph_formulation`, `intervals`** (2023): old experiments.
- **Sequence-to-graph** is already in-tree as `lasagna` (GFA in, GAF out, multi-threaded).

Overall: v1.0 gets a better heuristic for long reads and for reads shorter than the graph, and more cost models. **None of it adds anchoring, repeat awareness, gap normalization, order independence or work caps.** C/C++ bindings exist externally: `kokyriakidis/poasta-c` (issue #27).

## 8. Advantages and disadvantages for plegma

### Advantages
- Exact sequence-to-graph optimum under linear, affine or two-piece affine costs. This is the same objective plegma is evaluated on.
- Exactness is verified by a full-DP oracle plus a randomized differential test that includes repeat-alphabet graphs. The test methodology is reusable.
- Two-piece (convex) affine is supported, which matters for copy-number-sized indels.
- Greedy zero-cost extension plus opening gaps only at mismatches makes work proportional to the number of differences on similar sequences. That suits haplotypes that are near-identical outside the copy-number changes.
- Superbubble-exit dominance pruning is cheap: one O(V+E) index and bitmap scans.
- Memory-lean storage: 8×8 blocks, windowed rows, reached bits, score-only backtrace, reuse across alignments. Memory was 2.6× lower than SPOA in the paper.
- The bucket queue suits small integer costs. Output is deterministic.
- Ends-free budgets let plegma align a window with a free or fixed start and end exactly.
- Semi-global mode, FASTQ weights and consensus are all present.

### Disadvantages
- The default heuristic carries no information inside VNTRs (h ≈ 0 when the path-length range is wide), so the search falls back to Dijkstra.
- No seeding, anchoring or banding to bound the search. There are no caps or fallbacks, so the worst case is unbounded in practice and time is superlinear in length (13× time for 4× length).
- The chain heuristic needs every seed occurrence, which is about quadratic in tandem arrays.
- Worse than SPOA at high divergence.
- Order-dependent progressive alignment in input order, with no guide tree and no refinement. PO-PO alignment is not implemented.
- Co-optimal gap placement is decided by implementation order. There is no normalization or unit-phase awareness, so the same copy-number event can be represented inconsistently across haplotypes.
- The view and the bubble index are rebuilt in O(V+E) per added sequence. That is cheap relative to alignment, but repeated 456 times.
- Single-base nodes. GFA output is unchopped, one base per segment.
- Each edge keeps a `Vec<usize>` of sequence ids.
- `u8` costs and no substitution matrix.
- Consensus overhang bug (#37). README is stale. Proof documents referenced in the code are missing. Pruning exactness was broken and fixed during development.
- Single-threaded MSA construction.

## 9. Ideas to steal and pitfalls to avoid

### Ideas to steal
1. **Use an exact Dijkstra/A* WFA-style kernel inside sporos-bounded windows, in global mode.**
   - Once MUMs/SMEMs fix the window's entry and exit nodes and query offsets, run POASTA's loop on the subgraph between anchors with `Global` semantics.
   - The loop is a bucket queue, greedy DFA extension and gaps opened only at mismatch.
   - In non-repeat windows, [lo, hi] is narrow, so min-gap-cost is strong and the search is close to the diagonal.
   - Budget it, for example a cap on states expanded. On overflow, switch to the band-doubling DP below.
2. **Use the Ukkonen band-doubling certificate** from `banded-poa` (band_doubling.rs:137-151): accept the band-k result if `k ≥ forced + score/min_ge + 1`.
   - This gives a provable optimum with predictable, SIMD-friendly memory.
   - It is probably the better kernel inside long VNTR windows where A* has no useful heuristic **[inferred]**.
   - `forced` comes from window lo/hi versus query length, which the anchors already supply.
3. **Make the heuristic window-local.** Compute lo/hi (and min-gap h) to the *window exit anchor*, not the graph end.
   - Anchors shrink hi − lo, which directly strengthens min-gap-cost.
   - This is the cleanest reason sporos partitioning and an A* kernel compose well.
4. **Superbubble index in one rank pass**, with distance-to-exit computed as `dist_to_end(v) − dist_to_end(exit)` (index.rs:84-122, 180-219). Use it both for exit-dominance pruning and to pick partition boundaries: superbubble exits are natural cut points next to MUMs.
5. **Adopt the differential-testing harness as is:**
   - A full Gotoh oracle on small random DAGs built from AC-alphabet templates.
   - Assertions of admissibility and consistency for every edge.
   - Re-scoring every backtrace (differential_tests.rs; `test_support::score_from_alignment`, align/mod.rs:23-75).
   - Add VNTR templates: units 2-60 bp, varying copy counts.
6. **Storage tricks:**
   - Score-only backtrace (no pointers).
   - 8×8 blocks per rank group with a doubling window.
   - Reached bitsets with masked word scans (visited.rs:172-286).
   - Reset only the touched rows.
   - Keep scratch storage in the engine across windows. Each Mojo worker owns its engine.
7. **Gap-open-only-at-mismatch lemma and LIFO buckets** reduce queue traffic sharply on near-identical haplotypes.
8. **Two-piece affine** (independent tracks, aln_graph.rs:430-438) as the default for VNTR windows, so whole-unit copy-number gaps beat unit-splitting mismatch slides.
9. **Consensus fix from #37.** Trim terminal nodes with incoming support below a fraction of sequences, and do not give sentinel edges weight 0 semantics that reward overhangs.

### Pitfalls to avoid
1. **Don't rely on a length-only heuristic inside tandem arrays.** It goes to 0, and A* becomes Dijkstra with A* overhead. Anchor first, or use banded DP there.
2. **Don't index all k-mer occurrences inside arrays** for admissibility, as the chain heuristic does. Seeds from repeats are about quadratic. sporos MUMs are unique by construction, which avoids this, so keep repeat seeds out of any admissible-heuristic role.
3. **Don't let tie-breaking decide homology.** Add a canonical rule for co-optimal placements: left-normalize indels, or prefer gaps on unit boundaries consistent with the existing graph. This keeps copy-number events in the same place across haplotypes. POASTA's result depends on adjacency order (CLAUDE.md invariant).
4. **Don't use contiguous per-row windows when exploration is multi-diagonal.** A row touched at two far offsets allocates the span between them. In repeats, use sparse blocks or bands per diagonal cluster.
5. **Always cap work** with a state, time or memory budget, and define a fallback. POASTA has none.
6. **Don't assume per-step exactness gives a good MSA.** Order matters: add a guide order (for example length or similarity sorting, or anchor-consistent order) and optionally a re-alignment pass. POASTA has neither.
7. **Don't emit single-base GFA segments.** Unchop, and keep paths exact (#29).
8. **Don't rebuild global indices per sequence** if windows can be updated locally. With anchors, only the affected window's view and bubble data change.
9. **Pruning rules for gap states are subtle.** Two of POASTA's exactness bugs were in gap-piece pricing for pruning. Port the rules together with their regression cases (astar/mod.rs ~1301-1350) or skip pruning at first.
10. **Use wider costs than `u8`, and support a substitution matrix** (#24) if the evaluation scoring needs it.
