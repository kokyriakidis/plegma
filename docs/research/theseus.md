# Theseus (Marco-Sola Lab): source and paper review for plegma

- **Repo:** https://github.com/Marco-Sola-Lab/theseus-lib, cloned at HEAD `58440b4` (2026-06-19, "Fix github actions workflow file structure"). The upstream dev repo is `albertjimenezbl/theseus-lib`. The library is about 3.9k lines of C++23 (`theseus/`, `include/theseus/`), with two "not production ready" tools (`tools/theseus_msa.cpp`, `tools/theseus_aligner.cpp`). MIT licence.
- **Paper:** Jimenez-Blanco, Lopez-Villellas, Moure, Moreto, Marco-Sola, "Theseus: fast and optimal affine-gap sequence-to-graph alignment".
  - bioRxiv 10.64898/2026.02.12.705572 (Feb 2026).
  - Published in *Bioinformatics* 42(9) btag561 (PMC13544969).
  - Datasets: Zenodo 18482097 / 10.5281/zenodo.21132331.
- **Does it ship POA/MSA? Yes.** `TheseusMSA` builds a POA graph progressively, one sequence at a time. Outputs are MSA (FASTA), GFA, heaviest-bundle consensus, weighted majority-vote consensus, and DOT. It also has a sequence-to-graph mode that needs a user-given start vertex and offset.
- **Method:** I read the source only; I did not build or run it. Statements marked **[inferred]** are my reasoning, not something I read directly in the code or paper.

---

## 1. Graph representation and memory layout

### 1.1 Alignment graph: `theseus::Graph`
This is a pimpl wrapper around `GraphImpl` (`theseus/graph.cpp`).

- **Layout:**
  - Each node is an AoS record `struct Node { std::string sequence; std::vector<NodeId> in_nodes, out_nodes; bool alive; }`, held in `std::vector<Node> nodes_` (graph.cpp:363-370).
  - Each node therefore owns heap-allocated sequence bytes and two adjacency vectors. There is no CSR and no 2-bit packing.
  - Node IDs are `size_t`. Removed IDs are recycled through `free_node_ids_` (graph.cpp:80-100, 121-148).
  - Source and sink lists are maintained eagerly. `remove_id_from_vec` does a linear `std::find`, then swaps and pops (graph.cpp:150-163, 385-393).
- **Nodes are multi-character.** Typical nodes are long unitig-like segments, so WFA match runs (LCP) stay long inside a node.
- **Reverse views:** `node_rev(id)` returns a reversed `SequenceView` with in- and out-edges swapped (graph.cpp:289-300; graph.h:219-354). This lets the aligner run leftwards from a start point (`reverse_alignment`, impl.cpp:64-73). The reversal is not complemented. For GFA input, the tool adds a separate reverse-complement node for every segment (`name+`/`name-`) and both link orientations (tools/theseus_aligner.cpp:122-206). This was fixed after issues #1 and #2: the reverse was not complemented, and links could be crossed in one direction only.
- **Mutators used by POA:**
  - `expand_sequence` appends to a node (graph.cpp:102-108).
  - `split_sequence` splits a node (graph.cpp:110-119).
  - `substitute_out_edges` / `substitute_in_edges` move edges (graph.cpp:228-274).

### 1.2 Cycles
**Cyclic graphs are supported** in sequence-to-graph mode.

- The search is score-ordered and per-vertex. It never needs a topological order.
- A unit test aligns through a cycle: `tests/unit/seq_to_graph_tests.cpp`, edge `n4 -> n1`, expected paths `{n2, n4, n1}`.
- The paper claims the method "makes no assumptions about the graph topology". It benchmarks synthetic cyclic graphs against vg, which must unfold cycles.
- The POA graph built in MSA mode is a DAG by construction.

### 1.3 POA graph used in MSA mode (`theseus/msa.h`)
Theseus keeps **two parallel graphs**:

1. **A per-base POA graph** `POAGraph`.
   - Each `POAVertex` holds `char value`, `int weight`, and `std::vector<int> sequence_IDs, associated_vtxs (aligned nodes), in_edges, out_edges`, plus `associated_node_compact` (msa.h:50-76).
   - `_first_poa_vtx[c]` maps compact node `c` to its first POA vertex. The POA vertices of one compact node are contiguous IDs.
2. **A compacted alignment graph**, the `theseus::Graph` that the WFA kernel actually aligns against.

Construction:

- **Initial graph** (theseus_msa_aligner.cpp:37-53; msa.h:377-430): an empty source node 0, the first sequence as one big node 1, and an empty sink node 2. The first input sequence is the backbone. Sequences are then added in input order (tools/theseus_msa.cpp:224-235). There is no guide tree and no length sorting.
- **Adding an alignment** (`add_alignment_poa`, msa.h:290-370), after backtrace:
  - The compact-node path is converted to a POA-vertex path (`convert_path`, msa.h:260-286).
  - CIGAR ops are then walked:
    - `M`: increment the weight of the existing vertex and add the edge if it is new.
    - `X`: reuse the aligned sibling vertex with the same base if one exists, otherwise create a new aligned vertex (SPOA-style `associated_vtxs`) (msa.h:157-219).
    - `D` (query-only base): create a new vertex, appending to the compact node currently being grown via `expand_sequence` where possible (msa.h:324-350).
    - `I` (graph-only base): skip.
- **New edges split compact nodes.** When a new POA edge leaves or enters the middle of a compact node, `split_vertices` splits it (msa.h:86-143). Compact nodes only split and grow; nothing re-merges them. **[inferred]** After hundreds of divergent haplotypes, the compact graph fragments toward many short nodes, which erodes the LCP advantage.
- **Free-end shortcut:** the path end is always connected to the sink from the *end column* reached (msa.h:279-285, 358). With the semi-global end condition (§2.7), this adds a shortcut edge from mid-graph to the sink.

---

## 2. The algorithm

### 2.1 Scores and penalties
- **Scoring is a cost to minimise** (match ≤ mismatch, gap = o + L·e).
- **Tool defaults:** `m=0, x=2, o=3, e=1` (tools/theseus_msa.cpp:49-52). The paper does not state the benchmark penalties.
- **Non-zero match:** it is converted with the Eizenga transform. If match ≠ 0, then x' = 2x − 2m, o' = 2o, e' = 2e − m (internal_penalties.h:45-58).
- **Only single gap-affine works.**
  - `Penalties` declares Linear and DualAffine constructors (penalties.cpp:34-67).
  - The aligner has `TODO: Gap-linear and dual affine-gap` (impl.cpp:44).
  - I2/D2 wavefronts are allocated but never used (scope.h:205-214; vertices_data.h:79-91 commented out).
  - `gape2()` returns `gape_` (penalties.h:147, bug).
  - `InternalPenalties::_gapo2/_gape2` are never initialised.

### 2.2 How WFA is generalised to graphs: vertex-local diagonals
- **Each vertex v has its own local DP matrix** of query × `seq(v)`. A cell is `Cell { pos_t prev_pos; vertex_t vertex_id; int32 offset; int32 diag; Matrix from_matrix; }` (cell.h:65-96). That is 8+8+4+4+1 bytes, padded to **32 B per cell** **[inferred sizeof]**.
- **Coordinates:**
  - `offset` = query position i. This is a global coordinate, comparable across vertices.
  - `diag` = j − i, where j is the column inside the vertex (impl.cpp:253, 540, 599).
- **Entering a vertex:** a diagonal that reaches column j = |v| "jumps". Each out-neighbour w gets a new cell with `diag = -offset`, so j = 0 in w, with the same offset (impl.cpp:604-628 for M; impl.cpp:459-524 for M and I jumps out of I cells). Jump cells go into separate arrays `m_jumps_wf` / `i_jumps_wf`. Each vertex records the positions of its jump cells per score in a small circular scope (`_m_jumps_positions[score % nscores]`, vertices_data.h:84-91, 115-126).
- **Ops naming (careful):** Theseus `I` consumes a *graph* char (diag+1, offset+0). `D` consumes a *query* char (diag−1, offset+1). See impl.cpp:341 vs 389 and backtrace impl.cpp:719-730.
  - In the tests, a query longer than the graph yields trailing `D`s.
  - **[inferred]** This is the opposite of the SAM/GAF convention, yet it is written straight into the GAF `cg:Z` (impl.cpp:910-925).
  - Only I (graph-gap) needs cross-vertex jumps; a D never changes j.

### 2.3 Per-score step (`compute_new_wave`, impl.cpp:139-151)
For score s:

1. Expand and compact the invalid-diagonal segments of every active vertex (vertices_data.h:253-309).
2. For **every active vertex** (in activation order), run `process_vertex` (impl.cpp:120-136):
   - **`next_I`** (impl.cpp:330-374):
     - Sources: I[s−e] shifted (+1 diag), I-jumps at s−e, M[s−o−e], and M-jumps at s−o−e.
     - These are merged into the **scratchpad** (below) by max offset.
     - Survivors pass a validity and heuristics check and are copied into the score-s I wavefront.
     - Any I cell that reached j = |v| creates both an M-jump and an I-jump into each successor (`check_and_store_jumps`, impl.cpp:528-548).
   - **`next_D`** (impl.cpp:378-413): D[s−e] and M[s−o−e] (with M-jumps), shifted (offset+1, diag−1).
   - **`next_M`** (impl.cpp:417-455): M[s−x] and M-jumps at s−x (offset+1), plus I[s] and D[s] (no shift).
   - **Extend** every new M cell (impl.cpp:131-135).
3. Update the global heuristics, then s += 1 (impl.cpp:193-209).

**Scratchpad / "Unified Coordinate System"** (scratchpad.h:46-137; wavefront.h):

- One dense array of `Cell` is indexed by diagonal over [−|q|, max node length]. It is sized per alignment from the max node length (impl.cpp:85-96).
- `access_alloc(diag)` appends the diagonal to a touched list without a branch. It writes past `size` but inside capacity, then bumps the size if `offset == -1` (scratchpad.h:64-78).
- After each vertex and matrix, only the touched entries are reset (scratchpad.h:127-132).
- The merge is a branchless select `cell = (cell.offset < new.offset) ? new : cell` (impl.cpp:258-261).
- So each vertex's sparse cell list is merged into a dense diagonal array and written back as a compact list. The paper calls these steps densify, sparsify and extend; the code's naming is reversed.

**Scope:**

- **Circular buffer.** I/D wavefronts and per-vertex index ranges are kept for only `nscores = max(o+e, x) + 1` scores (impl.cpp:45-47; scope.h:65-93).
- **Kept forever.** The M and jump wavefronts are kept for the whole alignment in `BeyondScope`, for backtrace (beyond_scope.h:152-160). They are cleared only in `new_alignment`.

### 2.4 Extension across node boundaries
`extend_diagonal` (impl.cpp:580-631) uses an explicit stack (DFS):

1. Run a plain byte-by-byte LCP between the query and the node sequence (`LCP`, impl.cpp:552-565). There is no 64-bit block compare or SIMD.
2. Check the end condition.
3. If j hit the end of the node and the node has successors:
   - Invalidate that diagonal in the current vertex.
   - For every out-neighbour, activate it and check that the entry diagonal is still valid.
   - Push a jump cell and continue extending inside the successor.

Match runs therefore flow through arbitrarily many nodes in one score step. They fan out to every successor whose sequence matches. Vertices activated mid-score are extended immediately through the stack. Their `next_*` work starts at the next score, because the loop bound is captured at the start (impl.cpp:144-149).

### 2.5 Pruning that keeps exactness: diagonal invalidation
This is the core novelty, from vertices_data.h:54-92 and 318-418:

- **Why a diagonal can be invalidated.** When diagonal k of vertex v reaches the end of v at score s, it has already been propagated to every successor. Any later arrival on the same diagonal of v (score ≥ s) can only reach the same exit cell at a higher or equal cost, so it is dominated.
- **Neighbouring diagonals too.** Diagonal k±t could only reach the exit via a gap of cost o + t·e from the exited cell. So the invalid segment **grows one diagonal every e scores**, after an initial delay of o+e (M) or similar (I, D).
  - Separate segment lists are kept for M, I and D, each with its own rates (`invalidate_m_jump` vertices_data.h:351-375; `invalidate_i_jump` 318-342).
  - Each segment is `{start_d, end_d, rem_up, rem_down}`, decremented per score (`expand_invalid_vector`, 253-270).
  - Segments are sorted and merged each score (`compact_invalid_vector`, 192-245).
- **Use:** `valid_diagonal<M|I|D>(v, diag)` is checked before a cell is written. It is a **linear scan over that vertex's segments** (vertices_data.h:386-418).
- **Why it matters:** this stops cells re-entering the same (vertex, entry offset) and is what makes cycles terminate. The paper argues it does not lose optimality ("later exploration ... will be of a higher or equal score"). That is a proof sketch, not a formal proof.
- **[inferred] possible bug:** `VerticesData` is built with the *user* `Penalties` (impl.cpp:60; vertices_data.h:452), while all wavefront arithmetic uses `InternalPenalties`. With match ≠ 0 the Eizenga-scaled o and e differ from the user ones, so segment growth rates would be wrong. This could break exactness or pruning efficiency. It is harmless for m = 0, the default.

### 2.6 Heuristics (all optional; off by default)
In `include/theseus/heuristics.h`:

- **Lag pruning** (local, per cell; 132-138, 198-209):
  - Drop a cell if `max_offset − offset > 3·max(min(500, |q|/20), 50)`.
  - It applies only while the best offset advanced by more than that threshold within the last `lookback_lag` scores.
  - This is an X-drop measured in *query offset*, which is graph-agnostic.
- **Density drop** (global; 175-191):
  - After `s_min = max(min(1000, |q|/10), 100)·e` scores, abort if the max offset advanced less than `2·max(min(1000, |q|/10), 100)` over the last `s_min` scores.
  - It returns `END_UNREACHABLE`, and a partial backtrace from the best-offset cell is produced (impl.cpp:221-230, 636-664).
- **No max-score or step cap.** `THESEUS_STATUS_MAX_STEPS_REACHED` is defined but unused (alignment.h:54; heuristics.h:218-221 commented out).
- **Wrapper bug** (impl.h:81-91 vs theseus_aligner.cpp:78): `TheseusAligner::align` passes `(seq, start, off, 1, density_drop, lag_pruning, false, false)` to `align(seq, start_node, start_offset, weight, is_ends_free, reverse_alignment, density_drop, lag_pruning)`.
  - In seq-to-graph mode, `-d` sets the unused `is_ends_free` and `-l` turns on **reverse alignment**.
  - The heuristics are never enabled.
  - The MSA wrapper passes the arguments correctly (theseus_msa_aligner.cpp:75).

### 2.7 Exactness, alignment modes and end condition
- **The first hit is optimal.** The algorithm terminates at the first score s at which any cell reaches `offset == |q|` (`check_end_condition`, impl.cpp:569-576). Because scores are processed in increasing order, that first hit is optimal for the problem actually solved.
- **The problem actually solved: start fixed, query end-to-end, graph end free (semi-global on the graph end).**
  - The end vertex is *not* checked.
  - `_end_node` is set in MSA mode but never read (impl.cpp:173, 178).
  - `is_ends_free` / `_ends_free` is stored but never used (impl.cpp:187).
  - In MSA mode the start is the empty source node 0 at offset 0 (impl.cpp:169-180), so skipping a graph *prefix* costs gaps, but skipping a graph *suffix* is free.
  - The README states this design intent: "allowing to add partial sequences, as long as they start on either end of a backbone sequence".
  - The author's reply to issue #3 calls this "global alignment: given an initial position in the graph, it will align the full query".
- **No local mode, no query-ends-free mode, no true global (must-reach-sink) mode.**

### 2.8 Memory and traceback
- **Only M and jump cells are retained.** Every M cell, M-jump cell and I-jump cell of the whole alignment stays in `BeyondScope`, at 32 B each.
- **I and D cells are pass-through.** They live only in the circular scope. When an I or D cell feeds an M cell, it is copied whole, so the M cell inherits the I/D cell's `prev_pos` and `from_matrix`, which point back to the M cell where the gap opened (impl.cpp:300-327, 425-433).
- **Backtrace** (impl.cpp:700-765):
  - Follow `prev_pos` pointers.
  - Inside a vertex, a diagonal difference between consecutive M cells is a gap of that length (I or D), plus the matches between them. An equal diagonal means mismatch plus matches.
  - On a vertex change (jump), emit matches and then `|v_prev| − j_prev` I's for the part of the previous node that a gap skipped.
- **Memory class:** this is the "Singletrack" idea. The paper says Singletrack cuts memory 3×. The cost is still **O(total M cells) memory over the whole alignment**: no BiWFA and no O(s) memory mode.
- **Per-vertex overhead [inferred]:**
  - Each activated vertex keeps `2·nscores` `std::vector`s of jump positions plus three segment vectors.
  - Vertices are **never deactivated** within an alignment (`activate_vertex`, vertices_data.h:435-449; cleared only in `new_alignment`).
  - So every score iterates over *all vertices ever touched*, and calls `next_I/D/M` and expand/compact on each, even if they hold no live cells (impl.cpp:139-151).

---

## 3. SIMD, bit-parallelism, parallelism
- **None.** I grepped for intrinsics, OpenMP and threads and found nothing. The LCP is a scalar char loop (impl.cpp:561-564). The only micro-optimisations are the branchless max-merge and touched-list append in the scratchpad, plus a custom `Vector` that skips initialisation for trivial types (vector.h; cell.h:63-66).
- The paper's Discussion confirms the implementation is single-threaded, makes no heavy use of SIMD and has no intra-alignment parallelism, and lists these as future work.
- **[inferred]** Aligner objects own all their state, so running independent windows on separate threads with separate aligner instances would be straightforward. Nothing ships that does it.

## 4. Seeding / banding
- **None.** Sequence-to-graph mode requires the caller to provide `(start_vertex, start_offset, orientation)`. The tool reads them from FASTA headers `> v1 3 +` (tools/theseus_aligner.cpp:210-271). In effect, the caller supplies the anchor.
- There is no band. The only search-space limiters are WFA score ordering, diagonal invalidation, and the optional offset-based heuristics (§2.6). The paper positions Theseus as an alignment *kernel* behind an external mapper, using inputs extracted from vg map runs.

## 5. MSA / consensus output (msa.h)
- **MSA FASTA** (`poa_to_fasta`, msa.h:693-716; `poa_fasta_and_majority`, 452-685):
  - Build an augmented graph, adding edges between all aligned-node pairs of each edge's endpoints (O(E·k²)).
  - Topologically sort it with an iterative DFS.
  - Assign columns: a vertex takes the column of an already-placed aligned vertex, otherwise a new column.
  - Fill rows from each vertex's `sequence_IDs`.
  - Rows are named `Sequence_1..N` (no original names) plus a `Consensus` row. Source and sink columns are dropped.
- **Weighted majority voting:**
  - Coverage-aware: each sequence's start and end columns give the number of "valid" sequences per column, which handles partial sequences.
  - A column is called if its non-gap weight is more than half the valid weight; the base is the heaviest vertex (msa.h:620-684).
  - Outputs are the consensus, the gapped consensus and per-column weights.
- **"Heaviest bundle" consensus** (msa.h:723-761) is a *greedy walk*: from the source, step to the successor with the largest vertex weight. It is not Lee's DP heaviest bundle (SPOA/abPOA).
- **Per-sequence weights** are supported: `align(seq, weight)` adds `weight` to every vertex on the path (msa.h:310, 329). The tool hard-codes weight 1; header-weight parsing is commented out (tools/theseus_msa.cpp:96-100).
- **GFA output** prints numeric node IDs with `0M` links (impl.cpp:806-832). Issues #4 and #5 report empty segment names in v0.1, fixed upstream.

## 6. Behaviour on tandem repeats / VNTRs
**Exactness:**

- **Greedy extension is exact under WFA's rules.** With match cost 0 after the transform, taking matches greedily along a diagonal is optimal; this is standard WFA dominance. Running into "the wrong copy" is therefore not an *optimality* problem: Theseus still returns a minimum-cost alignment for its objective.
- **[inferred]** The real problems for VNTRs are four.
  1. **Many co-optimal placements.** A copy-number indel can sit at any copy boundary. The backtrace picks whichever source won the strict `<` max-merge in the scratchpad. The insertion order in the scratchpad is I/D-from-gap, then jumps, then M, so ties keep the first source written (impl.cpp:260). Gap placement is therefore deterministic but not normalised: neither left-aligned nor period-aligned. In progressive POA this tends to put the same CN indel in different columns for different haplotypes.
  2. **The free graph end is wrong for haplotypes** (§2.7). Suppose the graph is `X R^10 Y` and the query is `X R^5 Y` with a short Y. Matching `X R^5`, mismatching Y against the start of the 6th copy and then stopping for free costs about x·|Y|·(mismatch rate). The correct 5-copy gap costs o + 5|R|·e. Whenever |Y| is short, as with a MUM-bounded window whose right anchor is only tens of bp, the free-end alignment wins and is **wrong for global haplotype MSA**. The convert step then also adds a mid-graph edge to the sink (msa.h:281-285).
  3. **The heuristics are dangerous here.**
     - During a graph-only gap (extra copies in the graph), the query offset does not advance for G·e scores. Density drop would abort once `s_min` (≥ 100·e) scores pass with fewer than about 200 offsets of progress.
     - Lag pruning drops cells that lag in offset, which is exactly the correct "skip copies" path while wrong-copy paths race ahead.
     - Use exact mode only.
  4. **No dual-affine,** so long CN gaps are charged o + L·e with a single e.

**Cost scaling [inferred, order-of-magnitude]:**

- Work and memory follow WFA's O(n·s) time with an O(s²)-cells flavour, multiplied over vertices.
- A copy-number difference of Δ bp gives s ≥ o + Δ·e. With defaults (e = 1) a 3 kb difference means s ≈ 3,000. That is about s/e live diagonals by the end, and on the order of s²/(2e) ≈ 4.5 M M-cells × 32 B ≈ 150 MB retained for backtrace per alignment, before vertex and segment overheads.
- Inside an exact tandem array, every diagonal at a multiple of the period p extends for long runs. Extension work is O(n) per diagonal (offsets are monotone), so O(n·s/p) to O(n·s) character comparisons.
- After the graph fragments, the per-score loop over all activated vertices and the linear `valid_diagonal` scans add an O(s × V_active) term.
- The paper concedes this scaling: runtime and memory "scale with the optimal alignment score", and the advantage shrinks with divergence. On HIV-1, the most divergent set, Theseus is slower than SPOA (4.2×) and abPOA (9.2×) (Table 1 below).

**GitHub issues:** five, all closed, and none about repeats or performance.

- #1: links traversed one way only.
- #2: impossible CIGARs from a reverse-complement bug.
- #3, filed by kokyriakidis: per-read constraints for semi-global alignment. The answer was that only "global" from a given start is implemented.
- #4 and #5: empty GFA names.
- PR #6 (open) adds penalty tests.

## 7. Performance numbers from the paper
**MSA, journal version (PMC Table 1, seconds), single thread, Xeon Platinum 8480, 256 GB:**

| Dataset | SPOA | abPOA | POASTA | Theseus |
|---|---|---|---|---|
| HIV-1, 10 kb × 100 | 26 | **12** | 174 | 110 |
| SARS-CoV-2, ~30 kb × 2,732 | 7,664 | 205 | 65 | **33** |
| M. tuberculosis 50 kb × 342 | 2,571 | **124** | 6,512 | 370 |
| M. tuberculosis 250 kb | oom | oom | 8,690 | **910** |
| M. tuberculosis 500 kb | oom | oom | 29,024 | **2,453** |
| M. tuberculosis 1 Mb | oom | oom | 71,770 | **5,996** |
| Monkeypox ~200 kb × 100 | oom | oom | >48 h | **7,785** |

**Peak memory (GB):**

| Dataset | SPOA | abPOA | POASTA | Theseus |
|---|---|---|---|---|
| HIV-1 | 4.5 | 7.1 | 0.4 | 0.5 |
| SARS-CoV-2 | 42.6 | 13.5 | 2.4 | 0.8 |
| Mtb 50 kb | 116.9 | 66.0 | 26.5 | 53.1 |
| Mtb 250 kb | oom | oom | 43.9 | 52.2 |
| Mtb 500 kb | oom | oom | 46.3 | 52.7 |
| Mtb 1 Mb | oom | oom | 59.0 | 54.4 |
| Monkeypox | oom | oom | timeout | 204.9 |

**Headline claims and other results:**

- **MSA:** 1.6–17.6× faster than POASTA. The bioRxiv version says 2.0–12.0× and omits HIV-1. It is 7.3× faster than SPOA on average, but 4.2× *slower* on HIV-1. abPOA is faster on HIV-1 (9.2×) and Mtb 50 kb (3.0×).
- **Memory:** Theseus uses about 2× POASTA on Mtb 50 kb (53 vs 26.5 GB), and 205 GB on Monkeypox.
- **Accuracy:** TC and SP scores are shown in figures only. Theseus is "consistently among the most accurate". All tools except abPOA are optimal.
- **Read mapping kernel** (100/150/250 bp reads, GRCh37 + 1000G AF ≥ 1% subgraphs): 1.9–16.9× faster than the vg map alignment stage, SPOA, abPOA and POASTA. The numbers are in figures only.
- **Synthetic cyclic graphs:** 1.5–36.3× faster than vg, which unfolds cycles.
- **Not benchmarked:** GraphAligner, GWFA (Zhang et al., arXiv 2206.13574, edit distance only), Astarix, minigraph. Benchmark penalties and sequence order are not stated.
- **Relevance gap:** the datasets are low-divergence microbial and viral genomes, not tandem-repeat-rich human loci.

---

## 8. Advantages and disadvantages for plegma

### Advantages
1. **Exact affine-gap graph alignment whose cost scales with s, not with |G|·|q|.** Haplotypes identical or near-identical to an existing path cost little more than one LCP pass, and HPRC VNTR windows are dominated by such haplotypes.
2. **Works on multi-character (compacted) nodes,** with match runs crossing node boundaries through DFS. This is a natural fit for a compacted POA graph or for MC graph nodes.
3. **Handles cycles natively** through diagonal invalidation, with no unrolling. This could matter if plegma ever aligns to cyclic (collapsed-repeat) graphs.
4. **Progressive POA is already built in:** aligned-node reuse on mismatch, per-sequence weights, MSA FASTA, coverage-aware majority consensus, GFA.
5. **Small, readable implementation** (about 3.9k LOC). The kernel ideas (UCS scratchpad, pass-through gaps, jump cells) port easily.
6. **Supports reverse alignment** (leftward from a start), which is useful for extending both ways out of an anchor.
7. **Query-offset-based heuristics** (lag pruning and density drop) are graph-agnostic X-drop analogues.

### Disadvantages
1. **There is no global end condition.** The graph suffix is always free (impl.cpp:569-576) and `is_ends_free` is dead code. For haplotype MSA inside anchored windows this yields wrong alignments (§6.2). Only the start is anchored.
2. **There is no dual-affine or convex gap model.** VNTR copy-number indels inflate s linearly, and WFA cost grows about quadratically in s.
3. **Memory is O(total M cells) at 32 B per cell** with no BiWFA or O(s)-memory mode. The paper itself reports 53–205 GB on large MSAs.
4. **Per-score cost includes every vertex ever activated,** with linear invalid-segment scans. It degrades on fragmented POA graphs, which is exactly what hundreds of VNTR haplotypes produce.
5. **No SIMD, no 64-bit block LCP, no threads.** SPOA and abPOA beat it on the most divergent dataset (HIV-1).
6. **No seeding, banding or score cap:** a pathological window has no bailout except the heuristics, which are unsafe for VNTRs.
7. **Gap placement in repeats is arbitrary.** No normalisation means inconsistent MSA columns across haplotypes.
8. **POA order is input order with the first sequence as backbone,** with no guide tree. The greedy "heaviest bundle" is not the standard DP. MSA rows lose sequence names.
9. **Code-quality bugs to note:**
   - Wrapper argument order (theseus_aligner.cpp:78).
   - `gape2()` (penalties.h:147).
   - VerticesData uses user rather than internal penalties (impl.cpp:60).
   - I/D inverted relative to SAM in the GAF `cg` (impl.cpp:910-925).
   - GAF target length is summed from *node-name* lengths (impl.cpp:882-885).
   - `store_I_jump` recursion tests the *current* node's emptiness instead of the destination's (impl.cpp:519-520).
   - `_i_jumps_wf_pos` records the size of the unused I2 array (beyond_scope.h:89).
   - Unit tests pass `false` as a weight (msa_tests.cpp).

---

## 9. Ideas to steal and pitfalls to avoid for plegma

### Ideas to steal (concrete)
1. **A graph-WFA kernel for MUM/SMEM-bounded windows, made truly global.**
   - Start: the left anchor's end vertex and column, with query offset 0.
   - Accept: only a cell on the right anchor's start vertex/column with `offset == |q_window|`. That is one extra comparison in the end check.
   - Inside windows where s is small (SNPs, short indels, ±1 copy), this is far cheaper than full POA DP and is exact.
2. **Dispatch on a cheap lower bound for s, with a score cap.**
   - Per window, compute the min and max path length [Lmin, Lmax] of the graph sub-DAG between the anchors in O(V+E).
   - Then s ≥ o + e·dist(|q|, [Lmin, Lmax]) whenever |q| lies outside that range, with o charged only for a non-zero distance.
   - If the bound, or a hard `s_max`, is exceeded, fall back to a SIMD or banded DP kernel (abPOA/SPOA-style).
   - Theseus has no such cap; plegma needs one.
3. **Vertex-local diagonals with a global query-offset coordinate.** Use `diag = j − i` per vertex, enter successors at `diag = −offset`, and make every heuristic key on query offset. This is simple and graph-agnostic.
4. **A dense UCS scratchpad with a touched list** to merge sparse per-vertex cell lists by max offset, using a branchless select and resetting only touched entries. It maps well to Mojo `SIMD`/`UnsafePointer` buffers.
5. **Pass-through gap cells for backtrace.** Keep only M cells plus I-jump cells. Gap cells inherit the origin M cell's back-pointer, and gap length comes from the diagonal difference. Shrink the cell to about 16 B: 32-bit prev index, 32-bit vertex index, 32-bit offset, 16-bit diagonal or vertex-local column, and a 2-bit matrix tag. **[inferred]** This roughly halves Theseus's footprint.
6. **Diagonal invalidation as dominance pruning,** if plegma ever aligns to cyclic graphs or wants to avoid re-entering (vertex, entry offset). For a DAG window, a per-vertex bitset or interval set of "exited at offset" may be enough; keep segments sorted and binary-search them instead of Theseus's linear scan.
7. **Two graphs:** a per-base POA graph for MSA bookkeeping (aligned-node sets, weights, sequence IDs) and a compacted graph for alignment. Unlike Theseus, **re-compact** unbranched chains periodically so LCP runs stay long.
8. **Weights for deduplicated haplotypes.** In HPRC windows many haplotypes are identical. Align each distinct window sequence once and pass `weight = multiplicity` (Theseus's `align(seq, weight)` model). Keep the ID lists so rows can be expanded back to every haplotype on output.
9. **Coverage-aware majority consensus** (per-column valid weight from each sequence's start and end column) as a robust consensus option alongside a proper heaviest-bundle DP.
10. **Bidirectional extension from an anchor** using reversed node views (Theseus `node_rev`), as a step toward a BiWFA-style O(s)-memory meet-in-the-middle kernel within a window.

### Pitfalls to avoid
1. **Do not inherit the free-graph-end semantics.** Implement explicit end modes (global, query-ends-free, graph-ends-free) and test with short right flanks inside repeat arrays.
2. **Do not use offset-based X-drop or density-drop heuristics inside VNTR windows.** Graph-only gaps (missing copies) stall query-offset progress, and wrong-copy paths race ahead.
3. **Do not use single-affine gaps for copy-number events.** Support two-piece (dual) affine gaps (minimap2-like o2/e2), or charge CN indels per copy, so that s, and therefore WFA cost, stays small. Theseus declares dual affine but never implemented it.
4. **Expect no MUMs inside tandem arrays,** because copies are not unique. Windows spanning a whole array will have large s. Pair the WFA kernel with a fallback kernel and maybe period-aware anchoring (SMEMs at array boundaries) rather than relying on WFA alone.
5. **Normalise gap placement in repeats** with a deterministic, period-aware left- or right-shift, so successive haplotypes hit the same POA columns. Theseus's arbitrary tie-breaking scatters CN indels.
6. **Do not keep every activated vertex in the per-score loop.** Retire vertices with no live cells in the scope window and no pending jumps.
7. **Do not compare byte by byte.** Use 2-bit-packed or 64-bit-word XOR plus count-trailing-zeros LCP, and handle node boundaries by word-level compare up to the remaining node length.
8. **Avoid an O(all cells) backtrace store for large s.** Use BiWFA-style recursion, or at least cap memory and fall back.
9. **Choose the progressive order deliberately** (guide tree or longest/most-central first) instead of input order with the first sequence as backbone. Keep sequence names in the MSA output.
10. **Testing lessons from Theseus's bugs:**
    - Test non-zero match (the Eizenga transform paths).
    - Keep CIGAR I/D conventions consistent with SAM/GAF.
    - Fuzz wrapper argument orders.
    - Cross-check every alignment's recomputed score against the claimed score; Theseus ships `compute_affine_gap_score` (alignment.h:104-140) for this.

---

### Sources
- Code: https://github.com/Marco-Sola-Lab/theseus-lib (HEAD 58440b4). Issues #1–#5 and PR #6 were read with `gh api`.
- Paper (journal): https://academic.oup.com/bioinformatics/article/42/9/btag561/8742314 and https://pmc.ncbi.nlm.nih.gov/articles/PMC13544969/
- Preprint: https://www.biorxiv.org/content/10.64898/2026.02.12.705572v1.full
- Related: GWFA, https://arxiv.org/abs/2206.13574
