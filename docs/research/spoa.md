# SPOA (SIMD Partial Order Alignment): source review for plegma

- **Target:** https://github.com/rvaser/spoa, shallow clone at commit `c987896` (2026-10-16), version 4.1.6 (`meson.build:3`). MIT license. About 4.3k lines of C++17.
- **Clone location:** `/private/tmp/claude-501/-Users-kokyriakidis/8500d74c-03af-4a3f-8699-f1899147a1f9/scratchpad/src/spoa`. All `file:line` references below are relative to that directory unless another repo is named.
- **Downstream code read:** racon (lbcb-sci/racon, shallow clone in `.../scratchpad/src/racon`), smoothxg and pggb, and the Cactus clone. smoothxg, pggb and Cactus were already in `.../scratchpad/src/`, cloned by another task. I only read them.
- **Method:** I read the source and did not build or run it. Statements marked **[inferred]** are my reasoning from the code, not something stated in it.

Core files:

| File | Lines | Role |
|---|---|---|
| `include/spoa/graph.hpp` | 324 | Graph, Node and Edge API |
| `src/graph.cpp` | 694 | Graph construction, topological sort, consensus, MSA, subgraph, DOT |
| `include/spoa/alignment_engine.hpp` | 98 | Engine API, modes and gap models |
| `src/alignment_engine.cpp` | 113 | Factory: chooses gap model, SIMD first, falls back to scalar |
| `src/sisd_alignment_engine.cpp` | 928 | Scalar DP engine |
| `src/simd_alignment_engine_implementation.hpp` | 2065 | SIMD DP engine (SSE4.1/AVX2/SIMDe) |
| `src/simd_alignment_engine_dispatcher.cpp` | 119 | Runtime CPU dispatch |
| `src/main.cpp` | 359 | CLI, including GFA output |

---

## 1. Graph representation

### Data structures (`include/spoa/graph.hpp:40-100`)

- **`Node`** (`graph.hpp:54-58`) has these fields:
  - `id` (uint32): the creation index, which is also the index into `nodes_`.
  - `code` (uint32): an alphabet code.
  - `inedges` and `outedges`: `std::vector<Edge*>`.
  - `aligned_nodes`: `std::vector<Node*>`.
- **`Edge`** (`graph.hpp:84-87`) has these fields:
  - `tail` and `head` (`Node*`).
  - `labels`: `std::vector<uint32_t>`, the IDs of the sequences that traverse this edge.
  - `weight` (int64).
- **`Graph`** (`graph.hpp:312-319`) holds:
  - `num_codes_`, `coder_[256]` and `decoder_[256]` (the alphabet is learned dynamically).
  - `sequences_`: the first node of each sequence.
  - `nodes_` and `edges_`: each a `vector<unique_ptr<...>>`.
  - `rank_to_node_`: the topological order.
  - `consensus_`.
- **Memory layout:** pointer-heavy, an array of structs. Every node is a separate heap allocation, and each node holds three vectors, each another heap allocation. Edges are heap objects too. There is no compact CSR layout.
- **Alphabet:** `coder_` assigns codes in first-seen order (`graph.cpp:169-174`). Any byte is accepted, so `N` or lowercase letters become extra codes, and every new code adds a row to the query profile.
- **Copying and threading:** the graph cannot be copied, only moved (`graph.hpp:29-33`). It is not thread-safe. Cereal serialization is optional (`graph.hpp:193-283`).

### Sequence labels on edges

- `AddEdge` (`graph.cpp:82-92`) scans `tail->outedges` for an existing edge to the same head. If it finds one, it appends the new sequence ID to `labels` and adds to `weight` (`Edge::AddSequence`, `graph.cpp:61-64`). Otherwise it creates a new edge.
- Edge weight is `weights[i-1] + weights[i]`: both endpoint bases contribute (`graph.cpp:106`, `:238`). With FASTQ input, each base weight is its Phred score (`graph.cpp:143-145`). Otherwise the per-sequence weight defaults to 1 (`graph.cpp:124`).
- Paths are not stored explicitly. A sequence is recovered by starting at `sequences_[i]` and repeatedly calling `Node::Successor(i)` (`graph.cpp:23-31`).
  - `Successor` does a linear `std::find` over the `labels` vector of each out-edge.
  - Cost of one step: O(out-degree × labels per edge).
  - **[inferred]** With 456 haplotypes, an edge can carry 456 labels, so emitting all MSA rows costs roughly Σ path length × 456 comparisons. That is about 6×10⁹ for 456 × 30 kb. This is a real hidden cost.
- `Node::Coverage()` (`graph.cpp:33-48`) builds an `unordered_set` of all labels on in-edges and out-edges for every call.

### Aligned-node sets (columns)

- `aligned_nodes` is a clique of nodes that occupy the same alignment column but carry different bases.
- In `AddAlignment` (`graph.cpp:213-232`), when the sequence base is aligned to node `jt`:
  - If `jt->code` equals the base, `jt` is reused.
  - Otherwise, if some node in `jt->aligned_nodes` carries the base, that node is reused.
  - Otherwise a new node is created and linked bidirectionally to `jt` and to every member of `jt`'s aligned set.
- Each column therefore holds at most `num_codes` nodes, one per base.
- Insertions (`it.first == -1`) always create new nodes (`graph.cpp:210-211`). They are never merged with existing parallel insertion nodes unless the DP itself aligned the base to them.

### Adding a sequence (`Graph::AddAlignment`, `graph.cpp:156-248`)

- An `Alignment` is a `vector<pair<int32,int32>>` of `(node_id or -1, seq_pos or -1)` (`graph.hpp:23`).
- With an empty alignment, the whole sequence is added as a new linear chain (`graph.cpp:176-180`). This is how the first sequence becomes the backbone, and it is also how racon inserts its backbone.
- The unaligned prefix and suffix (outside `valid.front()..valid.back()`) are added as new unmerged chains (`graph.cpp:198-200`). The suffix is connected at `graph.cpp:242-244`. In local (SW) and overlap modes, the unaligned flanks become dangling branches.
- `TopologicalSort()` runs again after every addition (`graph.cpp:178`, `:247`). The sort is a full O(V+E) pass with no incremental update.

### Topological sort (`graph.cpp:250-304`)

- The sort is an iterative DFS with an explicit stack and three-color marks.
- Before a node is emitted, the DFS requires all predecessors of the node **and of every node in its aligned set** to be emitted.
- The node's aligned nodes are then emitted **contiguously** right after it (`graph.cpp:284-291`). The `ignored[]` flags stop aligned nodes from being emitted twice.
- This contiguity is what makes MSA column assignment trivial: `InitializeMultipleSequenceAlignment` (`graph.cpp:322-338`) walks `rank_to_node_` and gives the same column index to a node and its aligned set.
- DAG and sortedness checks are `assert`s only (`graph.cpp:282`, `:303`), and they are disabled in release builds (`b_ndebug=if-release`, `meson.build:9`).

---

## 2. Dynamic programming

### Scoring models (`alignment_engine.hpp:16-26`, `alignment_engine.cpp:33-74`)

- **Scores:** match `m`, mismatch `n`, and gap parameters `g`, `e`, `q`, `c`, all `int8_t`. Gap parameters must be ≤ 0 (`alignment_engine.cpp:47-56`).
- **Gap-length convention:** a gap of length `i` costs `g + (i-1)·e`, so a 1-base gap costs `g`. minimap2 charges O+E for a 1-base gap; SPOA charges only `g`. **[inferred]** Parameters copied from minimap2, as smoothxg does, are therefore slightly cheaper in SPOA.
- **Model selection** (`alignment_engine.cpp:58-67`):
  - `g >= e` gives **linear** gaps (with `e := g`).
  - Otherwise, `g <= q || e >= c` gives **affine** gaps.
  - Otherwise the model is **convex**: two-piece affine, `max` of the two affine functions, with matrices O and Q added.
  - CLI defaults are `5, -4, -8, -6, -10, -4`, which select convex (`main.cpp:208-213`).
- **Matrix count:** linear uses H only. Affine uses H, F (vertical, a deletion from the sequence) and E (horizontal, an insertion). Convex adds O and Q, for 5 matrices.

### Modes

There is an enum `kSW`, `kNW`, `kOV` (`alignment_engine.hpp:16-20`).

- **SW** (local, Smith-Waterman): the first row and column are 0 and every cell is clamped at 0 (`sisd_alignment_engine.cpp:185-192`, `:353-355`). The best cell anywhere in the matrix is the endpoint.
- **NW** (global, Needleman-Wunsch): the end cell must be a **sink node** (no out-edges) at the last sequence column (`sisd...cpp:356-358`).
- **OV** ("semi-global"): the mode is **asymmetric**.
  - The first column is 0, so leading graph nodes are free (`sisd...cpp:250-252`; SIMD `:633-635`).
  - Row 0 carries gap penalties, so the sequence prefix is **not** free (`sisd...cpp:233-246`).
  - The end can be at any sequence position, but only on a sink node (`sisd...cpp:359-360`).
  - The result is a dovetail: the start of the sequence overlaps the end of the graph. Trailing graph nodes are not free, and a sequence that starts before the graph is not handled symmetrically.
  - **[inferred]** For plegma's sub-window alignment between two anchors, only NW, or a properly symmetric free-end variant, is appropriate.

### DP over a DAG

- Rows are graph nodes in topological rank order, and columns are sequence positions. The matrix row index is `rank+1`, and row 0 is a virtual source (`sisd...cpp:322`).
- Each row combines **all predecessors** given by its in-edges:
  1. **First pass, first predecessor (`inedges[0]`):** fill the diagonal (match) term and the vertical term (F in affine, H+g in linear) from that predecessor's row (`sisd...cpp:330-334`).
  2. **Further passes, one per extra predecessor:** each pass takes an elementwise max into the same row (`sisd...cpp:336-348`). In-degree therefore multiplies the cost of the row.
  3. **Final left-to-right pass:** compute the horizontal term (E) and the final H (`sisd...cpp:350-362` for linear, `:528-543` for affine).
- First-column values (the cost of a graph prefix of gaps) are propagated through predecessors as the maximum over in-edges (`sisd...cpp:165-174`).

### SIMD engine (`simd_alignment_engine_implementation.hpp`)

**Dispatch and precision**

- **Instruction sets:** AVX2 (256-bit) or SSE4.1 (128-bit), selected at compile time (`:35`, `:130`). Through SIMDe the code can be ported elsewhere, including NEON, which is relevant to Apple Silicon (`:16-26`, `meson_options.txt`).
- **Runtime dispatch:** optional, through cpu_features or `cpuid` (`simd_alignment_engine_dispatcher.cpp:5-117`).
- **Fallback:** if no SIMD support was compiled in, `Create` returns null and the scalar engine is used (`alignment_engine.cpp:69-72`).
- **Precision is chosen per alignment from a worst-case score bound:**
  - `WorstCaseAlignmentScore(seq_len+8, graph_nodes)` (`alignment_engine.cpp:102-111`).
  - The engine uses int16 if the bound is above −32768+1024, int32 if it is above the int32 limit, and otherwise throws "possible overflow" (`simd...hpp:664-719`).
  - No int8 path exists. The Racon paper says 8-bit precision is insufficient.
  - **[inferred]** With CLI defaults (e = −6), int16 is only possible when `seq_len + graph_nodes` is at most about 5.3k. With smoothxg's parameters (e = −2, c = −1) the limit is about 15k.
  - **[inferred]** For plegma's 10-30 kb inputs SPOA always falls back to int32, which halves the lanes: 4 on SSE, 8 on AVX2.

**Vectorized dimension and layout**

- The **sequence (query) dimension** is vectorized, with each vector holding `kNumVar` **consecutive** sequence positions (`:479-486`).
  - This is the Rognes-Seeberg-style contiguous layout that the Racon paper cites. It is **not** Farrar's striped layout.
- Query profile: one vector row per alphabet code, `sequence_profile[code * W + j]` (`:477-487`).
  - Padding lanes past the sequence end get `-max(|m|,|n|,|g|,|q|)` (`:471-473`, `:481-482`).
- Column 0 is not stored in the vector matrix. It is kept in a separate scalar array, `first_column` (one per H/F/O), at `:419-447` and `:513-544`.

**Diagonal and vertical terms**

- The diagonal term shifts the predecessor's vector left by one lane (`_mmxxx_slli_si(H_pred_row[j], kLSS)`). The predecessor's last lane from the previous vector, `x`, is ORed in (`:779-795`).
- The vertical term is a plain lane-wise add of `g`.
- Each additional predecessor runs a second full sweep over the row, taking an elementwise max (`:797-821` linear; `:1154-1184` affine; convex at `:1627-1660`). Predecessor rows are combined by max, one full row pass per in-edge.

**Horizontal term: in-register prefix max**

- The E/insertion dependency within a vector is resolved with a log₂(k)-step **Hillis-Steele prefix max**, `_mmxxx_prefix_max` (`:84-92`, `:120-127`, `:178-185`, `:213-219`).
  - Each step is `a = max(a, (a + penalty·2^s) << 2^s lanes)`.
  - Precomputed `masks` (−∞ in the low lanes) and doubled `penalties` are built at `:743-759` and `:1100-1116`.
  - The carry across vectors is the last lane of the previous vector plus `g`, through `x` (`:824-840`).
- In the affine case, E is first seeded as `shift(H)+g−e+e` and then prefix-maxed with `e`; x carries `max(H, E−(g−e))` (`:1190-1206`).
- The convex case runs two prefix-max passes, one for E and one for Q (`:1665-1700`).
- The README says AVX2 is "marginally faster due to high latency shifts". A cross-lane shift on AVX2 takes `permute2x128 + alignr` (`:52-57`).

**Endpoint search**

- For each row, a running max vector is reduced with `_mmxxx_max_value` (`:241-250`). For SW and OV, `_mmxxx_index_of` then finds the column (`:261-276`, `:882-898`).
- **[inferred, possible bug]** `_mmxxx_max_value` starts from 0 (`:242`), and the affine and convex `score` accumulators also start at `zeroes` (`:1187`, the convex `score = zeroes`). In OV mode, a best score in a sink row that is negative is therefore reported as 0. The following `index_of(row, 0)` may then return −1, which makes the traceback empty (`:930`).
  - An empty alignment makes `AddAlignment` insert the sequence as a disconnected chain (`graph.cpp:176-180`).
  - SW clamps at 0 anyway, and NW uses `_mmxxx_value_at`, so this only affects OV with negative scores.

### Scalar engine (`src/sisd_alignment_engine.cpp`)

- The scalar engine has the same algorithm with int32 everywhere (`:13-14`, `:29-49`).
- Memory is the same: 1, 3 or 5 full int32 matrices of (nodes+1) × (seq_len+1) (`:92-115`).

### Memory: full matrices, always

- **SIMD:** `M_storage` holds `{1,3,5} × (nodes+1) × ceil(L/k)` vectors (`simd...hpp:418-447`). The scalar engine is the same at `sisd...cpp:92-115`. Full-matrix storage is required because the traceback recomputes moves from stored scores; it keeps no traceback-pointer bits.
- **Buffers only grow:** `Realloc` grows buffers and never shrinks them (`:414-447`). **[inferred]** An engine object keeps its peak footprint. Racon creates one engine per thread (`racon/src/polisher.cpp:179-183`), and pggb's help text tells users to "set lower [poa threads] if you OOM during smoothing" (`pggb/pggb:250`).
- **`Prealloc(max_len, alphabet)`** estimates the graph height as `max_len × alphabet` (`simd...hpp:381-394`, `sisd...cpp:72-75`).
- **[inferred] Footprint at int32** (bytes ≈ matrices × 4 × L × N, where N is the number of graph nodes):

| L (seq) | N (graph nodes) | linear | affine | convex (CLI default) |
|---|---|---|---|---|
| 500 (racon window) | 2k | 4 MB (int16: 2 MB) | 12 MB | 20 MB |
| 10 kb | 20k | 0.8 GB | 2.4 GB | 4.0 GB |
| 30 kb | 60k | 7.2 GB | 21.6 GB | 36 GB |
| 30 kb | 150k | 18 GB | 54 GB | 90 GB |

- VNTR graphs with 100+ haplotypes and copy-number variation easily reach N = 2-5 × L, because insertions of whole copies create new parallel nodes. **SPOA is not viable on whole 10-30 kb VNTR regions with hundreds of haplotypes.**

### Traceback

- **Scalar:** `sisd...cpp:372-462` (linear), `:553-679` (affine) and `:781-926` (convex).
- **Method:** traceback re-derives each move by checking equality against the stored H, F, E, O and Q values.
- **Move priority:**
  1. Diagonal through `inedges[0]` first, then through the other in-edges in order.
  2. Vertical (graph gap).
  3. Horizontal (sequence insertion).
- **Gap runs:** after a gap move, `extend_left` or `extend_up` loops follow the E or F chain to the end of the gap (`sisd...cpp:645-674`).
- **SIMD:** `simd...hpp:900-1066` and `:1264-1510`. One vector (`H`, `H_pred` for each predecessor, profile) is unpacked into small scratch arrays at a time, and it is reloaded only when the row changes or a lane boundary is crossed (`:1032-1034`). Border cases at `j_mod == 0` use the previous vector or `first_column` (`:970-990`).
- **NW ends:** the traceback finishes along row 0 or column 0 by emitting pure gaps (`:1042-1064`).
- **[inferred] Tie-breaking is biased** toward the diagonal and the oldest in-edge. Edges are appended in insertion order, so `inedges[0]` is usually the earliest sequence's path. Among co-optimal placements in a repeat, the traceback prefers matching through older paths and puts gaps as far left as possible. That is deterministic but arbitrary, with no repeat-aware normalization.

---

## 3. Banding, pruning, seeding

**There is none in SPOA itself.** No band, no X-drop, no seeds or anchors, and no early termination. Every sequence is aligned against the **entire** current graph in O(N·L).

The only restriction mechanism is **`Graph::Subgraph(begin, end, &mapping)` plus `UpdateAlignment`** (`graph.cpp:552-639`):

- `ExtractSubgraph` (`graph.cpp:552-573`) walks **backwards** from node `end` through in-edges and aligned nodes. It keeps nodes whose `id >= begin`.
  - The parameter names are swapped in the call at `graph.cpp:584`.
  - The filter compares **node IDs, not ranks**. It only makes sense when nodes 0..L−1 are a backbone added first, as in racon. **[inferred]**
- The subgraph is rebuilt with fresh IDs, its own edges (weights preserved) and aligned sets, and it is re-sorted (`graph.cpp:599-626`). Sequence labels are dropped (TODO at `graph.cpp:591`).
- `UpdateAlignment` maps subgraph node IDs back to the original graph (`graph.cpp:631-639`).
- Racon uses this for reads that do not span the whole window (`racon/src/window.cpp:93-108`).

**Consequences on long sequences:**

- **Memory:** quadratic, see the table above.
- **Time:** O(Σᵢ Nᵢ·Lᵢ·d̄/k), where d̄ is the mean in-degree and Nᵢ grows as sequences are added.
  - **[inferred]** For 456 × 30 kb with N around 100k, that is about 1.4×10¹² cell updates per region, or minutes to tens of minutes per region on a single thread even if memory were available.
- **Precision:** long inputs force int32, which halves throughput.
- **No banding:** a band would not save SPOA anyway. A sequence that differs by copy number needs a band as wide as the copy-number difference times the motif length.

---

## 4. Progressive strategy and order dependence

- **Order:** the library does no ordering. The caller chooses the order, and each sequence is aligned to the graph built from all earlier ones, then merged (README example; `main.cpp:282-321`).
- **Backbone:** the first sequence defines the backbone, either through an empty alignment or as the first chain.
- **No guide tree, no realignment, no iterative refinement.** The result depends on input order.
- **Strand:** `-s` aligns both strands and keeps the better one (`main.cpp:292-309`), which doubles the cost.
- **What downstream tools do about order:**
  - Racon puts the draft window first as the backbone, then sorts layers by start position (`racon/src/window.cpp:79-86`).
  - smoothxg sorts the path ranges in a block **longest first** (`smoothxg/src/blocks.cpp:206-219`; hard-coded `order_paths_from_longest = true`, `smoothxg/src/main.cpp:363`).
  - smoothxg also **deduplicates identical sequences**, adding each unique sequence once with weight = multiplicity (`smoothxg/src/smooth.cpp:720-745`, `:766-768`).

---

## 5. Consensus and MSA output

### Heaviest bundle (`graph.cpp:467-550`)

- **Forward pass in rank order** (`graph.cpp:476-490`):
  - Each node takes its **single heaviest in-edge**. Ties go to the predecessor with the higher accumulated score.
  - The node's score is that edge weight plus the predecessor's score.
  - The global max node is tracked.
- **Branch completion** (Lee 2003): if the max node is not a sink, `BranchCompletion` (`graph.cpp:512-550`) does the following, repeating until it reaches a sink (`graph.cpp:497-499`):
  - It invalidates competing tails of the max node's successors.
  - It rescores every node after it.
  - It moves to the new maximum.
- **Traceback** through `predecessors` gives `consensus_` (`graph.cpp:503-509`).
- **[inferred]** This is greedy and local. Each node picks its heaviest incoming edge, so in repeats the consensus can stitch together branches and produce a copy number that **no input haplotype has**.

### Other consensus variants

- `GenerateConsensus(min_coverage)` (`graph.cpp:378-387`) **drops** consensus bases whose `Coverage()` (distinct labels touching the node) is below the threshold. This trims; it does not re-path.
- `GenerateConsensus(summary, verbose)` (`graph.cpp:411-465`) returns per-consensus-position coverage. In verbose mode it returns a `(num_codes+1) × consensus_len` base and gap count matrix computed through MSA columns. Racon uses the non-verbose form to trim low-coverage ends (`racon/src/window.cpp:122-146`).

### MSA rows (`graph.cpp:322-367`)

- Column index = position in `rank_to_node_` after collapsing each aligned group into one column (`graph.cpp:326-333`).
- Each row starts as all `-`, and the sequence's path is written through `Successor(i)` (`graph.cpp:346-356`).
- `include_consensus` adds the heaviest-bundle row.
- Row width = number of column groups. Every insertion node that was never merged creates its own column.

### Graph outputs

- **DOT:** `PrintDot` in the library (`graph.cpp:641-681`). It shows consensus nodes and edges highlighted and aligned sets as dotted edges.
- **GFA:** only in the CLI, `PrintGfa` (`main.cpp:123-203`).
  - One `S` line per **single-base node**, with no unchopping.
  - `L` lines carry an `ew:f:` edge weight.
  - `ic:Z:true` marks consensus nodes.
  - `P` lines come from label walks, reversed for reverse-strand inputs, with an optional `Consensus` path.
- **Library graph access:** `nodes()`, `edges()`, `rank_to_node()`, `sequences()`, `consensus()`, `Successor()` and `aligned_nodes` (`graph.hpp:102-132`). This is how smoothxg builds its odgi block graph (`smoothxg/src/smooth.cpp:2583-2610`).

### Utilities

- `Subgraph` and `UpdateAlignment` (§3).
- `Clear()` (`graph.cpp:683-692`).
- Cereal `save`/`load` (`graph.hpp:194-282`).
- No pruning, no edge-weight thresholding, no node removal, and no graph-to-graph merging. Sequences can only be added.

---

## 6. Behaviour on tandem repeats

What in the design **causes** misaligned repeat copies:

1. **Order dependence plus greedy merging.** Each haplotype is aligned to a graph frozen from earlier haplotypes, and the merge is never revisited.
   - **[inferred]** Suppose haplotype A has n copies and B has n+1. B's extra copy becomes an insertion branch at whichever copy boundary the traceback picks.
   - A later haplotype C with n+1 copies may align its extra copy to B's branch if that scores best. If C's extra copy differs slightly, or tie-breaking picks another boundary, C opens a second parallel branch at a different copy.
   - Over hundreds of haplotypes this splits equivalent copies across columns ("column inflation"), or puts copy k of one haplotype in register with copy k±1 of another.
2. **The DAG cannot express "same motif, different copy count"** except as parallel branches. Because there are no cycles, every distinct copy-number layout adds new nodes that never merge with each other (insertions always create new nodes, `graph.cpp:210-211`).
3. **Arbitrary tie-breaking** (diagonal first, then `inedges[0]`; §2). In an exact tandem array, any copy is a co-optimal place for an indel. SPOA picks by edge insertion order and pushes the indel left. It does not normalize to a canonical repeat boundary.
4. **No anchoring.** Nothing pins unique flanks or unique sub-motifs before the DP runs.
   - In NW mode the alignment is globally optimal for one sequence against the graph. Within a long array, though, the optimal path for a sequence with a copy-number change can shift copy boundaries with no penalty.
   - In SW (smoothxg's default local mode) and OV modes, unaligned flanks become dangling branches that are never merged (`graph.cpp:198-200`).
5. **Scoring cannot separate a whole-copy indel from scattered gaps** except through the convex model.
   - The convex model, with a cheap second extension piece (e.g., `c = −1`), favours one long gap of a whole copy over several small ones. This is the right tool for VNTRs.
   - Linear gaps (racon's setting) encourage scattering.
6. **Greedy heaviest-bundle consensus** can mosaic the majority branch at each copy position into a copy number no input has. This matters if plegma ever uses a SPOA-style consensus as a backbone or guide.
7. **Memory and time blow-up on 10-30 kb** (§2-3). Downstream tools avoid it by chopping:
   - Racon uses 500 bp windows (`racon/src/main.cpp:46`).
   - smoothxg uses target POA lengths of 700/1100 in pggb (`pggb/pggb:39`), with sequences cut at 2 × target (`smoothxg/src/main.cpp:376`).
   - Chopping at fixed positions inside a VNTR creates block boundaries in the middle of the array that are not homologous across haplotypes. **[inferred]** This is a likely source of the very MC/pggb VNTR misalignments plegma targets.

What **mitigates** these problems:

- Deduplicating and weighting identical alleles (as smoothxg does) shrinks the problem.
- Ordering longest-first means the first backbone carries the most copies. Later, shorter alleles then align as deletions relative to an existing chain rather than creating insertions. **[inferred]** That cuts down parallel insertion branches.
- Convex gaps.
- NW mode with flanks included.

---

## 7. Downstream use

| Tool | How SPOA is used (verified in source) | API features relied on |
|---|---|---|
| **racon** | One engine per thread: `AlignmentEngine::Create(kNW, match=3, mismatch=-5, gap=-4)`, linear gaps, `Prealloc(window_length=500, 5)` (`racon/src/polisher.cpp:179-183`, `racon/src/main.cpp:46-53`). Per window, the draft segment goes in first as backbone with qualities. Layers are sorted by begin position. Reads that do not span the window use `Subgraph` + `UpdateAlignment`. The window is skipped if it has fewer than 3 sequences (`racon/src/window.cpp:65-149`). | Quality-weighted `AddAlignment`, `Subgraph`, `UpdateAlignment`, `GenerateConsensus(&coverages)` for coverage-based end trimming. |
| **smoothxg** (and pggb, where SPOA is the default) | Default params `1,4,6,2,26,1` (the minimap2 asm20 values) are negated into `m=1, n=-4, g=-6, e=-2, q=-26, c=-1`, which is **convex** (`smoothxg/src/main.cpp:322-327`, `smoothxg/src/smooth.cpp:2098-2103`). The mode is **local (SW) by default**; `-Z` selects global (`smoothxg/src/main.cpp:133-135`, `:487`; `smoothxg/src/smooth.cpp:752-755`). Sequences are deduplicated with weights and sorted longest first, with a small padding of flanks that is trimmed afterwards. pggb defaults are SPOA, local mode, and `-G 700,1100` (`pggb/pggb:39`, `:217-218`). | `AddAlignment(alignment, seq, weight)`, `GenerateConsensus()`, `GenerateMultipleSequenceAlignment(add_consensus)` for MAF output, and direct graph walking (`nodes()`, `outedges`, `sequences()`, `Successor`) to build the odgi block graph. |
| **Cactus / Minigraph-Cactus** | **No SPOA reference** in the Cactus clone (grep finds nothing). `cactus_bar` uses abPOA through `bar/inc/poaBarAligner.h:15`. The local clone lists `baseAligner="pecan|abpoa|minipoa"` (`poaBarAligner.h:37`). I did not verify whether that clone is pristine upstream. | n/a |
| Others (from memory, not verified here) | Medaka through pyspoa; Raven through racon. abPOA reports a Racon-abPOA variant. | |

---

## 8. Performance, complexity and benchmarks

- **Time per sequence:** O(Σ_v (indeg(v)+1) · L/k + N · (L/k) · log k). The log k term is the prefix max. Each topological sort is O(V+E). MSA generation costs O(Σ path · out-degree · labels) because of `Successor` (§1).
- **Space:** {1,3,5} × N × L × {2,4} bytes. There is no linear-space or banded alternative.
- **Racon paper (Vaser et al. 2017, Genome Res 27:737)** describes SPOA as follows:
  - "a SIMD version of the robust POA framework", "inspired by the Rognes and Seeberg Smith-Waterman intra-set parallelization approach".
  - "places the SIMD vectors parallel to the query sequence".
  - 8 lanes at 16-bit (4 at 32-bit), with 8-bit "insufficient".
  - Gap matrices "stored entirely in memory".
  - Horizontal processing "not performed using SIMD operations due to data dependencies", handled by "shifting and masking". This matches the prefix-max code.
  - Claimed complexity reduction: from O((2n_p+1)n|V|) to about O((2n_p/k+1)n|V|).
  - I did not find a SPOA vs poaV2 speedup figure in the main text.
- **abPOA paper (Gao et al. 2021, Bioinformatics 37:2209):**
  - Abstract: "up to 10 times faster" (bioRxiv preprint: "up to 15 times").
  - NanoSim sets: 2.7-9.5× faster than SPOA with adaptive banding.
  - The speedup grows with sequence length (tested 300 / 1000 / 5000 bp at 3-50× depth).
  - Consensus error rate is "comparable".
  - Peak memory is "similar" to SPOA with and without banding. abPOA also keeps full matrices.
  - abPOA "uses the same SIMD parallelization strategy as SPOA".
  - The text gives no absolute numbers; they are in the supplementary tables.
- **POASTA paper (van Dijk et al., Bioinformatics 2025, btae757):**
  - Notes that SPOA "still computes the full DP matrix and thus does not ameliorate demands on memory usage".
  - Reports POASTA at 4.1× average and up to 9.8× speedup with lower memory on bacterial gene sets, and megabase-scale graphs that SPOA and abPOA cannot build.

---

## 9. Advantages and disadvantages (judged against plegma's goal)

### Advantages

1. Simple, readable reference implementation of **exact** POA. NW alignment of one sequence against the graph is optimal under the chosen scoring.
2. **Three gap models, including convex (two-piece affine).** Convex is what VNTR copy indels need, with a cheap long-gap extension.
3. Clean SIMD design along the query axis:
   - Query profile.
   - One-lane shift with carry for the diagonal.
   - log₂k prefix max for horizontal gaps.
   - First column held outside the vectors.
   - Automatic int16/int32 selection from a provable score bound.
4. **Aligned-node clique + contiguous topological emission.** This gives MSA columns for free and caps each column at one node per base.
5. Per-base weights (qualities) and per-sequence weights, which make deduplication with multiplicity trivial.
6. A `Subgraph`/`UpdateAlignment` mechanism to align against a restricted region and map back.
7. Small dependency surface; MIT license; SIMDe portability (NEON on Apple Silicon is possible).
8. Battle-tested in racon, smoothxg and pggb at scale on short windows.

### Disadvantages

1. **Full O(N·L) score matrices**, 1/3/5 of them, int32 for long inputs. 10-30 kb × hundreds of haplotypes needs tens of GB (§2 table). Buffers never shrink.
2. **No anchoring, banding or seeding.** Every sequence is aligned to the whole graph, and run time grows with graph size.
3. **Strict progressive order dependence.** No guide tree, no refinement, no realignment of earlier sequences.
4. **Arbitrary tie-breaking** (oldest in-edge, diagonal first). Repeat-copy indels are placed without normalization.
5. **The DAG with no node merging after insertion** fragments copy-number variants into parallel branches and inflates MSA columns.
6. **Asymmetric "semi-global" mode** (free graph prefix and sequence suffix only), plus a possible OV-mode negative-score edge case in SIMD (§2, inferred).
7. **Local mode** leaves unaligned flanks as unmerged branches. This is smoothxg's default.
8. **Pointer-chasing graph:**
   - Per-node vectors.
   - Labels as vectors per edge.
   - `Successor` is a linear scan.
   - `Coverage()` allocates an `unordered_set`.
   - A full topological sort after every insertion.
9. Greedy heaviest-bundle consensus can produce copy numbers that no haplotype has.
10. Single-base GFA nodes with no unchopping, and GFA output only in the CLI. No graph pruning or simplification utilities.
11. The `Subgraph` filter is **node-ID based** and only correct for a backbone-first construction.

---

## 10. Ideas to steal and pitfalls to avoid

### Ideas to steal

1. **Column model.** Keep a per-column "aligned set" of at most 4-5 nodes, one per base. Emit each group contiguously in topological order so the MSA column index is just the group rank. In plegma, store `column_id` on each node directly instead of a clique of pointers, and keep the topological order over columns.
2. **Query-axis SIMD kernel, if plegma vectorizes this way:**
   - Precomputed query profile per base.
   - Diagonal = predecessor row shifted one lane, carrying the last lane across vectors.
   - Vertical = lane-wise add.
   - Horizontal = log₂k prefix max with doubled-penalty vectors and masks, plus a scalar first column.
   - Combine multiple predecessors by elementwise max over full row passes.
   - This maps directly onto Mojo `SIMD[DType.int16, W]` with `shift_left`/`shift_right` lane shifts. **[inferred]**
3. **Choose the precision from a worst-case bound per window.** After sporos anchoring, plegma's windows are small, so **int16** will almost always fit, which gives twice the lanes of SPOA's 30 kb int32 case. With difference encoding, int8 may even be possible.
4. **Convex two-piece affine gaps.** For VNTRs, set the second piece to roughly the motif-scale cost (high open, very low extend) so whole-copy indels are one gap.
5. **Restricted-graph alignment with an ID remap** (`Subgraph` + `UpdateAlignment`). This is exactly the operation plegma needs: align the haplotype segment between two MUM or SMEM anchors to the graph segment between the anchor columns, then map back. Define the region by **column/rank interval between anchor nodes**, not by node ID.
6. **Deduplication with weights** (smoothxg pattern). HPRC VNTR alleles are often identical across haplotypes. Align each unique allele once and carry its multiplicity as the weight.
7. **Longest-first ordering** (smoothxg) as a cheap baseline, so later alleles appear as deletions against existing chains rather than as new insertions. Better: order by similarity or a guide tree, or by copy number.
8. **Heaviest bundle with branch completion** as a cheap way to get a backbone or reference path if one is needed. Weight it by unique-allele multiplicity.

### Pitfalls to avoid

1. **Never allocate full N × L int32 matrices for whole regions.**
   - Run DP only inside anchor-bounded windows.
   - Store 2-4 traceback bits per cell, or use banding or Hirschberg-style linear-space recursion, instead of SPOA's "re-derive from stored scores".
   - Free or shrink buffers between regions.
2. **Do not let input order decide copy registration.** Anchor first with MUMs/SMEMs so homologous copies are pinned across haplotypes. Consider realigning early sequences against the final graph, as in iterative refinement.
3. **Deterministic, repeat-aware tie-breaking.** Do not inherit "oldest in-edge wins". When gaps are co-optimal inside a periodic window, normalize to a canonical position (e.g., left-normalize relative to motif phase), and break predecessor ties by weight or by anchor consistency.
4. **Merge parallel insertion branches.** After each window, collapse insertion nodes that carry the same sequence at the same column interval into a shared chain, or realign new insertions against existing insertion branches. SPOA never does this, which is why columns inflate.
5. **Avoid asymmetric or local modes for windows between anchors.** Use NW with both ends fixed at the anchors. Use free-end variants only at region edges, and make them symmetric.
6. **Avoid SPOA's graph data-structure costs:**
   - Use CSR or flat arrays for edges.
   - Store each sequence's path as an explicit node or column array instead of label vectors per edge searched linearly.
   - Use bitsets or counts for coverage.
   - Update the topological order incrementally within a window.
7. **Do not cut long sequences at fixed lengths** (the racon and smoothxg approach). Cut only at anchors, which are homologous by construction. Fixed-length cuts inside a VNTR are not homologous across haplotypes.
8. **Do not trust greedy consensus as a truth or backbone in repeats.** If one is needed, check that it is a real input path (or the closest one), not a mosaic.
9. **Watch the OV-style negative-score pitfall.** Any max-reduction should start at −∞, not 0, unless the mode is local (`simd...hpp:242`, `:1187`).
