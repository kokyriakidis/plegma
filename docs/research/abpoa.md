# abPOA: source review for plegma

Target: https://github.com/yangao07/abPOA, shallow clone at commit `2e095ba` (2026-09-17, "fix pypi installation"), version string 1.5.7 (`src/abpoa.c:19`). I read the source and did not build or run it. Paths below are relative to the repo root. Anything marked **(inferred)** is my own reasoning from the code, not something the code or docs state.

Papers and claims:
- Gao Y. et al., "abPOA: an SIMD-based C library for fast partial order alignment using adaptive band", *Bioinformatics* 37(15):2209–2211 (2021), doi:10.1093/bioinformatics/btaa963.
  - Adaptive band: `B_start = max{0, min{M_start, L−R} − w}`, `B_end = min{L, max{M_end, L−R} + w}`, `w = b + f·L` (b=10, f=0.01).
  - Speed: 2.7–9.5× faster than SPOA on NanoSim sets ("up to 10×"), with consensus accuracy similar to SPOA.
  - Memory: peak memory about the same as SPOA, with or without banding.
  - Limitations the paper states itself: the band is "expected" to cover the optimal path but this is not guaranteed, and results depend on input order.
- No separate paper exists for the v1.2+/v1.5 minimizer seeding and "progressive guide tree". A web search found none. Both are described only in the README ("Minimizer-based seeding mode" and "Minimizer-based progressive tree").

---

## 1. Graph representation

**Node struct** (`include/abpoa.h:91-105`, `abpoa_node_t`). Each node is one base, stored as an array of structs (`abpoa_graph_t.node`, `abpoa.h:107-112`).
- Edges are adjacency arrays on each node:
  - `in_id[]` and `in_edge_weight[]`
  - `out_id[]` and `out_edge_weight[]`
- `read_ids` is a `uint64_t**`: one bitset per **out-edge**, each `read_ids_n = ceil(n_seq/64)` words long (`abpoa.h:99`). For 456 haplotypes that is 8 words per edge.
- `read_weight[]`: per-read weight, used only with quality scores.
- `n_read` and `n_span_read`.
- `aligned_node_id[]`: the mismatch-sibling group, i.e. nodes that occupy the same MSA column with a different base (`abpoa.h:101`).
- Node 0 is a virtual source and node 1 a virtual sink (`abpoa.h:27-28`, `src/abpoa_graph.c:103-113`).
- Node ids are assigned in insertion order and never renumbered. Topological order lives in separate maps: `index_to_node_id` and `node_id_to_index`.

**Graph-level side arrays** (`abpoa.h:110`), indexed by node id:
- `node_id_to_max_pos_left/right`: the band's propagated max-score columns.
- `node_id_to_max_remain`: remaining length to the sink.
- `node_id_to_msa_rank`: MSA column.

**Adding an edge** (`src/abpoa_graph.c:480-556`):
- With `check_edge`, an existing edge is found by linear scan of in/out lists. Its weight increases by `w` (1, or a per-base quality).
- The read's bit is set on the out-edge bitset (`abpoa_set_read_id`, `src/abpoa_graph.c:465-469`).
- Bitsets are grown lazily when the total read count rises (`src/abpoa_graph.c:530-542`).

**Fusing an alignment into the graph** (`abpoa_add_subgraph_alignment`, `src/abpoa_graph.c:680-774`):
- Match: reuse the node.
- Mismatch: reuse an aligned sibling that has the same base (`abpoa_get_aligned_id`, `src/abpoa_graph.c:439-448`), else create a new node and add it to the sibling clique (`abpoa_add_graph_aligned_node`, `src/abpoa_graph.c:455-463`). The clique is maintained all-to-all, so cost is O(k²) per column.
- Insertion or clip: always create **fresh** nodes (`src/abpoa_graph.c:747-758`). An inserted base is never aligned to any other node, so a later read can merge into an insertion branch only by aligning to it in its own DP.
- Deletion: nothing is added.
- After every read the graph is fully re-sorted (`abpoa_topological_sort`, `src/abpoa_graph.c:767-768`) and `n_span_read` is updated over the whole topological range (`src/abpoa_graph.c:559-571`).

**Topological order** (`abpoa_BFS_set_node_index`, `src/abpoa_graph.c:221-266`):
- Kahn's algorithm, run as BFS from the source.
- A node is enqueued only when it **and all its aligned siblings** have in-degree 0. The siblings are then enqueued right after it, so aligned-node groups end up contiguous in the order.
- After sorting, in and out lists of every node are bubble-sorted by edge weight, descending (`src/abpoa_graph.c:192-219`). This is O(deg²) per node.

**"Remain"** (`abpoa_BFS_set_node_remain`, `src/abpoa_graph.c:268-309`):
- A reverse BFS from the sink computes `remain[v] = remain[heaviest out-neighbour] + 1`.
- So "max_remain" is the remaining length along the **heaviest-edge path**, not the longest or shortest path. This feeds the band and Z-drop.

**MSA columns** (`abpoa_DFS_set_msa_rank`, `src/abpoa_graph.c:359-410`):
- A stack-based (DFS) topological traversal with the same sibling gating as above. Siblings get the same rank.
- The column of a node is then `max(rank of node, rank of siblings)` (`src/abpoa_output.c:136-145`).

**Subgraph helpers** (`src/abpoa_graph.c:595-678`, `abpoa_subgraph_nodes`):
- Widen a node range until it is closed under in-edges and out-edges.
- Used by the public sub-graph API (`sub_example.c:118`), but **not** by the internal seeding path.

**Restoring a graph** (`abpoa_restore_graph`, `src/abpoa_seq.c:608-673`):
- From MSA FASTA, each column index becomes a node and siblings are created per column (`abpoa_fa_parse_seq`, `src/abpoa_seq.c:572-606`), so **column structure is preserved**.
- From GFA, only S and P lines are read and L lines are ignored (`src/abpoa_seq.c:467-570`). Segments are expanded into base chains. **(inferred)** Aligned-sibling relations are not reconstructed from GFA, because the P-line parser never calls `abpoa_add_graph_aligned_node`.

## 2. Dynamic programming

**Scoring** (`src/abpoa_align.c:88-99`, `src/abpoa_align.h:9-14`):
- Defaults: match 2, mismatch 4, O1=4, E1=2, O2=24, E2=1.
- The gap mode is implied by the parameters:
  - `O1 == 0` → linear.
  - `O2 == 0` → affine.
  - otherwise → convex, i.e. dual-affine: `min{O1+g·E1, O2+g·E2}`. This is the default (`src/abpoa_align.c:104`).
- An optional substitution matrix file is supported (HOXD70 is bundled). Cactus overrides the matrix with HOXD70, per issue #14.
- Optional `-G inc_path_score` adds `round(log(edge_w/node_w))`, floored at −20, to M and E transitions (`src/abpoa_graph.c:429-437`). The usage text marks it as "under development" (`src/abpoa.c:96-97`).

**Modes** (`abpoa.h:6-8`):
- Global (default), local, and extension.
- Local mode **turns the band off** (`wb=-1`, `src/abpoa_align.c:167`).
- Seeding and progressive order are used only in global mode (`src/abpoa_align.c:443`, `:507`).

**DP over the DAG** (rows = graph nodes in topological order, columns = query positions 0..qlen):
- The first row (the window's begin node) is initialised with global leading gaps (`src/abpoa_align_simd.c:635-688`).
- For each node, H/E are the max over **all in-window predecessors**. M comes from the predecessor's H shifted by one lane, E from the predecessor's E. Then F, the gap along the query within the row, is computed (`src/abpoa_align_simd.c:727-1074`; convex version at `:935-1074`).
- Naming: E = deletion, a vertical move along the graph. F = insertion, a horizontal move along the query.
- E is stored "pre-decremented" for successors (`:1070-1071`).
- Global end score: max over the in-neighbours of the end node, at column `min(qlen, dp_end)` (`simd_abpoa_global_get_max`, `:1092-1105`). If the band never reached `qlen`, the traceback silently appends a trailing insertion (`:121`, `:201`, `:315`).
- The subgraph is defined by an `index_map`: nodes in topological range `[beg_index, end_index]` that are **reachable forward from the begin node** (`:1257-1269`). Predecessors outside the map are dropped (`:550-559`).

**SIMD:**
- Vectorised along the **query** within one node row. The memory layout is plain sequential, not Farrar-striped (`simd_abpoa_var`, `:463-495`).
- The horizontal F dependency is resolved inside each vector by a log₂(lanes)-step shift-and-max prefix scan (`SIMD_SET_F`, `:691-725`, using `GAP_E1S[k] = 2^k·E1`). The carry between consecutive vectors goes through the `first` register.
- The query profile is precomputed per residue (`:531-539`).
- Lane width is 16-bit when `max(qlen·max_mat, max(qlen,gn)·E1+O1)` fits in int16, otherwise 32-bit (`:1293-1302`). **(inferred)** With match=2, a window over ~16.3 kb, or a graph over ~16k nodes, forces 32-bit lanes and halves throughput.
- Backends: SSE2/SSE4.1/AVX2/AVX512BW, with runtime dispatch (`src/abpoa_dispatch_simd.c`) and SIMDe for ARM.

**Adaptive band** (`src/abpoa_align.h:31-35`, `GET_AD_DP_BEGIN/END`):
- `beg = max(0, min(max_pos_left[v], qlen − (remain[v]−remain[end]−1)) − w)`
- `end = min(qlen, max(max_pos_right[v], qlen − (remain[v]−remain[end]−1)) + w)`
- `w = wb + wf·qlen` (`src/abpoa_align_simd.c:474`).
- After each row, the leftmost and rightmost columns holding the row max (`simd_abpoa_max_in_row`, `:1107-1119`) are propagated +1 to every out-neighbour as min-left / max-right (`simd_abpoa_ada_max_i`, `:1121-1130`).
- The row start is also clamped to be no earlier than the earliest predecessor start (`:744-747`).
- Only vectors overlapping the band are computed.
- **The band saves time only, not memory.** The full `gn × ceil((qlen+1)/lanes)` matrix is allocated: ×1 for linear, ×3 (H,E,F) for affine, ×5 (H,E1,E2,F1,F2) for convex (`src/abpoa_simd.c:52-83`). This is the cause of issue #14, where an allocation of 549 GB was attempted.
- **(inferred) Behaviour on a large copy-number difference:**
  - The band is the union of the "start diagonal", tracked by the propagated row maxima, and the "end diagonal", `qlen − remain`, plus w. So **one** large indel relative to the heaviest path is covered, because the band stretches between the two diagonals.
  - It fails in three cases:
    - (a) Two or more opposite indels, each larger than w, in the same DP. Example: +3 copies in one VNTR and −3 in the next, so net length is ≈0 and the end diagonal does not help.
    - (b) The query follows a non-heaviest branch whose length differs from the heaviest path, because `remain` is a heaviest-path length (`src/abpoa_graph.c:288-296`).
    - (c) The row max is misleading. In perfect repeats, ties widen the band, which helps. In imperfect repeats a wrong copy can win the row max and steer the band.
  - w is small: 20 bp for a 1 kb query, 110 bp for 10 kb, and about 15 bp inside a 500 bp seeded window.
  - On failure there is **no detection and no retry**. The result is a suboptimal alignment. In some band configurations it is a fatal `"sn_i < min_pre_beg_sn"` (`:804`, `:910`, `:1044`) or a backtrack error (issue #9: "bug related to the banding", sequences of 264 vs 443 bp).
  - Issue #10 shows misaligned long ONT reads that were fixed by widening `-f`.

**Z-drop** (`set_extend_max_score`, `:1082-1090`):
- Extension mode only, and disabled by default (`zdrop=-1`, `src/abpoa_align.c:105`). The CLI help line is commented out (`src/abpoa.c:101-104`).
- It breaks the row loop when `best − score > zdrop + E1·|Δremain − Δj|`, using the `remain` difference as a proxy for graph distance.
- `end_bonus` is parsed but unused ("TODO", `:1132`).

**Traceback** (`:116-458`):
- There is no traceback matrix. The traceback re-derives each move by testing score equalities against the stored H/E/F matrices, which is why all of them are kept in memory.
- Default priority is M, then E, then F, scanning predecessors in edge-weight order, which yields left-most gaps.
- `-R` puts gaps on the right and `-J` keeps the last gap at the end.
- Ties are therefore broken by predecessor weight order and gap-side preference. Issue #78: forward vs reverse-complement input give different MSAs because gaps always go left.
- The output is a graph CIGAR of packed 64-bit ops (`abpoa.h:45-50`, `src/abpoa_align.h:54-73`).

## 3. Seeding and partitioning (`-S`, off by default: `disable_seeding=1`, `src/abpoa_align.c:146`)

**Minimizers:**
- minimap2's sketch is copied verbatim (`src/abpoa_seed.c:46-168`), with k=19 and w=10 (`src/abpoa_align.h:16-17`).
- The sketch is forward-only unless `-s` (both strands).
- There is **no high-frequency or repeat minimizer filter anywhere**.

**Anchors are pairwise between consecutive sequences in processing order, not to the graph:**
- `abpoa_build_guide_tree_partition` (`src/abpoa_seed.c:717-756`) loops `tid = order[i-1]`, `qid = order[i]`.
- `collect_anchors1` (`:344-377`) merge-joins the two sorted minimizer lists and emits the full **cross product** of occurrences for every shared hash (`:353-369`).
- **(inferred)** In a tandem array with c copies, each minimizer yields about c² anchors.
- The author confirms the pairwise design in issue #37: "we just need the anchors between every two adjacent sequences."

**Chaining** (`abpoa_dp_chaining`, `:500-591`; the LIS chainer at `:593-701` is dead code):
1. Stage 1 is minimap2-like DP with hard-coded parameters (`:505`):
   - `max_bw=100`, `max_dis=100`, `max_skip_anchors=25`, `max_non_best_anchors=50`, `min_local_chain_score=100`.
   - Anchor gap score: `min(dq,dt,k) − (log2(|dq−dt|)/2 + 0.01·k·|dq−dt|)` (`:482-495`).
   - Anchors more than 100 bp apart, or with a diagonal shift over 100, cannot be linked.
   - Local chains are taken from chain ends with score ≥ 100. Chains that branch off an already-used chain are dropped (`:565-583`).
2. Stage 2 is DP over the local chains. It requires only collinearity, uses **no gap penalty**, and scores `chain_score[j] + score[i]` (`:398-431`, especially `:419`). Overlapping chains are scored by trimming the chain (`get_local_chain_score`, `:379-392`).
3. **(inferred) Possible bug:** at `:566`, `strand = anchors->a[i] >> 63` indexes the anchor array with the local-chain loop index `i` instead of the anchor index `j`. Local-chain strands can therefore be mislabelled when `-s` is on.

**Windowing** (`:434-468`):
- Walking the best global chain backwards, an anchor is kept only if it is at least `min_w + k` = 519 bp after the previously kept one in **both** sequences (`-n`, default 500, `src/abpoa_align.h:18`).

**Stitching** (`abpoa_anchor_poa`, `src/abpoa_align.c:209-310`):
- Each anchor's target position is mapped to a node through `tpos_to_node_id`, which is the previous sequence's position→node map filled by the previous `add_alignment` (`:299`). The two maps are swapped per read (`:305`).
- Each inter-anchor window is aligned with `abpoa_align_sequence_to_subgraph(beg_id, end_id)`, with both end nodes excluded (`:274`).
- The k anchor bases are **forced as k M-ops** onto the previous read's nodes, with no DP (`:280-285`).
- All window CIGARs are concatenated and fused into the graph **once** per read (`:299`).
- Strand comes from the anchor's strand bit XOR the previous read's strand (`:226`). Anchors are re-flipped when the previous read was reverse-complemented (`:236-262`).
- Incremental (`-i`) plus seeding is unimplemented: there are TODO stubs at `:453-460` and `:517-526`.

**Ambiguous seeds (repeats):** nothing handles them explicitly. Every occurrence pair becomes an anchor, and the chaining DP plus the gap-free chain-of-chains decides which one is used. **(inferred)** For VNTRs, this means:
- (a) Anchors pairing copy i with copy j (i≠j) are about as good as the true ones.
- (b) Where the copy-number difference gets absorbed, i.e. where the chain jumps diagonal, is essentially arbitrary and is decided per pair, not against the graph.
- (c) The k-mer forced match then **locks** that choice in, and later DP cannot undo it.
- (d) A read is anchored only to the immediately preceding read's path, so a read with a different copy number is pinned to that read's copies even if a better-matching haplotype is already in the graph.
- (e) The DP of a window can only use nodes reachable from the window's begin node.

Related issues:
- #66: seeding produced 52 GB of memory, a crash, or 10× larger output on a Cactus input.
- #78: on an HPRC graph region, the author says seeding plus progressive order "does not work well for this repeat region sequences."

## 4. Progressive strategy

- **Default order:** input order. `-L` sorts by length, descending, using an O(n²) swap sort with full sequence copies (`src/abpoa_align.c:374-391`).
- **Cactus practice:** passes sequences in descending length order (issue #48 comment).
- **`-p` "guide tree"** (`abpoa_build_guide_tree`, `src/abpoa_seed.c:244-337`). It is **an ordering, not a tree**, and no profile–profile merging happens:
  - Similarity is a multiset Jaccard of minimizers. Hits are `Σ min(cnt_i, cnt_j)` over shared minimizer hashes (`:252-291`).
  - Start with the most similar pair, then greedily append the sequence with the highest **sum** of similarities to everything already placed (`:310-332`).
  - **(inferred) Cost** is O(G·n²) for G distinct minimizers, because every hash group loops over all n² pairs (`:270-274`), plus O(n³) for the greedy step. At 456 sequences that is about 10⁵ pair updates per distinct minimizer.
- **Incremental addition:** every read goes sequence-to-graph and is fused immediately. There is no realignment, refinement pass, or iteration.
- **Order dependence:** explicit, and acknowledged by the paper and the author. Issue #78: "input order is very important for determining the number of gaps", because gaps already in the graph are free for later reads. Issue #48: homopolymer and tandem-repeat gap placement is inconsistent across rows, and `-p` "may not be able to differentiate them" because the similarity is only approximate.
- **SIMD determinism:** #66 notes output differs slightly between SSE2/AVX2/debug builds even without seeding.

## 5. Consensus and MSA output

**Heaviest bundle** (`abpoa_heaviest_bundling`, `src/abpoa_output.c:478-548`):
- A reverse-topological pass picks, per node, the **out-edge with maximum weight**, breaking ties by downstream score. `score = max_w + score[next]`.
- This is Lee's greedy heaviest-bundle, **not** a maximum-weight path DP.
- The path is then walked from the source (`abpoa_set_hb_cons`, `:376-392`).
- Per-base coverage and a logistic "phred" are attached (`:297-303`).

**Most-frequent alternative** (`-a 1`, `abpoa_most_frequent`, `:550-587`): per-column majority vote. A base is emitted only if its count ≥ the gap count (`:413`).

**Multiple consensus** (`-d 2..10`, `abpoa_multip_read_clu_kmedoids`, `:1136-1181`):
1. Build the RC-MSA.
2. Find candidate heterozygous columns where ≥2 alleles, gap counted as an allele, each have depth in `[min_het, n−min_het]`. Here `min_het = max(2, ceil(n·q)/2)`, with `q` defaulting to 0.25 (`:677-822`).
3. Merge columns with an identical read partition and keep a count.
4. Read distance is a weighted Hamming distance over the het columns: SNP weight 2, gap weight 1, times the count (`:824-863`).
5. k-medoids, seeded by the most frequent het partitions, at most 10 iterations (`:990-1134`).
6. k is decreased until every cluster has at least `min_het` reads and at least 80% of reads are clustered (`:1111`).
7. Heaviest bundle is then rerun with edge weights restricted to each cluster's read bitset (`get_inclu_edge_weight`, `:315-336`).

Notes on the clustering:
- **(inferred) Minor bug:** tie-breaking in `abpoa_update_kmedoids` only chooses between clusters 0 and 1 (`:1048-1051`), even for k>2.
- Issues #57 and #73 report missed clusters and alleles in short tandem-repeat data. #57 was fixed, #73 was a pre-1.5.5 limit to k≤2.

**RC-MSA** (`abpoa_generate_rc_msa`, `:150-193`):
- Rows are filled from each node's per-edge read bitsets. The column is the max rank across the sibling group.
- Consensus rows can be appended.
- Reverse-complemented rows are named `_reverse_complement` (`:80`).

**GFA** (`abpoa_generate_gfa`, `:196-295`):
- One single-base S line per node (no unchopping), L lines from in-edges, and a P line per read rebuilt from the bitsets.
- Reverse-complemented reads are emitted as reversed `-` walks.
- An optional consensus path can be added.
- Sibling and column relations are **not** written out.

## 6. Tandem repeats specifically

Design features that **cause** misaligned repeat copies:
1. **Minimizer seed ambiguity.**
   - No frequency filter, and the full cross product of occurrences is used (`src/abpoa_seed.c:353-369`).
   - Chaining has a 100 bp link and band limit and a gap-free chain-of-chains (`:419`, `:505`).
   - Anchors are forced as exact k-mer matches with no DP (`src/abpoa_align.c:280-285`), and only against the previous read (`:305`).
   - With ≥519 bp spacing, a multi-kb VNTR receives anchors *inside* the array whose copy pairing is arbitrary **(inferred)**.
2. **Band failure on large or multiple indels.**
   - The band is `w = 10 + 0.01·L` wide around the start and end diagonals. The end diagonal comes from the **heaviest path**.
   - Opposite copy-number changes in one DP, or a query that follows a minority-length branch, fall outside the band **(inferred)**. Failures are silent, or crash.
   - Workarounds the author suggests: increase `-b/-f` (#10), or disable the band with `-b -1` (#9).
3. **Greedy consensus collapse.**
   - Heaviest bundle follows the locally heaviest edge, so the consensus copy number is the one carried by the single heaviest edge at each branch, not necessarily any real haplotype **(inferred)**.
   - `-d` clustering relies on het columns, which misalignment inside the array pollutes: gaps count as alleles.
4. **Order dependence.**
   - Gaps already present in the graph are free for later reads (#78).
   - With leftmost gap placement, which copies count as "extra" depends on gap side and strand (#78, #48).
   - `-p` uses only a Jaccard score, which barely separates near-identical repeat haplotypes (#48).
5. **DAG-only graph.** Extra copies can only be insertion branches of fresh nodes (`src/abpoa_graph.c:747-758`). There is no repeat-aware or cyclic representation and no mechanism to merge equivalent insertion branches after the fact.
6. **Window confinement.** A window's DP sees only nodes reachable from the begin node within the topological range (`src/abpoa_align_simd.c:1257-1269`), so alignment to copies on another haplotype's branch can be impossible.

Design features that **avoid** some problems:
- Convex (dual-affine) gaps make long copy-number gaps cheap relative to affine: about 24+g vs 4+2g.
- The band covers the start-to-end diagonal span, so one large copy-number change per DP is handled.
- Row-max ties widen the band in perfect repeats **(inferred)**.

HPRC relevance:
- Issue #78 was filed by Glenn Hickey from Adam Novak's audit of **HPRC graph** alignments (CHM13 chr3:164033777-164033842). Cactus feeds abPOA with seeding and progressive order on.
- Issue #37 (Paten/Hickey) asks for **user-specified MUM anchors** "to facilitate the alignment of highly repetitive sequences (e.g. satellite arrays)". It remains open and unimplemented. This is exactly the plegma niche.
- Issue #48 (rlorigro) covers tandem-repeat and homopolymer gap inconsistency, suggests "seeding/chaining and locally computed guide trees", and is open.

## 7. Performance and complexity

- **Time per read:** Σ over nodes of (in-degree × band vectors). Worst case without a band is O(V·L·d/lanes). Plus an O(V+E) re-sort and O(Σdeg²) edge sort after every read (`src/abpoa_graph.c:192-219`, `:767-768`).
- **Memory:** O(V·L), regardless of the band:
  - 2 or 4 bytes × {1,3,5} matrices per cell (`src/abpoa_simd.c:52-62`).
  - **(inferred)** Example: convex at int32 for a 30 kb query against a 60k-node graph is about 36 GB.
  - Seeding cuts this to O(window_nodes × window_len) **only if anchors are found**. When they are not, memory explodes (#14, #66).
  - The single DP buffer only grows, never shrinks (`src/abpoa_simd.c:69-73`).
  - The `pre_index` arrays are re-malloc'd per row per alignment (`src/abpoa_align_simd.c:548-559`).
- **Guide tree:** O(G·n²) + O(n³) **(inferred)**.
- **Seeding:** O(Σ c²) anchors in repeats.
- **Benchmarks:** paper only. 2.7–9.5× faster than SPOA, similar accuracy and memory (100 clusters × {300, 1000, 5000} bp × {3, 10, 30, 50}× depth, NanoSim/PBSIM; plus Racon integration). Issue #66 has user numbers: 164–296 s and 6.9 GB without seeding, up to 52 GB with seeding on one Cactus block.

## 8. Advantages and disadvantages (judged against plegma's goal)

**Advantages**
1. A mature, fast SIMD POA kernel: query-vectorised rows, log-step prefix scan for F, and an adaptive 16/32-bit lane choice.
2. Dual-affine (convex) gaps by default, which suits copy-number indels.
3. The adaptive band covers both the start and end diagonal, so one large length difference per DP is OK at low cost.
4. A clean **subgraph alignment API** between two excluded anchor nodes, plus CIGAR concatenation and a single fuse per read. This is structurally what plegma needs between MUM anchors.
5. Per-edge read bitsets enable MSA rows, GFA paths, and cluster-restricted consensus with no per-read path storage.
6. Mismatch-sibling groups give a real column model, so RC-MSA output is cheap.
7. It can restore a graph from an MSA with columns preserved, which allows incremental alignment to an existing alignment.
8. It is battle-tested inside Cactus and Minigraph-Cactus, which gives a known baseline to beat on the same HPRC regions.

**Disadvantages**
1. Seeds are minimizers with no repeat filtering and a cross-product of occurrences, pairwise only to the previous read, forced as exact matches. This is the main cause of VNTR misalignment (#78 and #66, plus the code).
2. Chaining parameters are hard-coded (100 bp link and band), and the gap-free chain-of-chains places copy-number jumps arbitrarily.
3. The band fails silently, or crashes, on opposite or multiple indels and minority-branch lengths. There is no detection or fallback.
4. Memory is O(V·L) regardless of the band, which is impractical for 30 kb × large graphs unless seeding succeeds.
5. Strong order dependence. `-p` is a greedy Jaccard ordering, O(G·n²), weak on near-identical repeat haplotypes.
6. No refinement or realignment. Insertion branches are never merged and early errors are permanent.
7. Heaviest-bundle consensus is greedy and can produce a non-haplotype copy number. Multi-consensus relies on het columns polluted by misalignment.
8. The full topological re-sort and quadratic edge sort per read add overhead at 456 reads × tens of thousands of nodes.
9. `remain` is a heaviest-path length, so the band and Z-drop are biased toward the majority haplotype.

## 9. Ideas to steal and pitfalls to avoid

**Ideas to steal (concrete)**
1. **Two-diagonal band, made exact by anchors.** Inside a sporos-bounded window, precompute for each node the min and max path length from the window begin and to the window end. Allow columns `[max(minFromBeg, qlen − maxToEnd) − w, min(maxFromBeg, qlen − minToEnd) + w]`, intersected with abPOA-style propagated row maxima. This is sound for any copy number present in the graph, unlike the heaviest-path `remain` (`src/abpoa_align.h:34-35`).
2. **Band-hit detection with retry.** Not in abPOA. If the traceback touches a band edge, or the end cell is outside the band (abPOA's `best_j < qlen` path, `src/abpoa_align_simd.c:121`), redo the window with w doubled or unbanded. Small windows make this cheap.
3. **Kernel layout.** Rows = nodes, vectors along the query, an in-register prefix-max for F with `GAP_E·2^k` shifts and a carry register, and a 16/32-bit choice from a max-score bound (`:691-725`, `:1293-1302`). This maps directly onto Mojo SIMD.
4. **Subgraph DP between excluded anchor nodes**, concatenated CIGARs, and one fuse per read (`src/abpoa_align.c:268-299`). Keep this shape, but let sporos supply the anchors.
5. **Data structures:** read bitsets per out-edge for paths, MSA, and cluster consensus (`src/abpoa_graph.c:465-469`, `src/abpoa_output.c:106-123`), and sibling groups for columns.
6. **MSA-restore path** (`src/abpoa_seq.c:572-606`). Seed plegma's graph from the existing MC alignment of the flanks, or from a trusted reference MSA, and realign only the VNTR core.
7. **Dual-affine default**, with O2/E2 tuned so a full repeat unit costs less than scattered mismatches.

**Pitfalls to avoid**
1. **Never force anchors that are not unique.** Use MUMs, or SMEMs filtered by occurrence count, and verify uniqueness across the set (or the current partition) before treating them as hard matches. abPOA's cross-product of minimizer occurrences (`src/abpoa_seed.c:353-369`) is the anti-pattern.
2. **Anchor to the graph or profile, not just to the previous sequence** (`src/abpoa_seed.c:741-751`). Use multi-sequence MUMs or SMEMs that hit many haplotypes so window boundaries are consistent across all rows. abPOA's partitions differ per read pair, which breaks cross-haplotype column consistency **(inferred)**.
3. **Do not chain without gap cost.** abPOA's chain-of-chains (`:419`) lets a copy-number jump land anywhere. Plegma's recursive partitioning should place boundaries only on unique anchors flanking or inside arrays, and leave the copy-number decision to the DP inside the window.
4. **Do not allocate the full V×L matrix for a banded DP.** Store only band cells: per-row offset + length, as in abPOA's `dp_beg/dp_end`, plus compact storage. Or store traceback bits instead of five score matrices.
5. **Do not trust a silent band.** Detect and retry (idea 2). Avoid abPOA's fatal `"sn_i < min_pre_beg_sn"` paths.
6. **Reduce order dependence.**
   - Use a real guide tree, e.g. from sporos MUM/SMEM similarity or copy-number-aware distances.
   - Add at least one refinement round: remove each sequence and realign it.
   - Avoid abPOA's O(G·n²) pair counting at n=456.
7. **Normalise gap placement inside repeats** consistently on both strands, e.g. left-align to repeat-unit boundaries. Issue #78 shows strand-dependent MSAs from leftmost gap placement alone.
8. **Insertion branches need merging or realignment.** Fresh nodes per insertion (`src/abpoa_graph.c:747-758`) fragment extra copies across haplotypes into separate columns.
9. **Consensus is secondary for plegma, but do not use greedy heaviest-edge** if a representative is needed. Use a true max-weight path DP or a per-cluster path.
10. **Avoid a full topological re-sort per read.** Update the order incrementally per window, since inserted nodes lie between two known anchors.
11. Keep results deterministic across SIMD widths (#66).

Sources:
- Gao et al. 2021, *Bioinformatics* 37(15):2209 — https://academic.oup.com/bioinformatics/article/37/15/2209/6042680
- abPOA README (repo root).
- GitHub issues #9, #10, #14, #37, #48, #57, #66, #70, #73, #78 on yangao07/abPOA.
