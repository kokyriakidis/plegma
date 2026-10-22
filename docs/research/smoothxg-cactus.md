# smoothxg (PGGB), Minigraph-Cactus BAR/abPOA, and pgrealign, reviewed for plegma

Sources read (read only; nothing was built or run):

- smoothxg `9e8018d` (2026-08-21): `scratchpad/src/smoothxg`. Its submodule `sautocorr` was read from its own clone, `scratchpad/src/sautocorr`.
- pggb (HEAD): `scratchpad/src/pggb/pggb`, the shell driver. It shows the smoothxg flags PGGB actually passes.
- cactus `0657ef7` (2026-10-19): a sparse checkout of `bar/`, `src/cactus/`, `doc/` and `api/` in `scratchpad/src/cactus`.
- pgrealign `519ca6a` (2026-10-17): `scratchpad/src/pgrealign`. This is benedictpaten/pgrealign. cactus HEAD ships it as `cactus-pangenome --pgrealign`.
  - **This is the closest existing tool to plegma's goal.** It realigns the tandem repeats of an MC graph with abPOA and splices them back in with exact paths.
- GitHub issues: smoothxg #6, #98, #204; pggb #26, #462.
- The PGGB preprint (bioRxiv 10.1101/2023.04.05.535718, the same work as Nature Methods 2024 21:2008), smoothxg sections: main text and Supplementary B.1/B.2.

Paths below are relative to each repository root. "(inferred)" marks my reasoning rather than something the code or docs state.

---

## 1. smoothxg block decomposition

### 1.1 Prep: the order the blocks come from
- With no `-n`, prep runs `odgi chop` to 100 bp (`-X`, default 100: `src/main.cpp:70,365`). It then does a path-guided SGD 1D sort, `groom`, and a toposort (`src/prep.cpp:36-150`; `groom` at `:143`, toposort at `:147`).
- Blocks are contiguous runs of nodes in this 1D order. Block quality therefore depends on the stochastic PG-SGD sort.
  - ekg in smoothxg #6: "the prep step is stochastic ... so smoothxg is inherently inconsistent".
  - pggb forces `resume=false` after smoothxg because "smoothxg is not deterministic" (`pggb:732`).

### 1.2 Greedy block growth: `smoothable_blocks` (`src/blocks.cpp:7-327`)
- Every path step has a `seen` bit, one bitvector per path (`blocks.cpp:21-42`). A step belongs to exactly one block. This is what keeps the partition of every path exact and gapless.
- smoothxg walks handles in sort order (`:240`). For each handle it computes three things:
  - `sequence_to_add`: the handle length times its unseen steps (`:245-252`).
  - An estimate of the longest path inside the block (`:255-262`).
  - The longest edge jump in sort-order coordinates (`:264-282`).
- It finalizes the current block and starts a new one when any of these holds (`:284-299`):
  - `total_path_length + sequence_to_add > max_block_weight`. This is `-w`, which defaults to `poa_target_length * n_haps` (`main.cpp:377`).
  - `longest_edge_jump > max_edge_jump`. This is `-e`, default 0 = off.
  - `max_path_length > max_block_path_length`. This is the target POA length `-l` (`main.cpp:375`).
- `finalize_block` (`:106-230`) turns the block's handles into per-path ranges:
  - It collects every unseen step on the block's handles and sorts them by (path, rank) (`:109-127`).
  - It starts a new range when the path changes, or when the path position jumps by more than `max_path_jump` (`-j`). The jump is measured in path coordinates (`:134-152`).
  - It breaks ranges at steps already claimed by earlier blocks (`:154-177`), drops empty ranges (`:180-186`), marks the remaining steps seen and records lengths (`:189-201`).
  - It sorts ranges longest first (`:206-219`; `order_paths_from_longest` is hard-coded true at `main.cpp:363`).
  - It splits the block into the weakly connected components that its path ranges induce (`toposplit_block`, `:43-105`, called at `:224`).
- **A path that leaves the block and comes back yields several ranges of the same path in one block.** This is the mechanism that turns tandem copies into loops (§4).
- PGGB runs with `-j 0 -e 0` (`pggb:36-37,743-744`). Any path jump at all therefore starts a new range.

### 1.3 Breaking long blocks: `break_blocks` (`src/breaks.cpp:106-330`)
- This runs if any range is longer than `max_poa_length` (`-q`, default `2*target`, `main.cpp:376`) and the block has more than one range (`:211-219`).
- Repeat detection (`:225-271`):
  - It runs `sautocorr::repeat` on each range of at least `2*min_copy_length` bp. Its parameters: lags `min_copy_length..max_copy_length` (`-c` 1000, `-W` 20000, `main.cpp:285-286`), target z = 5, stride 50 (`main.cpp:457-458`).
  - `sautocorr` autocorrelates the raw ASCII byte values. It returns the last lag of the first run of z-scores above target, not the argmax (`sautocorr/sautocorr.cpp:47-69`).
  - If any repeat is found, `cut_length = mean(repeat lengths)/2` (`breaks.cpp:262-266`). Otherwise `cut_length = max_poa_length`.
- Chopping (`:272-310`): each range is cut independently, every `cut_length` bp of its own sequence, at node boundaries.
  - **All the pieces stay in the same block** (`block.path_ranges = chopped_ranges`, `:311`), re-sorted longest first (`:314-327`).
  - The cut points are not homologous across haplotypes (inferred from the code: each path is cut by its own offset).
  - When a repeat is found, even short ranges are chopped (`:274`).
- Splitting by identity (`:335-586`):
  - It deduplicates sequences in either orientation (`:340-371`).
  - Above a depth threshold it sorts them by (length, lexicographic) and greedily clusters them (`:373-528`). It tests each sequence against groups from the last one back. A length-ratio test (`-R`) stops early.
  - The identity test is either a mash estimate (`-E`/`-L`/`-D`, k=17) or a WFA gap-compressed identity (`-I`). The WFA penalties are 0/7/11/1 (`:198-203`) and the identity function is at `:72-102`.
  - Each group becomes its own block (`:536-566`).
- **Block splitting is disabled under PGGB defaults.** The gate is `min_dedup_depth_for_block_splitting != 0 && ...` (`breaks.cpp:373`). smoothxg's default is 0 = disabled (`main.cpp:100,313-314`), and `pggb` passes `-d 0 -D 0` (`pggb:749`).
  - So `-I` (set from wfmash `-p`/100 = 0.90, `pggb:486,741`) has no effect in current PGGB.
  - smoothxg #98 (2021) recommends `-I 0.7-0.9` to keep blocks homogeneous. Back then splitting was evidently active (inferred; the `-d` gate came later).

### 1.4 Multiple passes
- `-l` takes a comma list, and each value is one complete smoothxg pass (`main.cpp:287-292,374-1065`). PGGB defaults to `-G 700,1100` (`pggb:39`); the preprint describes 3 passes.
- A pass ends by writing a graph that the next pass re-preps and re-sorts.
- The stated purpose is to soften edge effects at block boundaries (PGGB preprint, main text).

### 1.5 Boundary consistency: summary
- Within one pass, boundaries are node boundaries of the chopped and sorted input graph.
- Every path is partitioned exactly once, guaranteed by the seen bits.
- Paths that run through the same node run start or stop at the same node. Paths that diverge do not (inferred).
- There is no cross-block node sharing. Adjacent blocks are separate POA graphs, joined only by the edges the paths imply (§3).

---

## 2. smoothxg POA usage

### 2.1 Engine
- The default is SPOA. `-A` selects abPOA (`main.cpp:130-132`; dispatch at `smooth.cpp:2075-2119`).
- Mode: local by default; `-Z` makes it global (`main.cpp:133-135,487`).
  - SPOA: `kSW` when local, else `kNW` (`smooth.cpp:752-755`).
  - abPOA: `ABPOA_LOCAL_MODE` / `ABPOA_GLOBAL_MODE` (`smooth.cpp:259-263`).
  - The call site has the comment "abPOA local mode is buggy, so when we use abPOA instead of SPOA, go global" (`smooth.cpp:2086`). The code still passes `local_alignment` unchanged: the comment and the code disagree.

### 2.2 Scores
- Default `-p 1,4,6,2,26,1`: match, mismatch, gap1 open/extend, gap2 open/extend. This is a convex two-piece gap (`main.cpp:322-361`).
- PGGB maps `-P asm5/asm10/asm15/asm20` to minimap2-style sets, with asm20 `1,4,6,2,26,1` as the default (`pggb:464-484`).
- The preprint describes the low-divergence set `1,19,39,3,81,1`.
- SPOA gets the penalties negated (`smooth.cpp:2101-2106`).
- With 4 parameters, SPOA copies gap1 into gap2 (affine), while abPOA sets gap2 = 0 (`main.cpp:348-360`).
- `-a` adaptive scores (`smooth.cpp:1982-2069`):
  - Estimates mash identities all-vs-all over sequences of at least `8*k` bp.
  - Takes the 30th percentile, floored at 0.7, and picks one of five score sets from asm5 to asm20.
  - Skipped above the padding depth limit.

### 2.3 abPOA parameters (`smooth.cpp:255-297`)
- `wb=311`: smoothxg always runs abPOA banded (`smooth.cpp:2090`, and `:266-270`).
- `wf=0.03`.
- `amb_strand=0`: sequences are pre-oriented.
- `disable_seeding=1`.
- `k=19, w=10, min_w=3313`. These are unused because seeding is off.
- `out_gfa=1`; MSA output only if MAF is requested.
- Per-base weights equal the dedup counts (`:332-337`).

### 2.4 Building each POA input (`smooth.cpp:177-244`, `:676-743`)
- Each path range is spelled with `poa_padding` bp of the path's own flanking sequence on each side (`append_to_sequence`, `smooth.cpp:75-126`). When the path ends, it is filled with `N`.
- The sequence is reverse-complemented if most of its bases are on reverse handles (`:208-210`).
- Exact duplicates, by xxHash of the padded sequence, become one sequence with weight = count. Names and strands are kept per duplicate (`:216-241`).
- Order: longest first, from the block sort.
- SPOA aligns the sequences one by one in that order and adds each with its weight (`smooth.cpp:759-769`). abPOA gets them all at once through `abpoa_poa` (`:337`), which is non-progressive (`progressive_poa` is commented out at `:294`).

### 2.5 Padding: the "poa padding" trick (`smooth.cpp:1946-1970`)
- `poa_padding = 311` when the block depth is at most `-Y`. Above `-Y` it is 0.
- It is then raised to `avg_seq_len * -O` (default 0.001) if that is larger.
- PGGB passes `-Y = 100*n_haps` (`pggb:748`). In practice every non-pathological block gets 311 bp of real flank on each side.
- After POA, the padding is cut from each path's node list: the first and last `padding_len` 1-bp POA nodes are dropped.
  - abPOA: `smooth.cpp:2520-2531`. SPOA: `:2611`.
  - MSA rows are trimmed to match (`:365-422`, `:790-847`).
- Then POA nodes and edges with path depth 0 are destroyed (`find_*_exceeding_depth_limits(…,1,max)`, `smooth.cpp:2559-2573,2639-2653`). This removes nodes that only padding used.
- Purpose, as the preprint puts it: "padding extends each block to the left and right", improving alignment at block edges. The extra context anchors the ends of the real block sequence, especially under local SW (inferred).

### 2.6 Consensus
- If consensus is not vanished (`-V`), a `Consensus_<block>` path is added.
  - SPOA: `GenerateConsensus()`, heaviest bundle (`smooth.cpp:771-780`).
  - abPOA: `abpoa_generate_consensus`, restricted to nodes that carry a real path step (`:346-351`, `:2535-2551`).
- The SPOA consensus path is not trimmed of padding (`smooth.cpp:2624-2637`). Its padding nodes therefore survive the depth-0 cleanup (inferred from the code).
- PGGB passes `-V` by default because consensus graphs are disabled over pggb #133/#182 (`pggb:62,706-710`).
- Consensus graphs (`-C`, `consensus_graph.cpp`) are a separate, optional product.

---

## 3. smoothxg lacing: rebuilding the graph with exact paths

### 3.1 Per-block graph (`smooth.cpp:538-626`)
1. The POA graph becomes an odgi graph with one 1-bp node per POA node (`build_odgi_*`).
2. The graph is unchopped (`:545`) and topologically ordered (`:557`).
3. The result is copied with **edges only from consecutive steps of the embedded paths**. The commented-out alternative "there were edges not supported by any path!" is at `:567-604`.
4. Paths are re-created in the original `path_ranges` order, with the consensus last (`:606-620`).
- The path name in a block is `<path>_<start>`.
- smoothxg records `{path, start_pos, end_pos, local path ordinal, block_id}` in an mmap multiset (`smooth.cpp:2283-2296`). Block graphs are stored zstd-compressed (`save_block_graph`).

### 3.2 Global lace (`main.cpp:599-768`)
1. It sorts the fragments by (path, start) (`main.cpp:611`).
2. It concatenates the nodes and edges of all block graphs with an id offset per block (`:642-665`).
3. For each original path, in parallel, it walks its fragments in order and appends each fragment's steps through the offset. It adds one edge from the previous fragment's last step to this fragment's first step (`:706-757`).
4. Gaps between fragments are only `assert(false)` (`:729-731`, `:760-762`). These are compiled out in Release builds (inferred).
5. It validates **every path spells exactly its original sequence**, from a separate sequence index, and exits on mismatch (`:770-803`). It then checks the path count (`:805-810`).
6. It re-walks all paths to make sure every consecutive step pair has an edge (`:996-1018`), then runs a global `odgi unchop` (`:1021`).
- Consensus paths are embedded per block, or merged across contiguous blocks with `-M` (`main.cpp:812-986`).

### 3.3 Fallbacks for long or deep blocks
- **There is no runtime fallback.** Expensive problems are only prevented by the cuts in §1.2 and §1.3: `-l` target, `-q` max, `-w` weight.
- The old "too deep, mirror the input subgraph" fallback (`MAX_POA_BLOCK_DEPTH 4000`) is commented out (`smooth.cpp:129`, `:1015-1085`).
- A single oversized POA still OOMs. smoothxg #6: SPOA's quadratic memory on a long VNTR block; the advice was to reduce threads with `-T`.

### 3.4 Known exact-path bug
- pggb #462 (open, 2025): smoothxg drops a path. "We lose the shortest path ... when its node is the last one in the GFA ... not connected to any other nodes." Only the path-count check at `main.cpp:805` catches it.

---

## 4. Known problems on tandem repeats, VNTRs and copy-number differences

### 4.1 Mechanisms visible in the code
- **Repeat copies become separate POA inputs in one block, which collapses them into a cycle (inferred from the code).**
  - A path that loops through the block's nodes (a VNTR already folded by seqwish) yields several ranges (`blocks.cpp:134-152`). The repeat-chopping step also cuts every range into period/2 pieces kept in one block (`breaks.cpp:262-311`).
  - Each piece is aligned onto the same POA graph, so all copies share nodes and the path traverses the block several times.
  - This is by design: the preprint says POA lets the graph "represent complex structurally-variable loci as simple loops through a single copy of duplicated sequences", and ekg in smoothxg #6 calls it a parameter "which defines the size of SDs/VNTRs to maintain as cyclic in the output".
  - For plegma this is exactly the wrong representation: copy-number alleles lose their linear MSA columns.
- **The pieces are not homologous.**
  - Cuts fall every `cut_length` bp of each haplotype's own sequence, not at homologous positions.
  - The period comes from ASCII-value autocorrelation with stride 50, and it ignores periods under `-c` = 1000 bp (`main.cpp:285`, `breaks.cpp:236-245`).
  - ekg in #6: "the exact repeat length is not stable, and so finding the start and end is not easy unless we use alignment. The only reasonable solution appears to be to run the self alignment ... and break them at obvious repeat ... starts". The current code still uses autocorrelation.
- **Blind cuts at `max_poa_length`.** For short-period VNTRs (period under 1 kb) longer than `2*-l` (1.4 / 2.2 kb under PGGB), ranges are cut blindly per haplotype into ~`max_poa_length` pieces (`breaks.cpp:268-269`, `:288-309`).
  - A haplotype with 3 more copies is cut at different phases, and the pieces are aligned together anyway (inferred).
- **Large copy-number differences in one block become long gaps against the longest-first sequence.**
  - Local SW, as PGGB runs it, can leave unmatched array ends unaligned or misaligned. Padding only stabilises the outer ends (inferred).
  - With splitting disabled (§1.3), very different alleles share one POA.
- **The block length targets (700/1100 bp) are far below typical VNTR array lengths (1-30 kb).** A VNTR region is cut into many blocks whose boundaries come from the 1D sort, not from repeat structure (inferred).

### 4.2 Issues and docs
- **smoothxg #6**: SPOA OOM caused by VNTRs.
  - ekg: "if we pick up a VNTR loop in the graph, the subsequences that traverse it all get concatenated", producing very long POA inputs.
  - "Fewer cycles are present in haplotypes that have a shorter VNTR allele length."
  - "There are loops within loops ... We can't just break if the self step count of a path increases."
  - This thread is where the autocorrelation repeat breaker came from.
- **smoothxg #204** (ekg): smoothxg's goal is a *local* partial order at scale `-G`, about 1 kb. Before smoothxg, short tandem repeats are "extremely dense and complex ... a single C might represent an entire 20bp homopolymer".
  - Values much larger than ~1 kb are "computationally prohibitive" because POA is quadratic.
- **smoothxg #98**: on human data, without identity clustering (`-I 0`), POA blows memory. Smoothing raised the maximum node degree to 9351, which ekg attributes to paralogous repetitive mappings.
- **pggb #26**: users see "over-alignment, mostly in aligning sequences to themselves (i.e. loops)".
  - ekg: "The alignments may not be able to consistently linearize highly-similar repeats in of a VNTR, because the direct relationships between the various copies may be ambiguous." The suggested levers are larger `-s/-p/-w/-e/-j`, and seqwish `-k` for STR looping.
- **PGGB preprint, Supp. B.2 and C**: before smoothxg, centromeres and satellites form high-degree tangles. Raising seqwish `-k` (47, 79) reduces complexity in alpha satellites.
  - The preprint also says filtering short matches leaves the graph "locally under-aligned", which smoothxg is meant to fix.

---

## 5. How Minigraph-Cactus uses abPOA in BAR, and why MC graphs mis-align VNTRs

### 5.1 Parameters in pangenome mode (cactus HEAD)
- **Base aligner.** `--pangenome` sets `baseAligner = pangenomeBaseAligner` = `abpoa` (`src/cactus/setup/cactus_align.py:387-390`; config `src/cactus/cactus_progressive_config.xml:300`). Progressive Cactus uses `pecan`.
- **Length limits.** `bandingLimit = --maxLen = 10000` (`cactus_align.py:399-405`; `refmap/cactus_pangenome.py:114`), with `caf maxRecoverableChainLength = maxLen/2 = 5000`.
- **CAF settings** (`cactus_align.py:407-418`):
  - `minimumBlockHomologySupport=0`.
  - `maxRecoverableChainsIterations=50`, with the code comment "helps quite a bit to reduce underalignment".
  - `bar minimumBlockDegree=1`.
  - POA seeding disabled.
  - `poaReferenceEvent` = the reference.
- **abPOA setup** (`bar/impl/poaBarAligner.c:390-425`, config `:390-411`):
  - Global mode only (the comment says "only global works").
  - Adaptive band `wb=1000`, `wf=0.1`.
  - A 5x5 matrix: `91 -114 -61 -123 ...`, HOXD70-like.
  - Convex gaps 400/30 and 1200/1.
  - Seeding off in pangenome mode.
  - Progressive guide tree on (Jaccard), disabled above 5000 rows.
  - `out_cons=0`.
  - Scores may be replaced by a `last-train` model (`src/cactus/paf/last_scoring.py:163-247`).
- **Sliding window.** `partialOrderAlignmentWindow=10000`, with a hard-coded 50% overlap (`poaBarAligner.c:1017-1023`); 5000 bp for genomes of at least 15 Gb.

### 5.2 How BAR chunks the sequence
- BAR aligns **adjacencies**: the sequences between CAF-aligned blocks, at each "end" of a flower.
- Each thread's adjacency is taken as a prefix of at most `bandingLimit` (10 kb) from the end (`get_adjacency_string_and_overlap`, `poaBarAligner.c:1473-1531`).
- If one dominant end covers all adjacencies and all are under 10 kb, it builds one MSA (`:1735-1760`).
- Otherwise it aligns each end separately: the left end sees 10 kb prefixes and the right end sees 10 kb prefixes of the reverse complements.
  - Where the two overlap, it trims both at the cut point that maximises summed column scores (`make_consistent_partial_order_alignments`, `:1336-1373`; `trim`, `:227-282`).
  - **The middle of an adjacency longer than 2x10 kb is never aligned** (inferred from the code).
- Within an MSA, the sliding window loads **up to `window_size` bases of each row from that row's own offset** (`:1075-1081`).
  - It aligns, then trims consecutive windows at the best column-score cut inside the overlap (`:1097-1150`), and concatenates (`:1170-1200`).
  - Empty rows are faked as one `N` (`:1083-1095`).
- Threads are sorted longest first (`get_end_sequences`, `:1684-1715`).
- The MSA becomes maximal gapless blocks, which become pinches back into the cactus graph (`create_alignment_blocks`, `:1533-1665`; pinch iterator `:1881-1928`).
- **POA orientation.** `poaOrientation="auto"` aligns so the reference reads forward. The config says the reason: abPOA/minipoa "put a gap that could go anywhere along a repeat at its left end, so this picks which end of the repeat it lands on" (`cactus_progressive_config.xml:294`).

### 5.3 Why MC mis-aligns VNTRs
Documented and inferred causes:
1. **Minigraph backbone.** Minigraph keeps only variation of at least 50 bp (PGGB preprint, App. E).
   - It "stores an allele as novel sequence whenever it cannot align it ... even when the allele is a copy of sequence the graph already has ... Cactus cannot merge two separate minigraph nodes later, so such an allele stays a giant insertion" (`doc/pangenome.md:192`).
   - Expanded VNTR alleles therefore often enter as unaligned insertions next to the reference copies. `rgfa-zip` only rescues confident (at least 5 kb, 95%) unique placements, which tandem copies rarely give (inferred).
2. **Anchors.** CAF anchors come from assembly-to-minigraph mappings. Inside arrays these mappings are ambiguous, so pinches can tie different copies together: cycles in non-reference paths, or shifted phase (inferred).
   - The reference is forced acyclic, never self-aligned (`doc/pangenome.md:134-139`). Self-alignment and cycles are only an option with `--collapse`.
3. **Window homology.** BAR's 10 kb windows are taken per row by offset. Haplotypes whose arrays differ by kilobases put non-homologous content into the same window, and the stitch is a heuristic cut (inferred from `poaBarAligner.c:1075-1081`).
   - Arrays longer than 20 kb between anchors are left unaligned in the middle.
4. **Gap placement.** POA left-places gaps in repeats (config comment above). Different flowers or orientations can place homologous indels differently (inferred, and the motive for `poaOrientation`).
5. **Clipping.** MC clips sequence not aligned to the minigraph paths, removing centromeres and satellites from the clip graphs (PGGB preprint App. E; `doc/pangenome.md:1463`).
6. **Folded repeats.** pgrealign's QC text notes that a graph "that has folded copies of a repeat into a cycle" can align each base with many others, and scores that as lost precision (`pgrealign/docs/method.md:769-778`).
   - Cactus now ships `--pgrealign` [EXPERIMENTAL] specifically to "realign each chromosome's tandem repeats" (`refmap/cactus_graphmap_join.py:430-441`; `doc/pangenome.md:567-606`). That is a direct admission that MC's BAR output in tandem repeats is unsatisfactory.

### 5.4 pgrealign, the closest competitor to plegma (brief)
- **Regions** (`docs/method.md:264-314`):
  - Each target interval is padded and anchored on two snarl boundary nodes of **one chain**, visited once and forward by the reference. Every walk that crosses the region must then pass both anchors.
  - Overlapping regions are merged, and `--max-span` caps the span.
- **Extract** (`:316-410`). Walk runs are classified as spanning, prefix, suffix, internal or complex.
  - Nodes visited in both orientations make the region `complex`, so inversions are skipped.
  - Distinct spanning sequences become weighted alleles, and partial ones become fragments with a side.
- **Realign** (`:412-648`), a medoid star with abPOA `-m 0 -r 1`, global:
  - The centre is chosen by weighted multiset-15-mer Jaccard. Copy number counts.
  - A homology gate (`--min-share` 0.5 distinct-15-mer share) marks unrelated or inverted alleles `not_homologous`.
  - Each allele is aligned pairwise to the centre, unbanded unless over `--max-cells`.
  - Indels are left-normalised, and the pairwise alignments are merged on centre columns. Each insertion slot's distinct inserts are aligned to each other, banded.
  - Fragments are placed exactly when they are substrings of an allele.
  - Scores are abPOA defaults x4, with E2 kept at 1, to favour one long gap over chance-match threading.
  - Region statuses: `too_big`, `timeout`, `aligner_failed`, `not_homologous` and others. Each keeps its old subgraph.
- **Replace** (`:650-724`):
  - Builds an induced graph with one node per distinct base per column, then unchops.
  - Numbers nodes topologically, column by column, and chops to 1024 bp or less.
  - Swaps out the interior and rewrites each run.
  - Checks every rewritten run spells the same upper-cased sequence. `verify` then compares every path.
- **Known gaps**:
  - No P-line rewrite, and lowercase (soft-masking) is lost.
  - Pairs of non-centre alleles are only projected through the centre.
  - abPOA reserves the full DP matrix even when banded: over 100 GB for two 100 kb alleles.

---

## 6. Advantages and disadvantages, judged against plegma's goal

Plegma's goal: correctly align VNTR alleles with copy-number differences across 34-456 haplotypes of 1-30 kb, and rebuild an exact-path graph.

### smoothxg: advantages
- Exact-path invariant enforced end to end: a full sequence check of every path and a path-count check (`main.cpp:770-810`).
- Clean separation: blocks become sequences, then POA, then a per-block graph, then lacing by recorded (path, start, end, block) fragments. Blocks are embarrassingly parallel.
- Deduplication with weights cuts POA work and steers consensus (`smooth.cpp:216-241`).
- Real-flank padding, then trimming, to stabilise POA ends (`smooth.cpp:75-126,1946-1970`).
- Graph edges come only from path steps. POA edges that no path uses, and padding-only nodes, are dropped (`smooth.cpp:589-604,2559-2573`).
- Optional adaptive scoring from estimated identity (`smooth.cpp:1982-2069`).
- Length-ratio and identity pre-clustering exists (`breaks.cpp:335-586`), though PGGB disables it.

### smoothxg: disadvantages
- Repeats are deliberately folded into cycles by putting multiple copies or pieces of one path into one POA (`blocks.cpp:134-152`, `breaks.cpp:262-311`). That is the opposite of a copy-number-resolving MSA.
- Repeat boundaries come from crude autocorrelation (ASCII values, stride 50, periods of at least 1 kb only). Cuts are per-haplotype offsets, not homologous positions.
- Block sizes of 700-1100 bp, with boundaries from a stochastic 1D sort, are much shorter than the arrays. They are not anchored on sequence homology, and runs are not reproducible.
- Local SW by default. Splitting is disabled in PGGB, so divergent alleles share one POA.
- There is no runtime fallback for expensive POAs (OOM, #6). Gap detection relies on `assert`. A path can be dropped (pggb #462).
- smoothxg operates on the whole graph and has no "take one region" interface.

### Minigraph-Cactus BAR: advantages
- Global abPOA with convex gaps, an adaptive band and a progressive guide tree. Pangenome-tuned (optionally trained) scores.
- The reference-orientation rule puts repeat gaps consistently on one side.
- Window stitching by column-score-optimal cut points (`poaBarAligner.c:227-282`). This is a reusable idea for overlapping windows.
- Results feed back as pinches into one consistent graph. The acyclic reference gives stable coordinates.

### Minigraph-Cactus BAR: disadvantages
- It only aligns between CAF anchors derived from minigraph mappings, which are ambiguous in arrays. Novel minigraph insertions can never be merged with existing copies (`doc/pangenome.md:192`).
- Only 10 kb is aligned from each end of an adjacency, so arrays over ~20 kb are partly unaligned. Windows are per-row offsets, not homologous.
- Tandem duplications are insertions (no `--collapse`) or ambiguous pinches. The reference is acyclic but other haplotypes can cycle.
- Tandem-repeat realignment was bolted on afterwards (pgrealign), which shows the core pipeline does not handle them well.

### pgrealign (for context): advantages and disadvantages
- **Advantages:**
  - Exactly plegma's interface: region, then alleles and fragments, then MSA, then splice with exact paths.
  - Copy-number-aware centre choice and a homology gate.
  - Explicit statuses that keep the old subgraph.
  - QC against optimal pairwise alignments.
- **Disadvantages:**
  - A star alignment only optimises centre-to-allele pairs.
  - Insertion slots relative to the centre (extra copies) are aligned to each other banded and in isolation.
  - Pairwise abPOA on 10-30 kb alleles costs quadratic memory and time unless banded.
  - Python plus subprocess abPOA.
  - Region anchors depend on vg snarls of the existing (possibly mis-built) graph.

---

## 7. Ideas to steal and pitfalls to avoid

### Ideas to steal
1. **A region interface like pgrealign's.**
   - Anchor each target on two boundary nodes of one snarl chain, visited once and forward by the reference. This guarantees every crossing walk passes both.
   - Pad the target into flank so the anchors sit outside the array.
   - Merge interfering regions, and record a status per region.
   - Classify walk runs: spanning, prefix, suffix, internal, complex. Treat nodes visited both ways as complex and skip them.
   - Make plegma's region extractor emit exactly the alleles (distinct spanning sequences with weights) and fragments (with side) that pgrealign uses. A plegma output can then drop into pgrealign's `replace`/`verify`, or into cactus `--pgrealign`, for head-to-head QC.
2. **Deduplicate with weights before any DP** (smoothxg, pgrealign). With 34-456 haplotypes there are usually far fewer distinct alleles.
3. **Copy-number-aware similarity.** Use multiset k-mer Jaccard (counts, not sets) for centre or guide-tree choice, with a homology gate (distinct-k-mer share, at least 0.5) to refuse unrelated or inverted alleles rather than thread them.
4. **Padding and anchor overlap.** Give each POA window real flanking sequence and trim it afterwards (smoothxg `smooth.cpp:75-126,2520-2531`). For plegma (inferred): extend each DP window into the flanking MUM/SMEM anchors, so window ends are forced through shared sequence. Then cut at the anchor midpoint, which is homologous in every row, unlike MC's per-row offsets.
5. **Overlapping windows stitched at a column-score-optimal cut** (cactus `trim`, `poaBarAligner.c:227-282`). Only stitch windows whose rows are homologous, which plegma's anchor partitioning should guarantee.
6. **Graph from MSA.**
   - One node per distinct base per column. Edges only from consecutive row bases.
   - Unchop under pgrealign's rule: merge u into v only if u has one successor, v has one predecessor, and no walk ends at u or starts at v.
   - Topological numbering, then a node chop of 1024 bp or less.
   - Drop anything no path uses (smoothxg `:2559-2573`; pgrealign `induce.py`).
7. **Lace by recorded fragments.** Keep (path, start, end, block or region id, orientation) for every replaced run. Rewrite reverse runs flipped.
8. **Verify everything.**
   - After every replace, check that every path spells its original sequence and that the path count is unchanged (smoothxg `main.cpp:770-810`; pgrealign verify).
   - Treat a failure as an input or logic error that stops the run, not as one hard region.
9. **Fall back by keeping the old subgraph** per region (`too_big`, `timeout`, `aligner_failed`, `not_homologous`, `complex`). Estimate memory before running. Never let one locus kill a contig.
10. **Left-normalise indels** in repeats (pgrealign `method.md:536-550`). Fix one gap orientation rule, such as reference-forward, as MC's `poaOrientation` does.
11. **Scoring for VNTRs.** Make long gaps relatively cheap so a copy-number difference becomes one long gap, not chance-match threading. pgrealign uses abPOA defaults x4 with E2=1; MC uses 400/30 and 1200/1. A convex (two-piece) gap model is what both converge on.
12. **QC harness.** For sampled allele pairs, measure the projected-alignment precision and recall against optimal pairwise alignments, before and after. Also measure `growth` (graph bases over longest allele) and forced-homology windows (pgrealign `method.md:738-811`). This is a ready-made scorecard for plegma vs MC vs pgrealign.

### Pitfalls to avoid
1. **Never put several pieces of the same haplotype into one POA problem** unless a cycle is wanted. That is how smoothxg folds VNTRs into loops. Plegma wants linear copy-number columns.
2. **Do not cut haplotypes at per-row offsets** (smoothxg `breaks.cpp:288-309`, MC `poaBarAligner.c:1075-1081`). Cut only at homologous anchors (MUM/SMEM) shared by all rows of a partition. If an anchor is missing in some rows, recurse or merge partitions; do not fall back to blind chopping.
3. **Do not estimate the period with crude autocorrelation.** If periodicity matters, get it from self-alignment or SMEM self-hits, as ekg himself concluded in #6.
4. **Do not rely on a stochastic graph sort for boundaries.** Derive them from sequence anchors, so results are deterministic and reproducible.
5. **Do not use local alignment for block interiors** when the ends are anchored. Use global or semi-global between anchors; smoothxg's local SW plus padding is a workaround.
6. **Do not aim a star alignment at large copy-number differences without care.** Non-centre pairs are only projected, and extra copies become separately aligned insertion slots. Plegma's progressive POA in windows should beat this on non-centre pairs; measure it with pgrealign's qc.
7. **Watch memory.** abPOA reserves the full matrix even when banded (pgrealign known gaps), and SPOA OOMs on long VNTR blocks (smoothxg #6). Plegma's windowed DP must bound cells per window and per thread, and check the bound before allocating.
8. **Path-handling edge cases.** Check them explicitly:
   - Isolated single-node paths (pggb #462).
   - Walks that start or end inside a region (fragments).
   - Walks that never touch an anchor (internal).
   - P lines vs W lines.
   - Lowercase or soft-mask preservation, which pgrealign loses.
   - Reverse-stored nodes.
   - Do not rely on `assert` for invariants; smoothxg's lacing gap checks vanish in Release builds.
9. **Keep consensus paths out of the output** unless asked. Strip padding from any consensus: smoothxg's SPOA consensus keeps padding nodes (inferred, `smooth.cpp:2624-2637`).
10. **Do not let one bad region abort the run** (smoothxg exits the whole run on one unchop or POA failure, `smooth.cpp:545-554`). Use per-region statuses instead.

---

### Key file index
- smoothxg:
  - `src/main.cpp`: options at 54-186, defaults at 274-377, lace at 599-768, validation at 770-810.
  - `src/blocks.cpp:7-327`.
  - `src/breaks.cpp:106-614`.
  - `src/smooth.cpp`: padding at 75-126, abPOA at 133-627, SPOA at 629-1089, per-block driver at 1930-2312, graph builders at 2442-2654.
  - `src/prep.cpp`.
- sautocorr: `sautocorr.cpp:5-82`.
- pggb: `pggb:36-54,464-486,734-755`.
- cactus:
  - `bar/impl/poaBarAligner.c`: params at 390-425, window at 995-1200, ends at 1336-1373, adjacency at 1473-1531, flower at 1735-1880.
  - `src/cactus/cactus_progressive_config.xml:258-411`.
  - `src/cactus/setup/cactus_align.py:379-421`.
  - `src/cactus/refmap/cactus_graphmap_join.py:430-441,1788-1920`.
  - `doc/pangenome.md:132-196,567-606`.
- pgrealign: `README.md`; `docs/method.md:264-869`.
- Issues: smoothxg #6, #98, #204; pggb #26, #462.
- Paper: [Building pangenome graphs (bioRxiv)](https://www.biorxiv.org/content/10.1101/2023.04.05.535718.full); [PubMed 37066137](https://pubmed.ncbi.nlm.nih.gov/37066137/).
