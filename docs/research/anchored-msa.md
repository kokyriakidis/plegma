# Anchor-based and repeat-aware alignment: prior art for plegma

Scope: prior art for plegma, a Mojo partial-order aligner (POA) for re-aligning VNTR regions of the HPRC Minigraph-Cactus graph. A region has 34-456 haplotypes, typically 1-30 kb long.

The design as corrected by the user:
- The POA graph is built **progressively, one sequence at a time**.
- Each new sequence is anchored **pairwise** against the current graph or its paths.
- Anchors are sporos pairwise MUMs, falling back to the largest SMEMs.
- An **optimal backbone chain** is found first. Inside long gaps the chain is recomputed with smaller anchors.
- The gaps between chained anchors are then filled with banded DP from karan.

Multi-MUMs (matches unique in every sequence) are out of scope, except where they inform pairwise use.

How to read this report:
- **Code reading** means I read the source in a local clone under `scratchpad/src/<name>`, given as relative paths.
- **Inferred** marks my own reasoning or a claim I did not check in code or paper.
- DOIs are given where I am confident of them. "DOI unverified" means I am not.

---

## 0. Summary

1. **Unanchored POA is already the accuracy baseline on these regions, and anchoring has made it worse.** The local evaluation repo `vntr-eval` scored candidates by cost/opt, the summed cost of the pairwise alignments the graph implies divided by the pairwise optima (1.0 is optimal).
   - abPOA with default flags, longest sequence first, reaches a median of 1.067 at hotspot VNTRs. Minigraph-Cactus (MC) scores 1.335.
   - abPOA `-S` (minimizer anchoring and partitioning) scores 1.894. Its MSA has 92,688 columns where default abPOA has 13,051, because haplotypes are left side by side.
   - Source: `vntr-eval/tools/POA.md`.
   - So **anchors must never make a decision that full DP would get right**. Their job is to bound memory and time, and to fix placement where real unique sequence exists.
2. **abPOA `-S` fails for three visible reasons** (`abpoa/src/abpoa_seed.c`):
   - Every pair of identical minimizers becomes an anchor, with no frequency filter.
   - Chaining refuses diagonal shifts over 100 bp (`max_bw = 100`) and predecessors more than 100 bp away (`max_dis = 100`).
   - It anchors against the **previous sequence in the guide order**, not against the graph.

   Each of these is fatal in VNTRs, where consecutive true anchors differ in diagonal by k·period. This mechanism is inferred from the code; the failure itself is measured.
3. **FMAlign2 is the closest published relative of plegma's idea**, and its code shows the gaps (`FMAlign2/src/mem_finder.cpp`):
   - It takes one level of MEM partitioning across all sequences.
   - Anchors are ordered by their mean position.
   - Each sequence picks the occurrence closest to that mean, which is a naive repeat rule.
   - There is no recursion. The segments go to MAFFT or HAlign.
4. **centrolign (Eizenga et al. 2026) is the closest relative for the graph and repeat parts.** It is a progressive graph-to-graph aligner for tandem arrays (`centrolign/include/centrolign/*.hpp`). It already does what plegma needs:
   - It indexes the **embedded path sequences** of each graph.
   - It keeps "minimal rare matches" occurring at most `max_count = 50` times in either graph.
   - It weights an anchor by `len / sqrt(count1·count2)` minus a concave length penalty.
   - It chains with a **3-piece (concave) affine gap** over DAG reachability, using a chain decomposition.
   - It **re-anchors between chained anchors** (`do_fill_in_anchoring`).
5. **Theory gives a clean guarantee for overlap-aware chaining, but only relative to the anchor set** (Jain, Gibney, Thankachan, JCB 2022, Theorem 2):
   - With their gap+overlap cost, the optimal chain cost equals the **anchored edit distance**, where a match not supported by an anchor costs 1.
   - If the anchor set is **all MEM occurrences ≥ L**: ED ≤ anchored ED ≤ ED + (matched bases of an optimal alignment that lie in exact runs shorter than L). This bound is inferred.
   - With MUMs only, or with SMEMs (which drop shorter on-diagonal MEMs covered by a longer off-diagonal one, typical in arrays), no such bound holds.
6. **What this means for plegma.** Uniqueness against the graph must be counted in **distinct graph locations (nodes)**, not text occurrences, as centrolign's colour-set-size counting does. Then:
   - Chain **all rare MEM occurrences** (sporos `-maxmatch --max-occ N`), weighted by rarity, together with MUMs.
   - Use overlap-aware precedence and a concave or period-aware gap cost.
   - Break ties between co-optimal chains with a deterministic left-normalisation, so equivalent copy-number gaps stack across haplotypes.
   - Recurse only in gaps larger than the DP budget, recomputing uniqueness **within the window**.
   - Fall back to full banded DP (or unit-level DP) inside arrays without anchors. Never force anchors there.
7. **Order matters as much as anchoring.** In vntr-eval, random input orders moved cost/opt by 0.11-0.33. Longest-first was best in 5 of 6 region and tool pairs. Inferred explanation: copy-number gaps become deletions against a graph that already holds all copies.

---

## 1. Empirical context from the local evaluation repos

Two local repos already study this problem, and the plegma plan should reuse them:
- `vntr-eval`: 149 region packages, the hap32 and full panels, MC baselines, and the metrics cost/opt and pair F1.
- `pgrealign` (benedictpaten): a medoid-star realigner on abPOA.

Findings relevant to anchoring:

| finding | number | source |
|---|---|---|
| abPOA default, longest first, hotspot VNTR cost/opt (median) | 1.067; MC 1.335; mafft L-INS-i 1.081 | `vntr-eval/tools/POA.md` |
| abPOA `-S` (minimizer anchoring) | 1.894 median; MSA columns 92,688 vs 13,051 at L012184 | same |
| Order sensitivity (random orders) | cost/opt spread 0.11-0.19 (abPOA), 0.23-0.33 (spoa); longest-first best in 5/6 | same |
| abPOA memory without `-S` | 18-29 B × max_len², about 9 GB at 20 kb | same |
| Unit-aware two-level aligner (units, then bases) | hotspot VNTR 1.048, 33/40 ≤ 1.1 | `vntr-eval/tools/UNIT_ALIGN.md` |
| Hard set, 45 chr20 regions with pair F1 < 0.8 | progressive abPOA 0.621 pair F1 vs medoid star 0.703; centrolign-to-medoid star 0.659; "no aligner rescues the hard set" | `vntr-eval/hardset/README.md` |
| Gap convention matters | star_gapnorm (left-normalised gaps) +0.014 pair F1 | same |

Lessons:
- The thing to beat is full-DP POA with a good order, not MC.
- Anchors that cut through arrays are actively harmful.
- Gap placement conventions are measurable in pair F1.

---

## 2. Survey of prior art

For each tool: **A** = anchors, **C** = consistency or chaining, **R** = recursion, **W** = window alignment, **Rep** = repeat handling.

### 2.1 Suffix-array / FM-index partitioning of MSAs

**FMAlign** (Liu, Zou, Xu, *Brief Bioinform* 2022, 10.1093/bib/bbab519; code `FMAlign/`, read)
- **A:** Minimizer-sampled seeds from one sequence, searched with an FM-index (FMtree) over the concatenated input. `common_seed.cpp` keeps only seeds present in all sequences (`if (seqsn < njob - 1) flag = false;`) and applies a frequency threshold.
- **C:** Seeds are chained into a common colinear set (`chain.cpp`).
- **R:** None.
- **W:** The segments between seeds go to MAFFT or HAlign.
- **Rep:** Only the frequency filter.

**FMAlign2** (Zhang et al., *Bioinformatics* 2024, 10.1093/bioinformatics/btae014; code `FMAlign2/`, read in full for the MEM path)
- **A:** `gsacak` builds a generalized SA, LCP and document array over all sequences.
  - `get_lcp_intervals` takes maximal runs with LCP ≥ L that contain an LCP equal to L.
  - It requires `right-left+1 >= min_cross_sequence`. That counts **suffixes, not distinct sequences**, so a repeat inside one sequence can satisfy it. Inferred from the code.
  - `interval2mem` left-extends while all occurrences agree, giving a multi-MEM maximal over that occurrence set.
  - MEMs that are more than 80% gap characters are discarded.
- **C:**
  - MEMs are sorted by **average position** over their occurrences (`sort_mem`).
  - In "accurate" mode each sequence keeps, per MEM, the occurrence **closest to the average position**. It then runs an O(m²) weighted non-overlap DP over the avg-pos order (`find_optimal_chain`).
  - Columns missing in any sequence are dropped. `min_seq_coverage` is forced to 1 in `find_mem`.
  - "Fast" mode does one DP on `avg_pos` with weight `len × occurrences`.
  - Small bug: backtracking `while (end_index > 0)` never emits anchor 0.
- **R:** None.
- **W:**
  - Gap segments between anchors are written to temporary FASTA files and aligned in parallel by an external MSA (MAFFT, ClustalO, HAlign3/4; `parallel_align`).
  - Anchors are pasted in ungapped.
  - `expand_chain` would Smith-Waterman a missing anchor into the gap of a sequence lacking it (partial anchors). With coverage forced to 1 this path is dead code. Inferred.
- **Rep:** "Closest to mean position" is the whole repeat model. It is wrong exactly when copy number differs, because the mean position of a copy shifts with copy number.
- Relevance: the same top-level idea as plegma, with none of the safeguards plegma needs.

**HAlign 3/4** (HAlign 3: Tang et al., *Mol Biol Evol* 2022, 10.1093/molbev/msac166; HAlign 4: Zhou, Zhang, Zou, Han, *Bioinformatics* 2024, 10.1093/bioinformatics/btae718)
- **A:** Center-star alignment. Exact matches between each sequence and the centre come from a suffix tree (HAlign 3), replaced by a BWT/FM-index in HAlign 4.
- **C:** Pairwise only, sequence vs centre.
- **R:** None reported.
- **W:** DP between matches. HAlign 4 uses the wavefront algorithm (WFA).
- **Rep:** None specific.
- Relevance: plegma's pairwise design is "center-star with a growing graph as centre". The medoid star in pgrealign is the same family.

**WMSA / WMSA 2** (Wei et al., *Bioinformatics* 2022, 10.1093/bioinformatics/btac658, DOI unverified; WMSA 2, *Brief Bioinform* 2023, bbad190)
- CD-HIT-style clustering, a center-star alignment within each cluster, and progressive profile alignment between clusters.
- No repeat model. Mentioned for completeness.

### 2.2 Whole-genome anchor aligners

**MUMmer / nucmer** (Delcher et al. 1999, *NAR*, 10.1093/nar/27.11.2369; MUMmer4: Marçais et al. 2018, *PLoS Comput Biol*, 10.1371/journal.pcbi.1005944)
- **A:** MUMs (`-mum`), MEMs unique in the reference (`-mumreference`, the default) or all MEMs (`-maxmatch`). sporos `mem` reproduces these exactly (`~/sporos/README.md`).
- **C:**
  - MUMmer 1 chains MUMs with a longest-increasing-subsequence (LIS) chain.
  - nucmer clusters anchors by diagonal difference and gap, then joins clusters.
  - Defaults recalled from memory: `-g 90` max gap, `-c 65` min cluster, `-b 200` break length.
- **R:** None as such. The gaps between anchors are extended by banded DP ("extend"), which is exactly plegma's fill step.
- **Rep:** `-mumreference` and `-mum` simply drop repeated matches. Repeats come out as separate overlapping alignments (show-coords, delta-filter).

**Mauve / progressiveMauve** (Darling et al., *Genome Res* 2004, 10.1101/gr.2289704; *PLoS One* 2010, 10.1371/journal.pone.0011147)
- **A:** Multi-MUMs (Mauve) and spaced-seed families present in subsets of genomes (progressiveMauve).
- **C:** Locally collinear blocks (LCBs) with breakpoint penalties.
- **R:** **Recursive anchoring.** Gaps between anchors are searched again with smaller seeds until they are small enough. This is the user's idea, applied to genomes.
- **W:** ClustalW or MUSCLE in the gaps.
- **Rep:** Uniqueness requirement. progressiveMauve scores anchors that are unique in subsets.

**M-GCAT** (Treangen & Messeguer, *BMC Bioinform* 2006, 10.1186/1471-2105-7-433)
- Multi-MUMs, recursively re-found inside the intervals between anchors, giving a hierarchy of anchors.
- I could not identify the "TPA" named in the brief. Not covered.

**Mugsy** (Angiuoli & Salzberg, *Bioinformatics* 2011, 10.1093/bioinformatics/btq665)
- All-vs-all nucmer, then a segment-based alignment graph with T-Coffee-like consistency, then LCBs.
- Consistency comes from merging pairwise anchors, not from multi-MUMs.

**Cactus pinch graphs / Progressive Cactus / Minigraph-Cactus** (Paten et al., *J Comput Biol* 2011, 10.1089/cmb.2010.0252; Armstrong et al., *Nature* 2020, 10.1038/s41586-020-2871-y; Hickey et al., *Nat Biotechnol* 2024, 10.1038/s41587-023-01793-w)
- Pairwise alignments are "pinched" into a graph, then filtered by the cactus structure ("melting" short chains).
- The final step (BAR) realigns bubbles with abPOA.
- MC is the graph being fixed. Its BAR abPOA settings are the `abpoa_mc` strategy in `vntr-eval`.

**SibeliaZ / TwoPaCo** (Minkin & Medvedev, *Nat Commun* 2020, 10.1038/s41467-020-19777-8; TwoPaCo, *Bioinformatics* 2017, 10.1093/bioinformatics/btw609)
- A compacted de Bruijn graph gives anchors; collinear blocks are found by walking bubbles; blocks are aligned with spoa.
- Tandem repeats collapse into cycles. This is useful as a reminder that a DAG POA cannot represent copy-number differences except as deletions. centrolign's `-c` mode recovers cycles afterwards.

### 2.3 Recursive anchoring for DNA (the user's idea, closest form)

**LAGAN / CHAOS / Multi-LAGAN / Shuffle-LAGAN** (Brudno et al., *Genome Res* 2003, 10.1101/gr.926603; CHAOS: Brudno & Morgenstern, *Bioinformatics* 2003, 10.1093/bioinformatics/btg1005)
- **A:** CHAOS local alignments, built from short degenerate seeds chained into short local hits.
- **C:** Sparse DP for the best colinear chain, giving a rough global map.
- **R:** **Recursion inside each gap between consecutive anchors** that is longer than a threshold, rerunning CHAOS with less stringent parameters (shorter seeds). This is exactly the "re-anchor inside partitions" step.
- **W:** Needleman-Wunsch restricted to a limited area (boxes) around the anchor chain.
- **Multi-LAGAN:** progressive along a tree. Profile-vs-profile anchors are the union of pairwise anchors projected into profile coordinates. This is the closest analogue to anchoring against a graph.
- **Rep:** None beyond chaining. Shuffle-LAGAN allows rearrangements.
- What plegma should copy:
  - Anchors are only "good enough to box the DP". The final alignment is optimal within the boxes.
  - Recursion is triggered by gap size, not by a fixed depth.

**DIALIGN** (Morgenstern et al., *PNAS* 1996, 10.1073/pnas.93.22.12098; DIALIGN-TX, *Algorithms Mol Biol* 2008, 10.1186/1748-7188-3-6; anchored DIALIGN, *AMB* 2006, 10.1186/1748-7188-1-6)
- **A:** Gap-free "diagonals" with probability-based weights.
- **C:** Greedy insertion of diagonals, keeping them consistent through transitivity frontiers.
- **R:** None.
- **W:** None. Gaps are what is left between diagonals.
- **Rep:** None.
- Useful idea: anchored DIALIGN accepts **user anchors** as hard constraints. That is the right interface for plegma's MUM backbone.

**MAFFT `--addfragments`, `--anchors`, long-sequence modes** (Katoh & Standley, *Mol Biol Evol* 2013, 10.1093/molbev/mst010; Katoh & Frith, *Bioinformatics* 2012, 10.1093/bioinformatics/bts578)
- `--add` and `--addfragments` align new sequences to an existing MSA. This is the "pairwise to the current alignment" design.
- LAST-based anchoring is used for long genomic sequences (`--addlong`, and the `--large` / `--memsave` modes).
- MAFFT's FFT step itself finds homologous segments by FFT correlation and runs DP between them (an anchor-and-fill design).
- No VNTR model. On these data, L-INS-i is good but slow (`vntr-eval`).

### 2.4 POA with anchors or windows

**abPOA `-S`** (Gao et al., *Bioinformatics* 2021, 10.1093/bioinformatics/btaa963; code `abpoa/src/abpoa_seed.c`, `abpoa_align.c`, read)
- **A:** Minimizers, k = 19 and w = 10 (`ABPOA_MMK`, `ABPOA_MMW`). `collect_anchors1` emits **every pair** of identical minimizers between the target sequence and the query, with no occurrence cap. In a tandem array this is about (copies)² anchors.
- **Target:** `abpoa_build_guide_tree_partition` anchors read `i` against **read `i-1` in the guide-tree order** (`tid = read_id_map[i-1]`), a linear sequence. The anchors are mapped into the graph through that read's `qpos_to_node_id`, swapped into `tpos_to_node_id` (`abpoa_anchor_poa`).
- **C:** minimap2-style DP chaining (`abpoa_dp_chaining`, `get_chain_score`):
  - Score: `min(dq, dt, k) − (log2(|dt−dq|)/2 + 0.01·k·|dt−dq|)`.
  - Hard limits: `max_bw = 100` on the diagonal shift, `max_dis = 100` on the predecessor search, `min_local_chain_score = 100`.
  - Local chains are then chained (`abpoa_dp_chaining_of_local_chains`).
  - Anchors closer than `min_w = 500 + k` are thinned (`ABPOA_MIN_POA_WIN = 500`).
- **R:** None.
- **W:** `abpoa_align_sequence_to_subgraph(beg_id, end_id)` aligns the query slice to the **closed subgraph** between two anchor nodes. `abpoa_subgraph_nodes` extends to the full upstream and downstream closure in topological order. Anchor k-mers are added as forced matches.
- **Rep:** None. In arrays, the wrong-copy diagonals are locally colinear and survive chaining.
- Useful ideas to keep:
  - Mapping anchors to nodes through a path's position-to-node array.
  - Window DP on a closed subgraph between anchor nodes.

**racon / spoa windows** (Vaser et al., *Genome Res* 2017, 10.1101/gr.214270.116)
- Reads are first aligned to a backbone (the draft). The backbone is cut into fixed 500 bp windows, and each read is cut at the backbone positions its CIGAR maps there.
- spoa runs per window with the backbone first, then the windows are stitched.
- The anchoring is the backbone alignment. It is a single linear target, as in abPOA `-S`.

**Flye** (Kolmogorov et al., *Nat Biotechnol* 2019, 10.1038/s41587-019-0072-8)
- The polisher aligns reads to the draft, splits them into "bubbles" between solid regions, and builds a consensus per bubble. My recollection is iterative edit proposals scored by a pair-HMM.
- Repeats are resolved earlier, at the repeat-graph level.
- Relevant only as "solid k-mers bound the windows".

**smoothxg / pggb** (Garrison et al., *Nat Methods* 2024, 10.1038/s41592-024-02430-3; code `smoothxg/src/blocks.cpp`, `smooth.cpp`, not read in depth)
- The graph is cut into blocks by path length and depth. The path segments of each block are re-aligned with abPOA or spoa, then laced back together with padding.
- The blocks are graph-defined, not anchor-defined.
- wfmash and sweepga: mapping, then plane-sweep filtering of pairwise mappings. Not a POA anchoring scheme (inferred from naming and docs).

**minimap2 chaining** (Li, *Bioinformatics* 2018, 10.1093/bioinformatics/bty191; Li 2021, 10.1093/bioinformatics/btab705)
- Chain score: f(i) = max(f(j) + α(j,i) − β(j,i), w_i).
  - α = min(min(Δx, Δy), w_i) rewards new matched bases.
  - β = γ(|Δx − Δy|), with γ(l) = 0.01·w̄·l + 0.5·log2 l. That is linear plus log in the indel imbalance.
- Gaps longer than G are forbidden, and the search is heuristic (max_skip, max_iter).
- Later versions add RMQ-based chaining for long gaps (`mg_lchain_rmq` in minigraph and minichain).
- The log term makes large indels cheap, which is good for copy-number differences. The band and skip limits are bad for them.

**POASTA** (van Dijk et al., optimal gap-affine POA by A*; bioRxiv 10.1101/2024.05.23.595521; code `poasta/`, not read)
- Exact sequence-to-POA alignment. A candidate optimal fill or verification engine for small windows.

### 2.5 Tandem-repeat-specialised aligners and VNTR decomposition

**centrolign** (Eizenga et al., bioRxiv 2026, 10.64898/2026.09.03.749043; code `centrolign/include/centrolign/`, headers read)
- **Target:** Progressive graph-to-graph alignment along a guide tree.
  - `PathMatchFinder` builds an ESA (enhanced suffix array) over the **embedded path sequences** of both graphs (`path_esa.hpp`).
  - `GESAMatchFinder` alternatively indexes all walks of a simplified graph (`gesa.hpp`, with a size-limit factor of 16).
- **A:** "Minimal rare matches": minimal sequences that occur in both graphs but at most `max_count = 50` times in either (`match_finder.hpp`, `esa.hpp`).
  - Counts use a colour-set-size or range-unique-query index over node labels (`use_color_set_size = true`).
  - As far as I can tell from the comments ("original node IDs of path nodes"), this counts **distinct graph positions**, so one location shared by many paths counts once. Inferred.
- **Weights** (`score_function.hpp`): the default `ConcaveLengthScaleInverseCount` is

  w = (len / (c1·c2)^0.5) − (len/1750)³·1750

  scaled by the partial-match fraction. Rare beats long, and very long matches are penalised.
- **C:** `SparseAffine` chaining, with a 3-piece affine gap (`gap_open {1.25, 50, 5000}`, `gap_extend {2.5, 0.1, 0.0015}`), giving a concave gap cost so huge copy-number indels are cheap per base.
  - Reachability comes from a chain decomposition of the DAG (`chain_merge.hpp`).
  - `autocalibrate_gap_penalties` rescales the gaps to the anchor score scale.
  - Matches are split at branch points (`split_matches_at_branchpoints`) to avoid reachability artefacts.
- **R:** `do_fill_in_anchoring`: after the chain, the subgraphs between anchors are extracted (`extract_graphs_between`), the excluded matches are redistributed to them (`divvy_matches`) under a pair budget, and they are re-chained (`fill_in_anchor_chain`). This is LAGAN-style recursion.
- **W:** `stitcher.hpp` fills between anchors (details not read).
- **Rep:** Rare-match weighting, a concave gap, and an optional cyclic graph (`-c`) for large tandem duplications.
- On the vntr-eval hard set, a centrolign star scored below the abPOA medoid star (0.659 vs 0.703), probably because its parameters are tuned for megabase arrays. Inferred.

**UniAligner / TandemAligner** (Bzikadze & Pevzner, *Nat Methods* 2023, 10.1038/s41592-023-01970-4)
- Parameter-free pairwise alignment for centromeres. Anchors are **rare matches** with frequency-penalised scores, chained, and then filled. centrolign credits it as its design basis.

**vamos** (Ren, Gu, Chaisson, *Genome Biol* 2023, 10.1186/s13059-023-03010-y; code `vamos/src/bounded_anno_local.cpp`, skimmed)
- Each VNTR allele is decomposed into a string of motifs from an "efficient motif set" by a wrap-around motif DP.
- Alleles are compared by edit distance over motif strings.
- spoa and abPOA are used for read consensus (`msa.cpp`).
- Gives a **unit alphabet** for unit-level alignment.

**TRGT** (Dolzhenko et al., *Nat Biotechnol* 2024, 10.1038/s41587-023-02057-3; code `trgt/src/hmm/builder.rs`, read)
- Builds a profile HMM of motif blocks: one per motif, with 3n + 1 states, plus a skip block and run start/end states.
- Viterbi labels **motif spans** in each allele. This is the motif tokenisation used for purity and visualisation.
- Directly reusable as a unit segmenter.

**TRF / TRAL / ULTRA / straglr**
- TRF: Benson, *NAR* 1999, 10.1093/nar/27.2.573.
- TRAL: Schaper et al., *Bioinformatics* 2015, 10.1093/bioinformatics/btv306, DOI unverified.
- ULTRA: Olson & Wheeler, ACM-BCB 2018, 10.1145/3233547.3233604.
- straglr: Chiu et al., *Genome Biol* 2021, 10.1186/s13059-021-02447-3.
- These are detection and genotyping tools (period, array extent, copy number). In plegma they serve as **array annotators** that switch on repeat mode for a window.

**Unit-level alignment with duplications**
- Benson, "Sequence alignment with tandem duplication", *J Comput Biol* 1997, 10.1089/cmb.1997.4.351.
- Bérard & Rivals, MS-Align, minisatellite map comparison, *J Comput Biol* 2003, 10.1089/10665270360688048.
- Sammeth & Stoye, *IEEE/ACM TCBB* 2006.
- These align **strings of repeat units** with unit duplication and contraction operations. The vntr-eval unit-aware aligner (`tools/realign_units.py`, `unit_dp.c`) is a working instance on these exact regions.

**VNTR pangenome papers**
- Lu, Chaisson et al., "Profiling VNTR variation across populations using repeat-pangenome graphs", *Nat Commun* 2021, 10.1038/s41467-021-24378-0 (danbing-tk).
- Ren et al. is vamos, above.
- I could not identify the specific "Park et al." VNTR-pangenome paper. Not covered.

### 2.6 Multi-MUM theory, used only where it informs the pairwise design

mumemto (Shivakumar & Langmead; MumemtoM, *Genome Res* 2026; `mumemto/README.md`):
- `-k` sets the minimum number of sequences containing a match ("partial multi-MUMs").
- `-f` caps the occurrences per sequence ("rare multi-MEMs").
- `-F` caps total occurrences.

The pairwise analogues are what plegma needs, and sporos already has them (`~/sporos/README.md`):
- `-mum`
- `-mumreference`
- `-maxmatch --max-occ N` (frequency-capped MEMs)
- f-SMEMs with `--max-occ`

---

## 3. Comparison table

| method | anchors | consistency / chaining | recursion | window aligner | repeat handling | fit for plegma |
|---|---|---|---|---|---|---|
| FMAlign | sampled seeds present in all sequences, FM-index | colinear chain | no | MAFFT / HAlign | frequency cut | low |
| FMAlign2 | multi-MEMs from generalized SA + LCP + DA | avg-position order; per-sequence nearest occurrence; O(m²) DP | no | external MSA per segment | nearest to mean position | idea only |
| HAlign 3/4 | exact matches to centre (suffix tree, then BWT) | pairwise to centre | no | DP / WFA | none | design analogue (star) |
| WMSA | clustering | center-star + profiles | no | DP | none | low |
| MUMmer / nucmer | MUM / mumreference / MEM | LIS (MUMmer 1), diagonal clusters | no | banded extend | drop repeats | fill model |
| Mauve | multi-MUMs, spaced seeds | LCBs + breakpoint penalty | **yes, smaller seeds in gaps** | ClustalW / MUSCLE | uniqueness | recursion model |
| M-GCAT | multi-MUMs | hierarchical | **yes** | external | uniqueness | recursion model |
| Mugsy | pairwise nucmer | segment graph + consistency | no | segment MSA | via nucmer | low |
| Cactus / MC | pairwise alignments | pinch graph + cactus filtering | no | abPOA (BAR) | melting | the thing being fixed |
| SibeliaZ | de Bruijn k-mers | bubbles → blocks | no | spoa | cycles | low |
| LAGAN / M-LAGAN | CHAOS local hits | sparse DP chain | **yes, shorter seeds in gaps > threshold** | limited-area NW | none | **core model** |
| DIALIGN | gap-free diagonals | greedy, transitivity | no | none | none | anchor-constraint interface |
| MAFFT add / anchors | FFT segments / LAST | progressive / add | no | DP | none | baseline |
| abPOA `-S` | minimizer pairs, uncapped | minimap2-like DP, bw 100 | no | POA on closed subgraph | none (fails) | anti-pattern + subgraph trick |
| racon | backbone alignment | linear | no | spoa per 500 bp window | none | low |
| smoothxg | graph blocks | n/a | no | abPOA / spoa per block | none | low |
| minimap2 | minimizers | linear + log gap, banded | no | ksw2 | occurrence cap | gap function |
| **centrolign** | **rare matches ≤ 50 in graph paths (node-distinct)** | **sparse concave-affine chain on DAG** | **yes (fill-in)** | stitcher | rarity weights, concave gap, cycles | **closest overall** |
| UniAligner | rare matches | chain | partial | DP | rarity | anchor weighting |
| vamos / TRGT | motif tokens (DP / HMM) | n/a | n/a | motif-level ED | motif model | unit fallback |
| MS-Align / Benson | repeat units | unit DP with duplications | n/a | unit DP | explicit | unit fallback |
| minichain | minimizers on graph nodes | haplotype-aware chaining, recombination penalty, path cover | no | WFA (minigraph) | n/a | graph chaining |

---

## 4. (a) Anchoring a new sequence against a partial-order graph

### 4.1 What to index

| option | how | uniqueness meaning | pros | cons |
|---|---|---|---|---|
| Previous sequence in order (abPOA `-S`) | sporos index of one earlier sequence | unique in that sequence | trivial; node map = that path | arbitrary; wrong copy number for the query |
| Consensus (heaviest bundle) | index the consensus string; node map = consensus path | unique on the consensus | one string; stable | misses non-consensus alleles; consensus can be a recombinant |
| **Closest path** (k-mer distance, copy-number aware) | index one chosen path, e.g. nearest by weighted Jaccard of 16-mers as in vntr-eval's `guide` | unique on that path | pairwise MUM theory applies directly; good anchors where the query resembles it | optimality only relative to that path |
| Medoid or longest path | as above, fixed per region | as above | matches the vntr-eval medoid star and longest-first evidence; holds the most copies | same |
| **All paths, node-distinct counting** (centrolign `PathMatchFinder`) | concatenate all embedded paths (or distinct alleles only); per text position store (path, offset, node, node_offset) | a match is unique if all its text occurrences map to **one start (node, offset)** | a substring repeated across 456 haplotypes through the same nodes counts once; anchors anywhere in the graph | sporos `--max-occ` counts text occurrences, so a post-filter is needed (below); duplicate branches with identical sequence look like repeats |
| Graph walks (GCSA/GESA) | index all walks, with prefix doubling | unique over the graph's language | includes recombinant paths | blow-up (centrolign caps at 16× graph size); spurious recombinant matches |
| Topological linearisation | concatenate node labels in topological order | none sensible | trivial | substrings across branch boundaries are not real sequence; reject |

### 4.2 Recommendation and how to compute node-distinct uniqueness

Start with **one reference path** (longest or medoid), with uniqueness defined on that path. Then move to **all distinct alleles in the graph** with node-distinct counting.

The node-distinct variant with sporos:
1. **Index.** Build the text from distinct path sequences (deduplicated, as pgrealign and vntr-eval do; hard-set regions have a median of 220 distinct alleles).
2. **Query.** Run all MEMs ≥ L with `--max-occ` set to `N_paths × c_max`, or uncapped, and filter afterwards.
3. **Map.** Convert each occurrence to `(start_node, start_offset)` with a `pos → node` array built while concatenating.
4. **Count.** Collect the distinct (node, offset) values per match. Since N ≤ 456, enumerate the SA interval and exit early once the distinct count exceeds c_max.
5. **Define:**
   - `count_G` = number of distinct graph start positions.
   - `count_Q` = occurrences in the query, from sporos `-mum` semantics or a self-index of Q.
   - Graph-MUM: count_G = 1 and count_Q = 1.
   - Rare: both ≤ c_max.
6. **Diverging walks.** A match whose text occurrences start at the same node but diverge into different walks is still **one graph location with several walks**. Keep the walks (centrolign's `match_set_t` stores `walks1` and `walks2`) and split the match at the branch point if the walks differ in length downstream. This follows centrolign's `split_matches_at_branchpoints`.

### 4.3 Mapping matches back to nodes

- Each anchor becomes a list of (query pos, node id) pairs, from the per-path `pos → node` array (abPOA's `tpos_to_node_id`).
- An anchor can span several nodes. Chaining needs **reachability between anchor end-node and next start-node**, and a graph distance for the gap cost.
- Path-specific distance: `dist2begin` along path k, as in minichain `MPC_index`. Reachability between paths: last2reach and a chain-decomposition index (minichain; centrolign `chain_merge.hpp`).

---

## 5. (b) Prior art for recursive pairwise anchoring

- **LAGAN / CHAOS:** anchor, chain, then rerun the seeding with shorter or less strict seeds **inside each gap longer than a threshold**. The final NW is limited to boxes around the anchors. Multi-LAGAN applies the same steps profile-vs-profile.
- **Mauve / M-GCAT:** recursive anchoring with smaller seeds inside the gaps between anchors.
- **centrolign fill-in:** after the global chain, extract the subgraphs between consecutive anchors and re-chain the previously excluded rare matches inside each, under a match-pair budget.
- **MUMmer / nucmer:** single level. Clusters are found, then extended by banded DP between anchors.
- **minimap2 / abPOA `-S`:** single level. abPOA partitions into POA windows of at least 500 bp between anchors and aligns each window to the closed subgraph between anchor nodes.
- **racon / Flye:** fixed windows on a backbone (racon) or solid-region-bounded bubbles (Flye). No recursion.

What recursion buys:
- **Uniqueness is relative.** A k-mer that occurs in every copy of the whole array can be unique inside a window that holds one copy. Re-counting inside the window is what turns SMEMs into MUMs.
- **Thresholds can shrink.** L_ℓ can decrease with depth (for example 51 → 31 → 21 → 15) while the false-match rate stays controlled, because the search space shrinks.

What recursion costs: every level freezes the level above, as greedy refinement does. So the top level must be conservative.

---

## 6. Optimal anchor chaining: theory, code, and what it guarantees

### 6.1 Colinear chaining with overlaps and gap costs (Jain, Gibney, Thankachan, RECOMB 2022 / JCB 2022, 10.1089/cmb.2022.0266; code `ChainX/src/include/algo.hpp`, read)

- **Anchors:** exact matches ([a..b] in S1, [c..d] in S2), plus sentinel anchors at both ends.
- **Strict precedence:** I′ ≺ I if a′ ≤ a, b′ ≤ b, c′ ≤ c, d′ ≤ d, with at least one inequality strict. Overlaps are allowed.
- **Weak precedence:** start points only.
- **connect(I′, I) = g + o**, where:
  - g = max(g1, g2), with g1 = max(0, a − b′ − 1) and g2 = max(0, c − d′ − 1);
  - o = |o1 − o2|, with o1 = max(0, b′ − a + 1) and o2 = max(0, d′ − c + 1). This is the overlap imbalance, the diagonal shift inside the overlap.
- **Theorem 2:**
  - For a fixed anchor set A, anchored edit distance = optimal chain cost (strict) = optimal chain cost (weak).
  - Anchored edit distance is ED in which a match not supported by an anchor costs 1.
- **Complexity:**
  - Exact: Õ(n), specifically O(n log⁴ n) strict, O(n log³ n) weak, and O(n log² n) for fixed-length anchors, with up to 4-D RMQs over six cases.
  - ChainX ships the practical **O(n·OPT + n log n) average-case** heuristic. In `compute_global`, `bound_redit` starts at 100. Predecessors are searched only within `j_a − i_a ≤ bound`, and the bound is multiplied by 4 until the cost is at most the bound.
  - `DP_global` is the O(|S1||S2|) reference: anchored ED with `matchAllowed` cells.
- **Empirical:**
  - Correlation with ED is above 0.9.
  - **MUM anchors saturate as L drops**, because short matches stop being unique.
  - **MEM anchors keep improving** as L drops, but their count and run time explode (Table 3: MEMs ≥ 7 take 2,652-9,221 s vs MUMs 2-5 s).

### 6.2 Chaining with overlaps revisited (Mäkinen & Sahlin, CPM 2020, 10.4230/LIPIcs.CPM.2020.25; arXiv 2001.06864)

- They maximise a **symmetric ordered coverage**: each anchor adds the minimum, over the two strings, of the newly covered length. The resulting chain induces an alignment with exactly `coverage` matched characters.
- This is an LCS restricted to the anchored matches.
- Complexity: O(n log² n) under strict precedence with the equal-match-length property, which MEMs have. O(n log n) under weak precedence or for k-mers.
- Shibuya & Kurochkin (WABI 2003) is the origin of the idea.

### 6.3 Guarantees plegma can rely on (inferred, from 6.1 and 6.2)

Let A be the anchor set and C* the optimal chain under ChainX's cost.
1. **Gap filling helps.** cost(final alignment after exact gap filling) ≤ cost(C*) = anchored ED(A). Exact DP inside each gap is at least as good as the greedy filling in the proof of Lemma 1. Plegma's alignment is therefore never worse than its chain cost, and it is optimal among alignments through C*'s anchors.
2. **All MEM occurrences ≥ L give an additive bound.**
   - Take an optimal alignment and any maximal run of exact matches in it of length ≥ L. It lies on one diagonal and extends to a MEM occurrence pair on that diagonal, so it is supported.
   - Hence ED ≤ anchored ED(A) ≤ ED + m_<L, where m_<L is the number of matched bases in runs shorter than L. With ED computed under the same unit cost, the bound is tight in the obvious sense.
   - This is the guarantee to aim for.
3. **MUMs, SMEMs and capped MEMs lose that bound.** No additive bound holds for MUM-only, SMEM-only or capped MEM sets.
   - SMEMs drop a MEM whose query interval is contained in a longer MEM's interval. In a tandem array, the longest match for a query copy is often **off-diagonal**, spanning a run of ref copies at a shifted phase, so the on-diagonal match disappears.
   - Caps drop all matches in arrays with more copies than the cap.
   - The loss is exactly the anchor support those dropped matches would have given.
4. **Pairwise chaining against one path is not optimal against the graph.** The bound holds relative to the indexed text. For the graph, the analogous object is chaining on a DAG (6.4).

### 6.4 Sequence-to-DAG chaining

- **Kuosmanen et al., RECOMB 2018** ("Using minimum path cover to boost DP on DAGs: co-linear chaining extended") and **Mäkinen et al., "Sparse dynamic programming on DAGs with small width"**, *ACM TALG* 2019, 10.1145/3301312.
  - Chaining over a DAG in O(k·n log n) with a minimum path cover of width k.
  - Reachability comes from a path cover with forward propagation (last2reach).
- **Rizzo, Cáceres, Mäkinen**:
  - "Chaining of maximal exact matches in graphs", SPIRE 2023, 10.1007/978-3-031-43980-3_29.
  - "Exploiting uniqueness: seed-chain-extend alignment on elastic founder graphs", *Bioinformatics* 2023, 10.1093/bioinformatics/btad225.
  - Together they define MEMs between a string and a graph and chain them. The founder-graph paper uses **graph-level uniqueness of node labels** for seeding. That is the formal version of section 4.2.
- **GraphChainer** (Ma et al., *Bioinformatics* 2023, 10.1093/bioinformatics/btad460): path-cover chaining for long reads on variation graphs.
- **Chandra & Jain**:
  - "Sequence to graph alignment using gap-sensitive co-linear chaining", RECOMB 2023, 10.1101/2022.08.29.505691.
  - "Haplotype-aware sequence-to-graph alignment", RECOMB 2024, 10.1101/2023.11.15.566493.
  - Code `minichain/graphUtils.cpp`, read:
    - `MPC()` uses the **GFA walks themselves as the path cover** when they exist (`is_hprc_gfa`), and a greedy minimum path cover otherwise.
    - `MPC_index()` precomputes `index`, `dist2begin`, `last2reach` and `Distance` per path.
    - `Chaining()` keeps a DP value per (anchor, path), `C[i][j]`, with one RMQ structure (an AVL tree) per path:
      - on the same path, the gap is linear in path distance and query distance;
      - switching path costs the **recombination penalty R**, through `C_m − recomb` and `C_p` for anchors reached via last2reach.
    - The cost is O(K·n log n) for K paths.
- For plegma the embedded paths are exactly the previously added sequences, at most 456 of them. **minichain's haplotype-aware formulation is directly usable.** Recombination should be allowed, but at a cost, since a new haplotype is often a mosaic of earlier ones.

### 6.5 minimap2's score (for reference)

- α = min(Δx, Δy, w) and β = 0.01·w̄·|Δx − Δy| + 0.5·log2|Δx − Δy|.
- The log term makes long indels cheap, which is good for copy-number differences.
- The heuristic limits (G, band, max_skip) are what break it in arrays. abPOA inherits a 100 bp band.

### 6.6 How multi-occurrence matches should enter the chain

- **MUMs** (count_Q = count_G = 1): single candidate, full weight. Treat them as near-hard constraints at the top level only when they are also **robust**. Robust means removing the anchor drops the chain score by more than a margin (a MAPQ-like confidence), or the second-best occurrence at a shorter length is far shorter.
- **Rare MEMs** (counts ≤ c_max): enter **every occurrence** as a separate candidate anchor, with weight discounted for ambiguity:
  - w = len·u, with u = (count_Q·count_G)^(−1/2) (centrolign's `pair_count_power = 0.5`);
  - or w = len − λ·log(count_Q·count_G) (inferred alternative).
  - The chain chooses among occurrences by colinearity. Ambiguity resolves itself when neighbouring anchors pin the diagonal.
- **SMEMs:** use them to *find* candidate intervals cheaply. For chaining, either enumerate all MEM occurrences in the SMEM's query interval (re-seeding, as BWA-MEM does), or treat each SMEM occurrence as a candidate with the same rarity discount. Never take "the longest SMEM" as a hard partition point inside an array.
- **Frequent matches** (> c_max): exclude them at this level. They come back in recursion, where they are re-counted in the window.

---

## 7. Proposed chaining formulation for plegma

### 7.1 Objects

- Target: the current graph G with embedded paths P_1..P_m (the sequences added so far). Query: Q.
- Anchor types at level ℓ, with minimum length L_ℓ (for example 41, 31, 21, 15):
  1. **G-MUM:** count_G = 1 (node-distinct) and count_Q = 1.
  2. **Rare MEM occurrence:** 1 < count_G·count_Q, with both counts ≤ c_max (for example 8-50). sporos `-maxmatch --max-occ` provides these, followed by a node-distinct post-filter.
  3. **Sentinels:** the source/sink nodes paired with Q's ends. The flanking anchor nodes of a region are guaranteed shared, since pgrealign regions sit between nodes every walk passes.
- Each anchor carries: (q_start, q_end), its start and end (node, offset), the set of paths through it, count_G, count_Q, and length.

### 7.2 Weight

w(a) = len(a)·(count_Q·count_G)^(−1/2) − β·max(0, len − L_long)³

Notes:
- The concave tail follows centrolign and only matters for very long matches, which are usually flanks. It can be off at first.
- MUMs get u = 1.
- Score in "matched bases", so that the gap cost is on the same scale. ChainX's unit cost is the cleanest base.

### 7.3 Precedence and overlap

- Strict precedence: start and end both non-decreasing in Q and along a path, with reachability in G.
- Allow overlaps. Charge o = |o_q − o_g|, the diagonal shift inside the overlap, as in ChainX and Jain et al.
- Overlaps occur constantly between MEMs of a periodic array (consecutive copies), so forbidding them would force artificial gaps.

### 7.4 Gap cost (indel imbalance d = g_q − g_g; substitution part s = min(g_q, g_g))

connect = σ·s + o + γ(d) + R·[path switch]

- **σ** is the expected cost per gap base that is not an indel. ChainX uses 1. Use 1 for the guarantee and tune afterwards.
- **γ(d)** is **concave and period-aware**:

  γ(d) = min_k(open_k + ext_k·|d|) + ρ·r(d)

  - The first term is a 2-3 piece affine, as in centrolign.
  - r(d) = min(|d| mod p, p − |d| mod p) is the distance from |d| to the nearest multiple of the local period p. p comes from the region's catalogue motif (pgrealign catalogs / `hardset/chr20.tsv` have period and motif) or from a TRF-style self-alignment.
  - A shift by whole copies is cheap. A shift by a partial copy pays ρ per base.
- **Implementation:** a concave piecewise-linear γ is the minimum of affine pieces. "Maximise f − min_k(piece_k)" equals "max over k of (f − piece_k)", so each piece is one RMQ structure, which is centrolign's approach. The ρ·r(d) term is periodic, not affine. It can be handled by bucketing anchors by diagonal residue mod p (p RMQ structures, or p ≤ ~100 buckets), or applied only as a re-scoring term in the heuristic ChainX loop. Inferred design.
- **R** is the recombination penalty for changing path, as in minichain. A first version can be a single-path chain (R = ∞) against the chosen reference path.

### 7.5 Algorithm per level

1. Sort anchors by query start.
2. DP as in minichain, with one RMQ structure per (path, gap piece). Or use the ChainX O(n·OPT) banded heuristic first, which is simpler and exact when the guess-and-quadruple loop terminates.
3. Tie-breaking (essential for VNTRs): among equal-score predecessors, prefer the one that puts the indel imbalance **leftmost**, meaning the smallest graph coordinate of the gap. Implement it as a lexicographic score (score, −gap_position) or a tiny ε·position penalty. Every haplotype's copy-number gap then lands at the same equivalent position, and POA merges them.

### 7.6 Recursion (LAGAN-style)

1. For each pair of consecutive chained anchors (a_i, a_{i+1}), take the query interval I_q and the subgraph G_i: the nodes reachable from a_i.end that can reach a_{i+1}.start (abPOA's closure trick).
2. If max(|I_q|, span(G_i)) ≤ W_DP (the window budget for karan banded DP), or |d| plus a margin fits the band, stop.
3. Otherwise, recompute matches **restricted to the window**, with L_{ℓ+1} < L_ℓ and counts recomputed inside the window. Chain them globally between the two fixed endpoints, then recurse.
4. Stop when there are no anchors, or at a minimum length L_min (about 13-15 for 30 kb windows, where random matches become likely).
5. If a window is still too big and is inside an annotated array, switch to the **repeat fallback**:
   - (i) full banded DP with band ≥ |d| + slack, which is the accuracy baseline; or
   - (ii) unit-level DP: tokenise the query and graph paths with TRGT-style motif HMMs, vamos-style DP or vntr-eval's `unit_dp`; align the token strings with unit duplication operations (MS-Align / Benson); then fill each unit pair with karan.

### 7.7 Complexity

- Anchors per level: about n_MUM + Σ_rare(count_Q·count_G) ≤ n·c_max².
- Exact chaining: O(n log² n) to O(n log⁴ n) per level for one path; O(K·n log n) for haplotype-aware chaining with K paths.
- The ChainX heuristic is O(n·OPT + n log n) average.
- A recursion over disjoint windows costs at most O(depth) times one level.
- Perfect arrays (period p, array length A) produce about (A/p)² MEM occurrence pairs; for example A = 10 kb and p = 30 give about 1.1·10⁵. The cap c_max removes them at the top level, which is intended.

### 7.8 What it guarantees

1. Each level's chain is optimal for its anchor set under the stated cost (exact algorithm), or optimal when the ChainX loop terminates (heuristic).
2. With σ = 1, ρ = 0, linear γ and uncapped MEMs, the chain cost equals the anchored ED. The final alignment then satisfies ED ≤ cost ≤ ED + m_<L for the indexed path (sections 6.1 and 6.3).
3. Caps, SMEMs, period terms and recursion give up the additive bound. What remains:
   - the final alignment is optimal given the frozen anchors;
   - it is never worse than the chain cost.
4. There is no guarantee relative to the multiple-alignment objective. Progressive POA is order-dependent (section 8).

---

## 8. (c) Order of addition and choice of anchor target

1. **Copy-number placement is an equivalence problem.** In a perfect array, every placement of a k-copy deletion has the same edit cost. The "right" placement is a convention, unless interior unit variants make one placement strictly better.
   - Pairwise optimal chains are therefore non-unique.
   - Different haplotypes will choose different equivalent positions unless the tie-breaking is deterministic **in graph coordinates**.
   - Inconsistent choices create parallel deletion edges and unaligned copies. That raises induced pairwise cost and graph size. vntr-eval's star_gapnorm effect (+0.014 pair F1) is the measured symptom.
2. **Longest-first makes every later copy-number difference a deletion against existing nodes.** Inferred explanation of the vntr-eval result. It holds when the graph already contains the maximum copy number. A later, longer sequence would need insertions, which POA represents as new unaligned branches that later sequences may or may not join. Ordering by **estimated copy number**, descending, is the natural refinement.
3. **What to anchor against:**
   - Against the **previous** sequence (abPOA `-S`), the anchor diagonal depends on the predecessor's copy number. Anchors in arrays then jump between equivalent phases from one sequence to the next.
   - Against a **fixed reference path** (longest or medoid), all sequences use one coordinate system. Equivalent gaps land consistently given deterministic tie-breaking.
   - Against **all paths with node-distinct uniqueness**, anchors can land on any allele. The chain may recombine across paths, which is right for mosaics but adds freedom in arrays. Keep R high inside arrays.
4. **Guide-tree order (nearest first)** keeps each new sequence close to some existing path. That favours anchoring against the closest path, but it was worse than longest-first in vntr-eval (default/guide 1.073 vs 1.067 median; random-order spread larger).
5. **Weights and bundles.** POA's heaviest-bundle consensus reflects the order. Anchoring against the consensus inherits whatever was added first. Avoid it for anchoring in arrays.

---

## 9. (d) Critical assessment of the pairwise MUM → SMEM recursive design

### Where it works
- **Unique flanks.** Region boundaries come from pgrealign anchors, and the unique sequence around arrays yields long MUMs. They fix the global frame and bound the DP to the array plus margins.
- **Arrays with interspersed unique variants.** Unit variants that are rare across copies give G-MUMs or rare MEMs inside the array. With rarity weights and colinearity, these pin the phase and fix copy placement.
- **Compound arrays** (two motifs, or STR+VNTR). The junctions between sub-arrays are usually unique and become natural partition points.
- **Long regions (20-45 kb)** where abPOA ran out of memory or timed out. Flank and junction anchors split the problem into windows that fit karan.

### Where it fails
1. **Perfect tandem arrays have no unique interior anchors.** Every substring shorter than the array occurs in every copy: count_G ≥ copies.
   - MUMs vanish.
   - The "largest SMEM" is a match across many copies at an arbitrary phase. As an anchor it fixes the copy-number gap at the array end it touches, which is arbitrary.
   - The recursion finds nothing new, because inside a window the copies are still copies.
2. **SMEMs in arrays are ambiguous and off-diagonal.** The supermaximal match for a query copy often spans a shifted run of reference copies. It is long, so it outweighs the true-diagonal shorter MEMs, which SMEM enumeration drops. This is the same mechanism that sinks abPOA `-S`.
3. **Copy-number differences shift anchors.** The diagonal between consecutive true anchors differs by k·p. Chainers with bands (minimap2, abPOA: 100 bp) or linear gap costs mis-chain or break.
4. **Top-level freezing.** One wrong top-level anchor inside an array cannot be undone by later levels.
5. **Uniqueness against the wrong target.** Text-occurrence counting over all paths makes every shared substring look repeated, 456 times. Counting against a single path makes uniqueness depend on that path's copy number.
6. **Progressive inconsistency.** Even with perfect pairwise chains, equivalent placements differ between haplotypes unless the tie-breaking is global (section 8).
7. **Theory gap.** "Optimal backbone" from MUMs and SMEMs has no approximation guarantee to ED (section 6.3). The guarantee needs all MEM occurrences.

### Concrete fixes from the literature
- **Node-distinct uniqueness over all paths**: centrolign PathESA with colour-set-size counting; founder-graph uniqueness (Rizzo et al.).
- **Rare MEMs with occurrence caps, all occurrences as candidates, rarity-discounted weights**: centrolign, UniAligner, sporos `--max-occ`, mumemto `-f` analogue.
- **Overlap-aware chaining with gap costs** (ChainX), on the DAG with path switching (minichain).
- **Concave or period-aware gap cost.** Large copy-number indels must be cheap, as in centrolign's 3-piece affine and minimap2's log term, with no band limit.
- **Recompute uniqueness inside windows; shrink L with depth**: LAGAN, Mauve, centrolign fill-in.
- **Robustness filter on top-level anchors** (score-drop or uniqueness margin). Allow only flank and junction anchors at level 0.
- **Deterministic left-normalised tie-breaking** in graph coordinates (vntr-eval gap_norm).
- **Order by length or copy number**, descending (vntr-eval).
- **Inside anchorless arrays, do not anchor.** Use full banded DP, which is the measured best, or unit-level alignment (vntr-eval unit-aware 1.048; MS-Align; vamos/TRGT tokenisation).

---

## 10. Experiment plan, in order

All on `vntr-eval` packages:
- the 11-region pilot, then all 149 regions (hap32, then the full panel);
- then the 45-region hard set;
- then the held-out chr6.

Metrics, using `tools/evaluate.py` and `hardset_strategies.py`:
- cost/opt (all pairs);
- pair F1 (left-normalised optima);
- graph nodes and sequence;
- time and peak memory.

Baselines:
- MC;
- abPOA default/longest;
- spoa longest;
- abPOA `-S`;
- the pgrealign medoid star and `select`;
- unit_aware;
- centrolign star.

1. **E1. Anchor census (no alignment).**
   - **Setup:** for each region, use the reference path R (longest, medoid) and each query. Compute with sporos: MUMs, `-mumreference`, all MEMs, frequency-capped MEMs (c ∈ {2, 8, 50}) and SMEMs, at L ∈ {15, 21, 31, 51}.
   - **Labels:** label each anchor as inside or outside the catalogue array, and whether it lies on **some** optimal pairwise alignment. Use forward/backward DP with the edlib or karan exact score to mark co-optimal cells.
   - **Report:** counts, query coverage, the fraction of co-optimal anchors per type and position, and the off-diagonal fraction for SMEMs.
   - **Decides:** which anchor types are safe at level 0 and the c_max and L schedule.
2. **E2. Pairwise chain quality vs optimum.**
   - **Chainers:** ChainX cost (port the `compute_global` loop: about 50 lines, O(n·OPT)), minimap2-like, and centrolign-like weights.
   - **Anchor sets:** {MUM}, {MUM + SMEM}, {MUM + rare MEM}, {all MEM}.
   - **Measure:** anchored ED / ED, and the fraction of chained anchors that are co-optimal.
   - **Checks:** verify Theorem 2 numerically on small regions against `DP_global`, and the inferred additive bound in section 6.3.
3. **E3. Gap cost.**
   - **Variants:** linear, 3-piece concave (centrolign parameters rescaled), and period-aware ρ·r(d) with the catalogue period.
   - **Metric:** E2's metrics, at hotspot VNTRs with large copy-number differences.
4. **E4. Recursion schedule.**
   - **Variants:** depth 1 vs LAGAN-style, with L schedules {51 → 31 → 21 → 15} and in-window re-counting. Budget sweep W_DP ∈ {1, 2, 5, 10} kb.
   - **Measure:** pairwise cost ratio vs full DP, time and memory.
   - **Decides:** the window budget for karan.
5. **E5. Pairwise end-to-end with karan fill.**
   - **Setup:** chain + banded fill vs karan or edlib full DP.
   - **Target:** cost ratio ≤ 1.001 outside arrays and ≤ 1.01 overall, at lower time and memory.
6. **E6. Anchor target in the progressive POA** (order fixed to longest-first).
   - **Targets:** previous sequence (abPOA `-S` analogue), consensus, fixed longest path, fixed medoid path, nearest path, and all paths with node-distinct counting (with and without the recombination penalty).
   - **Metrics:** cost/opt and pair F1 on the 11 pilot regions, then 149.
7. **E7. Order.**
   - **Variants:** longest-first, copy-number-descending, guide (Prim from the medoid), given, and 3 random seeds. Measure the spread.
   - **Expectation:** anchoring against a fixed path should shrink the order spread, to below 0.05 random-order spread. Inferred hypothesis.
8. **E8. Tie-breaking and normalisation.**
   - **Variants:** left vs right vs none, and unit-variant-aware.
   - **Metrics:** pair F1, graph size and cost/opt. star_gapnorm is the reference.
9. **E9. Array fallback.**
   - **Setup:** in windows with no safe anchors, compare full banded DP, unit-level DP and forced SMEM anchors (negative control).
   - **Expectation:** forced SMEM anchors reproduce the abPOA `-S` failure.
10. **E10. Scale.**
    - **Setup:** the full panel (up to 456 haplotypes, or distinct alleles with weights), regions > 20 kb (L001909 and the four > 27 kb hard-set regions).
    - **Target:** fit in 12 GB where abPOA ran out of memory.
11. **E11. Hard set and held-out chr6.**
    - **Acceptance:** beat abPOA default/longest on hotspot cost/opt (1.067 median) and the medoid star on hard-set pair F1 (0.703), with no control regressions above 0.005.

---

## 11. References (code files read)

- `FMAlign2/src/mem_finder.cpp`, `FMAlign2/src/sequence_split_align.cpp`, `FMAlign2/include/mem_finder.h`
- `FMAlign/common_seed.cpp`, `FMAlign/README.md`
- `abpoa/src/abpoa_seed.c` (guide tree, `collect_anchors1`, `abpoa_dp_chaining`, `get_chain_score`), `abpoa/src/abpoa_align.c` (`abpoa_anchor_poa`), `abpoa/src/abpoa_graph.c` (`abpoa_subgraph_nodes`), `abpoa/src/abpoa_align.h`
- `ChainX/src/include/algo.hpp`, `ChainX/src/chainx.cpp`
- `minichain/graphUtils.cpp` (`MPC`, `MPC_index`, `Chaining`), `minichain/README.md`
- `centrolign/include/centrolign/{match_finder,anchorer,score_function,chain_merge,esa,path_esa,gesa,induced_match_finder}.hpp`, `centrolign/README.md`
- `trgt/src/hmm/builder.rs`
- `vamos/src/bounded_anno_local.cpp`, `vamos/src/msa.cpp`
- `mumemto/README.md`
- `vntr-eval/tools/POA.md`, `vntr-eval/tools/UNIT_ALIGN.md`, `vntr-eval/hardset/README.md`, `vntr-eval/README.md`
- `pgrealign/README.md`
- `~/sporos/README.md` (anchor types available)

All repos were cloned at depth 1 under `/private/tmp/claude-501/-Users-kokyriakidis/8500d74c-03af-4a3f-8699-f1899147a1f9/scratchpad/src/`. vntr-eval, pgrealign, smoothxg, spoa, poasta, theseus and cactus were already present from another task. They were read only; nothing was built or run.

Papers are cited inline with their DOIs. The JCB paper on chaining with overlaps (Jain et al. 2022) and the CPM paper (Mäkinen & Sahlin 2020) were read in PDF. The following were not checked against the paper:
- the exact minimap2 constants in later versions;
- the TRAL and WMSA DOIs;
- the details of Flye's polisher;
- UniAligner's exact anchor score.
