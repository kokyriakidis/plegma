# plegma: design

What plegma must do, what its rivals do well and badly, the algorithm, and the experiments that decide its open
choices. The evidence is in `docs/research/`: one report per rival (abPOA, SPOA, POASTA, Theseus, smoothxg with
Minigraph-Cactus and pgrealign), a survey of anchored and repeat-aware alignment (`anchored-msa.md`), and a briefing
on the test bed (`vntr-eval.md`). Numbers quoted here come from those reports, which cite their sources; anything
marked *inferred* is reasoning, not measurement.

## 1. The job

Re-align the VNTR regions that the HPRC Minigraph-Cactus (MC) graph aligns wrongly, and rebuild each region's graph
with every path spelling its haplotype exactly. A region holds 34 (hap32) to about 456 (full panel) haplotypes. At
the hotspot VNTRs the median region spans 5.7 kb on CHM13 (up to 55 kb, full-panel alleles up to 115 kb), the
median repeat unit is 32 bp, and the median array has 153 copies. Every haplotype starts and ends on the same
anchor bases, so the alignment is global.

plegma is to be the fastest, the leanest and the most accurate POA aligner. Accuracy is measured with vntr-eval:

| measure | MC | best so far | target |
| :-- | --: | --: | --: |
| cost/opt at hotspot VNTRs, median (unit cost) | 1.566 | 1.048 (unit-aware), 1.067 (spoa, longest first) | <= 1.1 in every region, below 1.048 |
| non-repeat SV controls | 1.000 | 1.000 | exactly 1.000 |
| matched controls | 1.102 | 1.015 | no region worse than MC |
| `kmer_frac_extra` at hotspots | 0.60 | 0.40-0.52 | toward the controls' 0.05 |
| hard-set pair F1 (pgrealign qc), median | | 0.703 (medoid star) | above it |

cost/opt alone does not predict genotyping (vntr-eval Stage 4). What does:
- pairwise consistency with each pair's own optimal alignment, especially precision;
- each copy-number indel kept as one block and stacked in the same columns across every haplotype that carries it;
- no junctions the haplotypes do not need.

plegma is judged on all of these.

Speed and memory to beat: abPOA takes 0.1 s per hap32 region (median; 1.2 s at hotspots) but reserves about 20 B
per DP cell even when banded. It refused 11 full-panel regions for memory. SPOA needs 80 B × L². pgrealign's star
took 2,137 s wall for chr20's 14,511 regions. Removing the quadratic term is where plegma can win outright.

## 2. The rivals

| | abPOA | SPOA | POASTA | Theseus | smoothxg | pgrealign |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| kernel | SIMD DP, adaptive band | SIMD DP, full matrix | exact A*, gap-affine | exact graph WFA | SPOA/abPOA blocks | abPOA pairwise star |
| exact | no (band may miss) | yes | yes | yes, but graph end free | no | per pair |
| memory | full matrices even banded | full matrices | grows with search | every stored cell | as its POA | as abPOA |
| seeds | minimizers vs previous sequence | none | none (chain heuristic opt-in) | none | none | none |
| order | input or greedy similarity | input | input | input | longest first | medoid star |
| repeats | seeds and band fail in arrays | order-dependent ties | heuristic blind in arrays | free end mis-aligns | folds copies into cycles | cross-slot scatter |

**abPOA** (`research/abpoa.md`).
- Strengths:
  - a fast SIMD kernel: rows are nodes, lanes run along the query, and a prefix scan handles insertions, in
    16-bit lanes or 32-bit when needed;
  - dual-affine gaps by default;
  - alignment to the subgraph between two anchor nodes, then one fusion per read.
- Weaknesses:
  - The `-S` seeding is the failure plegma exists to fix. Every pair of equal minimizers becomes an anchor, with no
    frequency cap.
    - Anchors come from the previous sequence only, never the graph.
    - The chain rejects diagonal shifts above 100 bp, and its second stage has no gap cost.
    - Anchors are forced in as exact matches.
    - On vntr-eval, `-S` scores 1.894 against plain abPOA's 1.067.
  - The band follows the heaviest path, so it can miss several opposite indels without telling anyone. Issues #9
    and #10 report band failures, #14 memory.
  - It always allocates full matrices.
  - Its consensus is greedy.
  - Issue #37, open: Paten and Hickey asked for user-supplied MUM anchors in satellite arrays. That is plegma's
    niche.
- Take:
  - the kernel layout;
  - an exact band per window, from min/max path lengths between its anchors;
  - band-edge detection with a retry.
- Avoid:
  - forcing non-unique anchors;
  - chaining without gap costs;
  - matrices for the whole graph.

**SPOA** (`research/spoa.md`).
- Strengths:
  - simple, exact local, global and overlap modes;
  - linear, affine and convex gaps;
  - columns kept as aligned-node groups, at most one node per base;
  - `Subgraph` plus an ID remap for windows.
- Weaknesses:
  - full score matrices, about 7-36 GB at 30 kb against a 60k-node graph (*inferred*);
  - no banding or seeding;
  - insertions never merged;
  - ties broken by the oldest edge, so gap placement depends on input order;
  - no stored paths.
- Take: the column model, and score width chosen per window from a bound.
- Avoid: full matrices, and order-dependent tie-breaking.

**POASTA** (`research/poasta.md`, v1.0, released 2026-10-18).
- Strengths:
  - exact A* under linear, affine or two-piece costs;
  - greedy match extension, the graph version of WFA;
  - a bucket queue, and a score-only backtrace;
  - superbubble pruning;
  - an oracle plus randomised differential tests;
  - an unmerged `banded-poa` branch with a certified band-doubling Gotoh.
- Weaknesses:
  - its heuristic charges only path-length mismatch, which is zero across copy-number-variable arrays, so the
    search there degrades to Dijkstra;
  - no caps of any kind;
  - global indexes rebuilt per sequence;
  - input order only;
  - ties decide homology.
- Take:
  - its exact search inside anchored windows, with the heuristic measured to the window's exit anchor;
  - the band-doubling certificate;
  - its test harness.

**Theseus** (`research/theseus.md`).
- Strengths:
  - exact graph WFA whose cost grows with divergence, not graph size;
  - multi-base nodes;
  - cycles allowed;
  - diagonal invalidation as dominance pruning;
  - has a progressive MSA mode;
  - beats POASTA 1.6-17.6× on viral and microbial MSAs.
- Weaknesses:
  - the graph end is always free, so in `X R^10 Y` against `X R^5 Y` it mis-aligns the short right flank rather than
    pay for the copy-number gap;
  - single-affine gaps only;
  - no SIMD, no threads, no score cap;
  - stores every cell;
  - several bugs (argument order, GAF ins/del swapped).
- Take: a graph-WFA kernel for windows whose score is small, with a hard score cap and a fallback.

**smoothxg, MC and pgrealign** (`research/smoothxg-cactus.md`).
- smoothxg:
  - It folds VNTR copies into cycles by design: several pieces of one haplotype go into one POA.
  - It cuts haplotypes at their own offsets, not homologous ones.
  - Block boundaries come from a stochastic sort.
- MC's BAR:
  - It aligns only 10 kb from each end of a gap between anchors.
  - Its 10 kb windows are per-row offsets, not homologous.
  - Minigraph keeps expanded copies as insertions Cactus cannot merge.
- pgrealign, now `cactus-pangenome --pgrealign`, is plegma's real rival on this job: a medoid star with abPOA, left
  normalisation, a homology gate and exact-path splicing.
  - Its weakness is cross-slot scatter: extra copies inserted in different centre slots never merge. That accounts
    for 89.6% of its 3.3 Mb excess on chr20.
- Take:
  - pgrealign's region interface and status codes;
  - deduplication with weights;
  - graphs built only from path steps;
  - full path verification after every splice.

**The common faults** plegma must avoid:
- order dependence;
- tie-breaking that is not in graph coordinates;
- insertion branches that never merge;
- quadratic memory;
- anchors that are trusted without paying for the diagonal shift they imply.

## 3. What the test bed already shows

From the pairwise MUM census (`research/vntr-eval.md` 5.4, scripts in `experiments/e0-mum-census/`), hap32 at
hotspot VNTRs:
- Pairwise MUMs inside arrays are common: 99% of CHM13 pairs have one of at least 20 bp, median 18 per pair.
- Most are false homology: only 14-21% lie on any optimal alignment, against 81-89% in the flanks. Each false one
  forces about 1,400-1,600 extra edits.
- Chaining MUMs MUMmer-style (heaviest chain) is worse than MC: 1.65, and 1.83 with recursion.
- Charging each link its diagonal shift fixes most of it: 1.026-1.034 at the top level, 1.040-1.063 with recursion.
- The windows left inside arrays stay at kb scale: median 0.8-1.3 kb, p90 3.5-5 kb, maximum 18 kb.
- Order: longest first was best for progressive POA in 5 of 6 tests, CHM13 first was worst, and random orders move
  cost/opt by 0.11-0.33.

So anchors must be scored by what they cost, not by how long they are. The DP inside arrays must handle kb-scale
windows with large, known diagonal shifts.

## 4. The algorithm

Sequences join the graph one at a time. Each is anchored pairwise against what is already aligned. The backbone is
chained optimally from sporos's matches first. Only then are the gaps between anchors filled with karan, and the
result fused into the graph.

### 4.1 The graph

- One node per base.
- Each column is a group holding at most one node per base, as in SPOA and abPOA.
- Paths are stored explicitly, as arrays of nodes, never rebuilt from edge labels.
- Each edge carries a bitset of the haplotypes that use it.
- The topological order is kept incrementally and updated per window, not re-sorted per sequence.
- Identical input sequences are aligned once and carry their count as a weight.
- Output:
  - MSA rows, named as the input;
  - GFA with compacted nodes and exact paths;
  - pgrealign's `msas.jsonl.gz`.

### 4.2 Order and the anchor partner

- **Start the graph:** the starting sequence is the first thing E6 decides. The candidates are the longest
  sequence, the medoid (least haplotype-weighted multiset 15-mer Jaccard distance), and the highest estimated copy
  number.
- **Add the rest:** in descending length, or Prim order nearest the graph.
- **Anchor partner:** each new sequence is anchored against one path:
  - first, a fixed reference path (the start sequence), which keeps every haplotype in one coordinate system;
  - later, the nearest path;
  - later still, every path with node-distinct uniqueness, following minichain's haplotype-aware chaining with a
    recombination penalty.
- Never anchor against the consensus inside arrays, and never against the previous sequence.

### 4.3 Anchors (sporos)

For a query Q and partner path P at level ℓ with minimum length L_ℓ:
- **MUMs:** `Reference.mems(mode="mum")`, unique in P and in Q. They enter at full weight.
- **Rare MEMs:** `mode="maxmatch"` with `max_occurrences = c_max`.
  - Every occurrence is a candidate.
  - Its weight is length × (count_P · count_Q)^(−1/2), as centrolign does.
  - The chain picks among occurrences by colinearity.
- **SMEMs** find candidate intervals cheaply but never become partition points by length alone. In an array the
  longest SMEM for a copy usually sits on a shifted diagonal. A SMEM's interval is re-seeded as all MEM
  occurrences inside it, as BWA-MEM does.
- **Frequent matches** (above c_max) are dropped at this level. They return in the recursion, counted again inside
  the window.

### 4.4 The optimal backbone: chaining

Colinear chaining with overlaps and gap costs, after ChainX (Jain, Gibney and Thankachan 2022). Linking a′ to a
costs:

    connect(a′, a) = σ·s + o + γ(d) + R·[path switch]

where:
- s = min(g_Q, g_P) is the substitution part of the gap;
- o = |o_Q − o_P| is the diagonal shift inside an overlap, since overlaps are allowed (consecutive copies always
  overlap);
- d = g_Q − g_P is the indel imbalance;
- γ is concave, a minimum of two or three affine pieces, so one long copy-number gap is cheaper than many small
  ones. A *period-aware* term ρ·min(|d| mod p, p − |d| mod p) makes a shift by whole copies cheap; p comes from
  the repeat catalogue or SMEM self-hits.

The chain maximises the summed anchor weights minus the summed link costs.

Ties are broken in graph coordinates: among equal chains, the indel goes leftmost. Equivalent copy-number gaps then
land in the same columns for every haplotype, which is what pairwise consistency and genotyping reward.

Guarantee:
- With σ = 1, ρ = 0, linear γ and all MEM occurrences of length at least L, the chain cost equals the anchored edit
  distance (ChainX, Theorem 2).
- The final alignment then lies between ED and ED + m_<L, where m_<L counts matched bases in exact runs shorter
  than L (*inferred* bound; E2 checks it).
- Caps, SMEMs and the period term trade that bound for speed and repeat sense. The alignment always stays optimal
  given its anchors, and never worse than the chain.

Algorithm:
- First, ChainX's O(n·OPT) guess-and-quadruple loop. It is exact when it terminates, and about 50 lines.
- Then one range-maximum structure per gap piece and per path, as minichain and centrolign do, when the anchor
  count calls for it.

### 4.5 Recursion

For each pair of consecutive chained anchors:
1. Take the query interval and the subgraph between the anchors' nodes.
2. If the window fits karan's budget W_DP, or its band (|d| plus a margin) is cheap, stop.
3. Otherwise:
   - re-seed inside the window with a smaller L_{ℓ+1} (for example 41 → 31 → 21 → 15);
   - count uniqueness again in the window;
   - chain between the two fixed ends;
   - recurse.

Top-level anchors must be robust, because each level freezes the one above. At level 0 an anchor counts as robust
only when removing it costs the chain more than a margin. In practice level 0 is mostly flanks and array junctions.

Where an array window has no safe anchors, do not force any. Fill it with full banded DP first, the measured best.
Unit-level DP comes later: tokenise copies, align the unit strings, and fill each unit pair with karan.

### 4.6 Filling the gaps (karan)

- **Each window is a global alignment pinned at both anchors,** inside the band `Band(min(0, d) − w, max(0, d) + w)`.
  Its centre is the diagonal shift d known from the chain, and w is the slack for substitutions and small indels.
  Band doubling gives exactness: widen until the cost proves the band sufficient. That is karan's `max_cost` and
  POASTA's certificate.
- **Costs:** two-piece gap-affine (`Costs.two_piece`), with the second piece cheap enough that a copy-number
  difference is one gap. The candidates are pgrealign's scaled abPOA set and MC's 400/30, 1200/1. `Ties.LEFT`
  matches the chain's rule.
- **Each window is aligned against a path segment:** karan is sequence to sequence. Every distinct path sequence
  through the window's subgraph is a candidate. All of them go through karan's batch call and the cheapest is kept.
  - The new sequence may switch haplotype from window to window, approximating sequence-to-graph alignment.
  - That is exact when the best route through the window is an existing haplotype.
- **Later, a graph kernel inside the window,** so that extra copies inserted at different slots can merge
  (pgrealign's cross-slot scatter). The candidates are an abPOA-style SIMD DP with the exact per-window band, and a
  POASTA/Theseus-style exact search with a score cap.

### 4.7 Fusing into the graph

- Walk the window CIGARs along the chosen path's nodes:
  - a match joins the existing node;
  - a mismatch joins or extends the column's group;
  - an insertion becomes new nodes;
  - a deletion becomes an edge.
- After each window, merge parallel insertion branches that spell the same sequence.
- After the last sequence, verify that every path spells its input exactly and that the path count is unchanged.

### 4.8 Refinement

One remove-and-realign pass per sequence, against the finished graph, removes most of the order dependence.

## 5. Changes to sporos and karan

- **karan:** moved to the shared nightly (done, `016ca6f`). Nothing else is needed for the first version:
  - `Band(low, high)` already takes an asymmetric band;
  - `Costs.two_piece`, `Ties.LEFT`, `max_cost` and batches exist.
  - Later, for the window kernel: a reusable `Aligner` per thread to keep memory flat, and a GPU batch for many
    windows.
- **sporos:**
  - `Reference.mems` already gives MUMs, `mumreference` and capped MEMs between a path and a query.
  - What plegma will need:
    - (a) cheap seeding of many small windows without building a full index each time, through a
      substring-of-a-text API or a direct exact-match finder for short pairs;
    - (b) MEM occurrence counts in both sequences on each match, so weights need no second pass;
    - (c) later, node-distinct counting over many paths.
  - Each is added when an experiment needs it, with tests in sporos first.

## 6. Experiments

All runs use the vntr-eval packages: an 11-region pilot, then all 149 (hap32, then the full panel), the 45-region
hard set, and held-out chr6.
- Scoring: `tools/msa_graph.py` + `tools/evaluate.py --no-truvari` (cost/opt, affine cost/opt, `kmer_frac_extra`,
  SV pieces), plus qc-style pair F1 and precision, gap runs per row, time and peak memory.
- Baselines: MC, abPOA longest first, SPOA longest first, abPOA `-S`, the pgrealign medoid star, and unit-aware.

| | experiment | decides |
| :-- | :-- | :-- |
| E0 | MUM census (done, above) | anchors must pay for diagonal shifts |
| E1 | anchor census per type (MUM, mumreference, MEM capped at 2/8/50, SMEM) × L ∈ {15, 21, 31, 51}: share on some optimal alignment, inside and outside arrays, SMEM off-diagonal share | level-0 anchor types, c_max, L schedule |
| E2 | pairwise chain quality: ChainX cost, minimap2-like, centrolign-like weights, on {MUM}, {MUM + SMEM}, {MUM + rare MEM}, {all MEM}; check Theorem 2 and the additive bound on small regions | the chainer |
| E3 | gap cost: linear, concave, concave + period term | γ and ρ |
| E4 | recursion: depth 1 vs LAGAN-style, W_DP ∈ {1, 2, 5, 10} kb | window budget |
| E5 | pairwise end to end, chain + karan fill vs karan unanchored | gate: cost ≤ 1.001 × optimum outside arrays, ≤ 1.01 overall, faster and smaller |
| E6 | anchor partner in progressive POA: fixed longest, fixed medoid, nearest path, all paths | §4.2 |
| E7 | order: longest, copy number, medoid + Prim, CHM13, three random seeds | order; anchoring to a fixed path should shrink the random spread (*inferred*) |
| E8 | tie-breaking: left, right, none | gap stacking (pair F1, gap runs) |
| E9 | anchorless arrays: banded DP vs unit-level DP vs forced SMEMs (negative control) | the array fallback |
| E10 | scale: full panel, regions above 20 kb, L012272 at 115 kb | memory within 12 GB where abPOA failed |
| E11 | hard set, chr6, then pgrealign's replace/verify/qc on chr20 | acceptance: beat 1.067 and 0.703, no control worse by 0.005 |

E1 to E5 are pairwise and need no graph. They come first, and build the anchor and fill layers that the POA then
calls.
