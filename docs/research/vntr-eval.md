# vntr-eval and pgrealign: a briefing for plegma

Sources: `benedictpaten/vntr-eval` at `da5c634` and `benedictpaten/pgrealign` at `519ca6a`, cloned into
`scratchpad/src/`. File paths below are relative to those clones. Statistics marked **[computed]** were
computed here by `scratchpad/scripts/region_stats.py` (sequence statistics) and
`scratchpad/scripts/mumcheck.c` + `mum_driver.py` + `mum_summary.py` (the pairwise MUM experiment). No code
from either repository was executed. Statements marked **[inferred]** are my reading, not a claim made in
the repositories.

---

## 0. TL;DR

- **Test bed.** 149 self-contained region packages (`regions/<id>/`). Each holds `hap32.fa` (CHM13 + GRCh38 +
  32 sampled HPRC haplotypes, anchor to anchor, so every sequence shares the same first and last base:
  **global alignment**), `hprc.fa.gz` (the full HG002-free panel, up to 457 haplotypes), `mc.gfa` (the
  Minigraph-Cactus baseline), `truth.fa` (HG002, evaluation only), `calls.tsv`, `region.json`. The 40
  `hotspot_vntr` regions are the test set; three control strata must not regress.
- **Headline metric.** `all_cost_over_opt` = sum over all hap32 pairs of the **unit (Levenshtein) edit
  cost** of the pairwise alignment the graph implies (shared nodes are forced matches, the rest aligned
  optimally) over the sum of the pairs' optimal edit distances. A gap-affine version (mismatch 4, gap of k =
  6 + 2k) is reported beside it. MC: 1.566 at hotspot VNTRs (2/40 at or below 1.1), 1.102 at matched
  controls. Target: hotspots **<= 1.1**, controls unchanged. Best so far on hap32: the repeat-unit-aware
  aligner, 1.048 (33/40 <= 1.1).
- **But cost/opt does not predict calling.** The later, whole-contig work (Stage 4, pgrealign) found that
  k-mer redundancy, gap-run fragmentation and above all **pairwise-consistency precision/F1 against
  direct optimal pairwise alignments** track `vg call` SV accuracy; cost/opt is blind to a deletion cut
  into pieces. The production method became **a medoid star**: every distinct allele aligned pairwise
  (abPOA, unbanded, global) to the panel medoid, left-normalised, merged on the medoid's columns, with
  same-slot insertions aligned to each other. It beat every progressive aligner, the unit-aware aligner,
  FAMSA, mafft, MST threading and backbone variants on calling.
- **Plug-in point.** pgrealign's `replace` reads only `id`, `status`, `rows` from `msas.jsonl.gz`, and
  `star.star_msa` takes any `align(seqs, banded)` callable. Either is a clean slot for plegma.
- **Pairwise MUMs in VNTR arrays [computed].** They are *not* rare: a hotspot haplotype vs CHM13 has a
  median of 18 MUMs (>= 20 bp) inside the array. But most of them are **not on any optimal alignment**
  (14-21% of array MUMs are, against 81-89% of flank MUMs), and forcing the MUMmer-style heaviest chain
  makes the pairwise cost *worse than MC* (1.65-1.83 pooled). A chain that charges for diagonal shifts
  brings the loss down to 3-6% (1.03-1.06), but the unanchored windows left inside arrays stay at
  kb scale (median ~1 kb, p90 ~4-5 kb, max ~18 kb). Flanks are safe anchors (section 5.4).

---

## 1. The test data

### 1.1 What the regions are (`regions/README.md`, `docs/findings.md`)

In the HPRC v2.1 Minigraph-Cactus (MC) graph, ~1% of VNTR regions hold about half of short-read
`vg call`'s SV false positives. There the graph does not tangle (no cycles); it **misaligns** the panel:
repeat arrays of different haplotypes are stacked out of phase in one giant bubble of ~3.6 bp nodes
(848 nodes/kb vs 113 at matched VNTRs), expansions are placed at different ends of the array and repaid
by compensating deletions (e.g. a haplotype 9 bp longer than CHM13 encoded as 1,430 bp inserted +
1,421 bp deleted), and some haplotype families run through parallel copies 96.7% identical to the CHM13
sequence they bypass (`docs/findings.md` section 2).

Each region is cut between two **anchor nodes**: CHM13 nodes that every panel path passes once, untouched
by HG002 truth variants, at least 200 bp outside the locus (median flank 222 bp; up to 42 kb). So all
sequences of a region start with the same left-anchor base(s) and end with the same right-anchor
base(s): "a global (end-to-end) alignment is the right model. Do not use local alignment"
(`regions/README.md`).

### 1.2 Strata (`regions/README.md`)

| stratum | n | role | definition |
|---|---|---|---|
| `hotspot_vntr` | 40 | **test** | VNTRs (period >= 7) with the most `vg call` SV FPs (>= 12) |
| `hotspot_other` | 24 | test (secondary) | non-VNTR loci with >= 3 FPs: mostly STRs (period < 7), a few LC / satellite / non-repeat |
| `control_vntr_matched` | 40 | must not regress | FP- and FN-free VNTRs matched to the hotspots on length and period bin |
| `control_vntr_correct` | 25 | must not regress | VNTRs with a truth SV that vg calls correctly |
| `control_nontr_sv` | 20 | must not regress | non-repeat loci with a correctly called truth SV |

150 loci were selected; L000601 could not be packaged (a ~94 kb copy-number locus with no usable
anchors). Nine pilot regions are flagged in `notes` (`pilot:`): misaligned (L009656, L014297, L015415,
L005990), under-aligned (L012184, L011138, L012272), alleles absent from the panel (L002013), negative
control (L016870). All loci are in the HG002 T2T-Q100 v1.1 SV benchmark, chr1-22 + X.

### 1.3 Files per region (`regions/README.md`)

| file | content |
|---|---|
| `hap32.fa` | **the input.** One record per hap32 panel path spanning anchor to anchor; one line per sequence, upper case. CHM13 first (`CHM13#0#chrN`), then `GRCh38#0#chrN[offset]`, then `recombination#k#chrN#f`. Identical sequences kept as separate records (allele counts matter). `len=` in the description. |
| `hap32.fragments.fa` | panel path pieces that enter but do not span (`:enters_L`, `:exits_R`, `:internal`); optional |
| `hprc.fa.gz` | the same span for the full HG002-free HPRC v2.1 eval panel: CHM13, GRCh38, then `hprc#1..` (no sample names), with `group=`, `group_size=`, `same_as_hap32=` |
| `mc.gfa` | baseline M0: MC's **full-panel** alignment restricted to the hap32 paths (S, L, one P line per hap32 record). Not an alignment of the 34 sequences alone. |
| `truth.fa` | `HG002#1`, `HG002#2` over the span. **Evaluation only; never use for building or tuning.** |
| `calls.tsv` | truth / vg / PanGenie records over the span, with truvari status (context only) |
| `region.json` | the `regions.tsv` fields + anchor sequences, flank lengths, `tr_annotation` (census period, motif, copies, GC, CHM13 and truth tandem scans, RepeatMasker), hap32 stats, MC graph stats, provenance; truth-derived fields are held out |
| `regions.tsv` | the index: `region_id, stratum, contig, core_start, core_end, span_start, span_end, anchor_left, anchor_right, period, motif, copies, n_hap32, n_hprc, vg_fp, vg_fn, pg_fp, pg_fn, in_benchmark, notes` (CHM13v2.0, 1-based inclusive; core = the VNTR interval, span = anchor to anchor = exactly the CHM13 sequence in hap32.fa) |
| `vntr_ge1kb.bed` | 3,111 benchmarked VNTRs >= 1 kb (no sequences): the scale-up target, 46% of all vg SV FPs |

Both panels are HG002/HG003/HG004-free and exclude the `gref_CHM13` cover path. Do not add those
samples from anywhere.

### 1.4 Sequence statistics [computed]

`scratchpad/out/region_stats.tsv` has one row per region. Cells are **min / median / p90 / max over
regions**. "full" = `hprc.fa.gz` records (includes CHM13 and GRCh38).

| stratum | n | span (= CHM13) bp | core bp | period | copies | distinct hap32 seqs | hap32 median len | hap32 max len | hap32 max-min | full-panel records | distinct full seqs | full max len | distinct full-panel bp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hotspot_vntr | 40 | 1,419 / 5,709 / 15,612 / 55,376 | 890 / 4,824 / 12,312 / 18,884 | 8 / 32 / 120 / 192 | 6 / 153 / 470 / 996 | 13 / 30 / 34 / 34 | 1,566 / 4,664 / 11,768 / 51,676 | 1,982 / 8,518 / 20,475 / 55,376 | 1,114 / 5,304 / 17,796 / 32,300 | 373 / 455 / 457 / 459 | 166 / 375 / 442 / 450 | 3,903 / 12,428 / 42,648 / 115,456 | 0.27 / 1.95 / 4.57 / 21.6 Mb |
| hotspot_other | 24 | 859 / 2,652 / 4,339 / 6,619 | 422 / 1,916 / 3,855 / 6,054 | 0 / 4 / 6 / 6 | 68 / 587 / 1,514 / 1,933 | 12 / 29 / 32 / 34 | 919 / 2,459 / 3,268 / 4,062 | 1,295 / 3,488 / 6,543 / 11,775 | 223 / 2,076 / 4,612 / 10,080 | 340 / 454 / 457 / 457 | 114 / 280 / 418 / 442 | 1,751 / 4,286 / 11,716 / 29,811 | 0.17 / 0.77 / 1.36 / 1.87 Mb |
| control_vntr_matched | 40 | 1,313 / 3,669 / 7,379 / 19,353 | 856 / 3,129 / 6,780 / 10,047 | 7 / 36 / 177 / 1,990 | 4 / 88 / 458 / 698 | 2 / 8 / 20 / 33 | 1,313 / 3,586 / 7,419 / 19,029 | 1,313 / 4,082 / 7,881 / 19,353 | 0 / 750 / 2,624 / 17,270 | 278 / 455 / 458 / 458 | 34 / 82 / 165 / 411 | 1,505 / 4,342 / 10,710 / 21,229 | 0.05 / 0.29 / 1.36 / 2.82 Mb |
| control_vntr_correct | 25 | 1,172 / 2,722 / 4,587 / 6,844 | 672 / 2,111 / 3,073 / 4,296 | 10 / 42 / 782 / 1,532 | 2 / 49 / 145 / 211 | 2 / 8 / 15 / 28 | 558 / 2,747 / 3,683 / 4,339 | 1,284 / 2,931 / 4,616 / 6,844 | 70 / 527 / 3,226 / 4,931 | 343 / 453 / 456 / 458 | 15 / 59 / 182 / 273 | 1,340 / 3,633 / 6,631 / 15,836 | 0.02 / 0.21 / 0.54 / 0.85 Mb |
| control_nontr_sv | 20 | 658 / 816 / 1,129 / 2,183 | 201 / 206 / 529 / 1,664 | - | - | 2 / 4 / 10 / 12 | 685 / 906 / 1,675 / 2,271 | 779 / 1,051 / 1,860 / 3,714 | 50 / 308 / 921 / 2,886 | 269 / 453 / 457 / 458 | 7 / 20 / 30 / 51 | 780 / 1,056 / 1,860 / 3,714 | 7 / 15 / 66 / 99 kb |
| all | 149 | 658 / 3,009 / 7,419 / 55,376 | 201 / 2,462 / 6,679 / 18,884 | | | 2 / 13 / 32 / 34 | 558 / 2,919 / 6,982 / 51,676 | 779 / 3,714 / 11,726 / 55,376 | 0 / 1,578 / 9,361 / 32,300 | 269 / 454 / 457 / 459 | 7 / 141 / 422 / 450 | 780 / 4,410 / 20,298 / 115,456 | |

- Totals: 5,026 hap32 records, 65,199 full-panel records (match the README). 136 regions have 34 hap32
  records; 13 have fewer (24-33). The README says 17 regions have `n_hap32` < 34; I count 13 from
  `hap32.fa` (small discrepancy, not investigated). 8 regions have a non-empty `hap32.fragments.fa`. One
  region has an N run in hap32 (L012184, 3,146 bp of N in one haplotype).
- **Periods (census, `region.json`):** hotspot VNTRs: 8 at 7-15 bp, 11 at 16-30, 14 at 31-60, 7 at
  61-200 (median 32 bp, max 192); copies median 153 (max 996, L009000: 14 bp x 996). hotspot_other are
  all period < 7 (STRs; median 4), with up to ~1,900 copies. Matched controls span 7-1,990 bp periods.
- **Hotspot alleles are long and length-polymorphic.** Median longest hap32 allele 8.5 kb; 18/40
  hotspots have a hap32 allele > 10 kb and 5/40 > 20 kb. The longest hap32 allele is a median 1.29x
  CHM13's length (max 5.1x). Within-region hap32 length spread is a median 5.3 kb.
- **The full panel is much bigger and has outliers.** A hotspot has a median 375 distinct full-panel
  sequences (vs 30 in hap32), and the longest full-panel allele is often far longer than any hap32
  allele: L012272 115 kb (hap32 max 55 kb), L015347 61 kb, L001084 56 kb, L007009 43 kb, L007040 43 kb,
  L012184 39.6 kb (hap32 max 9.7 kb, CHM13 4.8 kb), L001909 38 kb, L008976 38 kb.

### 1.5 The hard set (`hardset/`)

`hardset/chr20.tsv` lists 45 chr20 tandem-repeat regions of the **whole-contig pgrealign build** (not the
149 packages) whose realigned **pair F1 < 0.8** (pgrealign `qc`; median over all 14,511 chr20 regions is
1.0): 27 VNTR, 7 STR+VNTR, 8 STR, 1 each LC, LC+VNTR, SAT; median 220 distinct alleles (5-452), longest
allele median 4.3 kb (1.3-45.8 kb). 8 of them hold 9 of the 10 remaining local ONT phase-flip stretches.
**The allele FASTAs (`work/hardset/chr20/fasta/<id>.fa`) are git-ignored and not in the clone**; only
the table (ids, CHM13 coordinates, boundary node ids, motif, allele counts, F1 before/after) is. To use
the hard set plegma needs the HPRC chr20 graph and `pgrealign extract` (or `tools/hardset.py`).

Baseline: original graph pair F1 median 0.647 (LQ 0.509); medoid star 0.703 (LQ 0.628), 1,961 s total,
<= 601 s per region. Strategy benchmark in section 3.4.

---

## 2. The metrics

### 2.1 Stage 0 graph metrics (`tools/METRICS.md`, `tools/evaluate.py`)

A candidate is a **GFA 1.0** with one P (or W) line per `hap32.fa` record, same name, spelling exactly its
sequence (case-insensitive; N runs must be nodes of their own). Extra paths are ignored. Node ids,
lengths and compaction are free. Invalid graphs (name mismatch, wrong spelling, unknown node) get
`status: invalid`, exit 3, no metrics. The graph should be a DAG on the used handles.

**`alignment.all_cost_over_opt` ("cost/opt").** For each unordered pair of hap32 paths A, B (561 pairs
for 34; above `--max-pairs` 600 a seeded sample that keeps every CHM13 pair):
- the handles A and B share are forced matches (free);
- between consecutive shared handles, the two segments are aligned optimally;
- `cost` = sum of those segment edit distances: the cheapest alignment that stacks every shared node;
- `dopt` = unconstrained optimal **global unit-cost (Levenshtein)** edit distance (Myers/Hyyro
  bit-vector, `ed_unit` in `tools/c/vntreval.c`);
- `all_cost_over_opt` = sum(cost) / sum(dopt) over all pairs (pooled, not a mean of ratios).

**Gap-affine variant** (`affine_all_cost_over_opt`): mismatch 4, a gap of k bases 6 + 2k, exact banded
Gotoh with a provably-global doubling band, on all pairs within a 5e10-cell budget (else CHM13 pairs + a
sample). So: **the headline is unit cost; the affine model is (X=4, O=6, E=2), match 0.**

Other alignment fields: `all_excess_per_kb` (1000 x sum(cost - dopt) / sum of mean pair lengths; MC
182 vs 4.4 at controls), `all_U_per_kb` (homologous bp left on parallel nodes), per-pair ratio quantiles,
`all_frac_pairs_optimal`, and `ref_*` versions over CHM13 pairs only.

Other sections:
- **size**: nodes, edges, nodes/kb, node length stats, fraction of 1-bp nodes, non-reference bp,
  distinct walks vs distinct sequences (same sequence on different walks = misalignment), top-level
  bubbles (cuts in topological order = vg's top-level snarls), SV bubbles, largest bubble.
- **inflation**: each path decomposed against CHM13 the way vg deconstruct writes it; `sv_pieces_*`
  (records with |len diff| >= 50 per path; MC ~8-9 at hotspots vs 1 at controls), `indel_bp_ratio`
  (graph indel bp / net indel bp; 1 = no compensating INS/DEL pairs), `anchored_over_free`.
- **redundancy**: `kmer_frac_extra` = over all 21-mers, (graph positions - max copies in one path) /
  positions. MC 0.60 at hotspots vs 0.05 at controls. **It predicted read placement best** (Spearman
  0.75 with change in core MAPQ<5 in Stage 2, vs 0.26 for cost/opt).
- **truth**: `h1/h2_d_graph` = exact min edit distance of each HG002 haplotype to any source-to-sink path
  (recombinant paths allowed); `truvari.*` = the truth written through the graph, scored raw,
  `refine` (phab) and whole-span phab. Needs the big data (`VNTR_EVAL_DIR`); report `needs_big_data`
  without it, or pass `--no-truvari`.
- **reads** (`--reads`, Stage 1): placement redundancy of the region's existing short reads (needs the
  big data).
- **`--panel FASTA`**: scores a full-panel graph against `hprc.fa.gz` on 300 sampled pairs.

### 2.2 The target and the MC baseline (`docs/results.md`, `results/stage01_tables.md`)

Target (`docs/plan.md`): **hotspots move to the matched controls' distribution, cost/opt <= 1.1, and the
controls do not change.** Later added: k-mer redundancy toward the controls' ~0.05.

Median cost/opt per stratum (regions <= 1.1 in brackets), hap32 arm (`docs/results.md` section 1):

| method | hotspot VNTR (40) | hotspot other (24) | matched ctrl (40) | correct ctrl (25) | non-repeat SV (20) |
|---|---|---|---|---|---|
| **MC (`mc.gfa`)** | **1.566 (2)** | **1.625 (3)** | **1.102 (19)** | **1.015 (17)** | **1.000 (20)** |
| mafft FFT-NS-2 | 1.327 (4) | 1.313 (6) | 1.100 (20) | 1.057 (15) | 1.000 (20) |
| mafft L-INS-i | 1.111 (18, n=39) | 1.164 (8) | 1.028 (33) | 1.015 (22) | 1.000 (19) |
| mafft E-INS-i | 1.182 (9, n=37) | 1.221 (6) | 1.030 (33) | 1.007 (22) | 1.000 (20) |
| mafft G-INS-i | 1.110 (15, n=39) | 1.163 (8) | 1.027 (34) | 1.034 (21) | 1.000 (19) |
| abPOA (global, longest first) | 1.083 (25, n=37) | 1.150 (9) | 1.024 (33) | 1.006 (22) | 1.000 (20) |
| spoa (global, longest first) | 1.067 (20, n=27) | 1.125 (10) | 1.016 (34, n=39) | 1.005 (24) | 1.000 (20) |
| abPOA, Cactus BAR settings | 1.170 (11) | 1.284 (6) | 1.031 (31) | 1.007 (22) | 1.000 (20) |
| repeat-unit-aware | **1.048 (33)** | **1.053 (20)** | 1.015 (37) | 1.015 (25) | 1.000 (19) |

Other MC baseline numbers at hotspot VNTRs: gap-affine cost/opt 1.475 (good methods 1.09-1.15),
excess 182 edits/kb (controls 4.4), U 40.5 bp/kb, SV pieces per path ~8-9.2 (floor ~1), `kmer_frac_extra`
0.60 (good methods 0.40-0.52; controls 0.05), nodes/kb 764 (MSA graphs 265-381; MC's count includes cuts
made by the other ~420 haplotypes). Full MC scored on its own full panel: 1.513.

### 2.3 Calling-level metrics (Stages 2-4)

- **Stage 2** (19 regions): local giraffe re-map: core MAPQ<5 fraction (MC 0.817, L-INS-i 0.671, HG002's
  own graph 0.518 = the tandem-ambiguity floor).
- **Stage 3**: local `vg call`, scored as the **called haplotype pair's edit distance (ED) to HG002**
  (`tools/score_haplotypes.py`, best pairing), 3 replicates, plus raw / refined / phab truvari.
- **Stage 4 / pgrealign**: patch whole chr20 (and held-out chr6), re-map, re-call; truvari SV F1 raw and
  `refine`, aardvark small variants, per-region paired bootstrap of refined FP+FN.
- **pgrealign `qc` pair F1** (used for the hard set and whole-contig QC): for up to 20 random allele pairs
  (+ CHM13 vs 20 heaviest), the aligned base pairs of the MSA's projected pairwise alignment vs those of
  the pair's own optimal (unbanded, left-normalised abPOA, pgrealign scores) alignment: precision =
  tp/msa, recall = tp/opt, **F1 = 2tp/(opt+msa)**; plus `growth` (graph bases / longest allele),
  `forced` (40-column windows with >30% mismatches = forced homology), `ref_sv` (>= 50-column gap runs).
  This penalises co-optimal gaps placed elsewhere, unlike cost/opt.

**Important lesson (Stage 4):** "Stage-0 metrics do not predict calling here" (`results/stage4_chr20.md`
4p). cost/opt is blind to fragmentation ("one deletion and the same deletion in pieces cost the same
edits", 4t). Better predictors: per-region **pair precision** vs direct optimal pairwise alignments (4u,
monotone dose-response with error), and **gap runs** of the hap32 rows (4t). The adopted method was chosen
on pairwise consistency F1 and calling.

### 2.4 Exactly how to evaluate a plegma MSA (from the vntr-eval root)

```bash
make -C tools/c                                     # optional; evaluate.py builds tools/bin/libvntreval.so itself
python3 tools/evaluate.py regions/L014297 regions/L014297/mc.gfa --no-truvari      # baseline -> results/mc/L014297.json

# arm 1: hap32 alone
plegma regions/L014297/hap32.fa > L014297.msa.fa          # rows named as hap32.fa (first word), '-' gaps
mkdir -p candidates/plegma                                 # msa_graph.py's CLI needs the directory to exist (tools/POA.md)
python3 tools/msa_graph.py L014297.msa.fa regions/L014297/hap32.fa candidates/plegma/L014297.gfa [--engine native]
python3 tools/evaluate.py regions/L014297 candidates/plegma/L014297.gfa --no-truvari
#   -> results/plegma/L014297.json  (method name = the GFA's directory name; --out to redirect, --skip truth to skip the truth DP)

# arm 2: full panel, projected onto hap32
python3 tools/panel.py union regions/L014297 work/union/L014297.fa        # distinct seqs u0001.. (CHM13 first) + .map.tsv
plegma work/union/L014297.fa > L014297.full.msa.fa                          # keep the u0001.. names
python3 tools/panel.py project L014297.full.msa.fa work/union/L014297.map.tsv regions/L014297/hap32.fa \
    candidates/plegma__all/L014297.msa.fa
python3 tools/msa_graph.py candidates/plegma__all/L014297.msa.fa regions/L014297/hap32.fa candidates/plegma__all/L014297.gfa
python3 tools/evaluate.py regions/L014297 candidates/plegma__all/L014297.gfa --no-truvari
# full-panel graph itself: tools/panel.py panel-graph FULL.msa.fa MAP.tsv regions/<id> OUT.gfa;
#                          tools/evaluate.py regions/<id> OUT.gfa --panel regions/<id>/hprc.fa.gz

python3 tools/evaluate_all.py plegma --jobs 2 --threads 2                   # every candidates/plegma/<id>.gfa; no --no-truvari flag,
#                                                                          so truvari reports needs_big_data without VNTR_EVAL_DIR
python3 tools/summarise.py --methods mc,plegma --md results/summary.md     # per-stratum table vs mc, Wilcoxon
```

- `msa_graph.py`: column induction (rows with the same base in a column share a 1-bp node), then unchop
  (`vg mod -u` equivalent), N treated as its own symbol, all-gap columns dropped, every path re-checked.
  `--engine native` needs no vg. `--merge-blocks` optionally collapses SNV mesh.
- Requirements: Python >= 3.9 stdlib, `make` + a C compiler. No vg, no mafft. Typical evaluation 3-15 s
  per region, 1-3 min for the largest (L015347, L012272); optima are cached in
  `$VNTR_WORK/cache/evaluate/<id>/optima.json`.
- The gold standard for calling (Stage 3-4) needs the big data (`~/PycharmProjects/vg-call-eval`) and
  cannot be run from the clone. **[inferred]** For plegma's first loop, report cost/opt, affine
  cost/opt, `kmer_frac_extra`, SV pieces, and also compute qc-style pair F1 against direct pairwise
  optima (pgrealign `qc.region_qc` can do this from Python on package dicts; see section 4.4).

---

## 3. What has been tried

### 3.1 Stage 0-2 realigners on the 149 packages (`docs/results.md`, `tools/POA.md`, `tools/UNIT_ALIGN.md`)

Every MSA was turned into a graph by `tools/msa_graph.py`; 3,488 graphs scored, 0 invalid.

| method | configuration | notes |
|---|---|---|
| mafft FFT-NS-2 | fast progressive | regresses 10/40 matched controls (worst L000709 1.00 -> 2.97); recreates parallel copies (L012184 1.86, L002143 9.5) |
| mafft L-INS-i / E-INS-i / G-INS-i | iterative, consistency | L-INS-i: the only method that improved Stage 2 read placement at every hotspot (8/8, MAPQ<5 0.817 -> 0.671), fewest SV pieces (3.7-4.1), lowest-ish k-mer redundancy (0.42); **damages non-repeat SVs** (L006556 1.00 -> 1.19, L002292, L003883). Cannot run the full panel at hotspots (finished 4/1/1 hotspot VNTRs; predicted up to 244 GB, 515 h). |
| abPOA 1.5.7 | `-m 0 -r 1`, defaults (M 2, X 4, O 4,24, E 2,1, band `-b 10 -f 0.01`), **longest sequence first** | 0.1 s median on hap32. `-S` (minimizer seeding) **fails on VNTRs**: anchors chain wrongly, haplotypes left side by side (92,688 columns vs 13,051 at L012184). Unbanded no better. `-p` guide tree: order-insensitive but worse. Regresses L011224 (10 bp CA repeat). Memory 18-29 B x max_len^2. |
| spoa 4.1.5 | `-l 1 -r 1` (global!), defaults, longest first | best POA on hap32 cost/opt; memory 80 B x max_len^2, so > 12.2 kb sequences refused (14 hap32 regions, 26 full-panel) |
| abPOA, Cactus BAR settings (`poa_abpoa_mc`) | HOXD70-like matrix, `-O 400,1200 -E 30,1 -b 1000 -f 0.1 -p`, 10 kb end rule emulated | already beats MC on the full panel (1.515 -> 1.269); so most of MC's damage comes before BAR (CAF pinches, inferred by the authors) |
| repeat-unit-aware (`unit_aware`, M4) | `tools/realign_units.py` + `tools/unit_dp.c`, see 3.2 | **best Stage 0** (1.048, 33/40); the only method that finishes the full panel everywhere (149/149, median 13 s, max 1,722 s, 9.4 GB) |

**Input order is part of a POA method** (`tools/POA.md`): across random orders cost/opt moves by 0.11-0.19
(abPOA) and 0.23-0.33 (spoa) on hap32, up to 0.94 on the full panel; a bad order can be worse than MC
(spoa random order at L005990 1.419 vs MC 1.211; at L014297 2.58 vs 1.30). **Longest first was best in
5/6 region-tool pairs on hap32 and 4/6 on the full panel**; CHM13-first ("given") is worst; Prim
nearest-first from the medoid ("guide") is close. Explanation offered: once the longest array is in the
graph, later haplotypes align to it with deletions instead of each inserting extra units somewhere.
In the full-panel iteration (`results/stage4_chr20.md` 4h) frequency-first, shortest-first and seeded
orders all hurt (cost/opt 1.3-1.5 mean).

**Equal-cost alignments disagree**: abPOA, spoa and L-INS-i share only 18-31% of their CHM13-vs-haplotype
residue pairs at 4 of 6 pilot hotspots although their cost/opt is within 0.07 (`tools/POA.md`).

**The panel question.** Aligning the full panel and projecting onto hap32 (the fair comparison with MC)
keeps 98% of unit-aware's hotspot gain (1.059), 91% of abPOA's (1.130), 88% of spoa's (1.115), but the
controls keep only 23-66% of their small gain. Projection **fragments deletions**: hap32 rows inherit
junctions supported by absent panel alleles (L014297 spoa: indel runs per haplotype 28.5 -> 45.4, mean
37 -> 24 bp).

Runtime on hap32 (`results/runtime_summary.md`, shared, loaded 10-core Mac; wall clock): medians abPOA
0.1 s, FFT-NS-2 0.5 s, spoa 0.9 s, Cactus abPOA 1.1 s, unit-aware 2.0 s, L/E/G-INS-i 3-4 s (hotspot
medians 1.2 / 5 / 5.2 / 4.3 / 5.5 / ~47 s). Full panel (104 shared regions): abPOA 0.5 s, FFT-NS-2 2.5 s,
FFT-NS-i 5.9 s, spoa 6.2 s, unit-aware 6.6 s, Cactus abPOA 14.4 s; hotspot medians abPOA 25 s, spoa 91 s,
unit-aware 132 s, L-INS-i 463 s (on 4 regions).

### 3.2 The repeat-unit-aware aligner, exactly (`tools/realign_units.py`, `tools/unit_dp.c`, `tools/UNIT_ALIGN.md`)

1. **Repeat model from CHM13 + region.json only.** Period candidates: annotated period, motif length,
   CHM13 tandem-scan periods, and their divisors; scored by **lag identity** (share of i with s[i] =
   s[i+p]) on the CHM13 core; the smallest within 0.02 of the best wins. Unit = best-agreeing p-window,
   refined to the consensus of a wrap-around alignment of the CHM13 core (1-3 rounds), rotated so CHM13's
   array starts at unit position 0. Rejected (-> fallback) when non-TR, < 3 units, identity < 0.60, or
   array < 50% of the core.
2. **Segmentation DP (C):** each distinct N-masked sequence is aligned to the model *free prefix ->
   last 300 bp of CHM13 left flank -> unit as a cycle (any number of times, partial last unit, bypassable)
   -> first 300 bp of right flank -> free suffix*, gap-affine (match +2, mismatch -4, open -6, extend -1).
   Traceback cuts the array at every wrap, so every base gets a unit position: all haplotypes are phased
   by the same unit.
3. **Symbols:** each distinct unit string is a symbol (no merging); similarity enters via the cost.
4. **Unit-level MSA: exact sum-of-pairs progressive aligner (C).** Substituting unit a for b costs
   ed(a,b), gapping a unit costs its length (linear gaps, so profile SP is exact). UPGMA guide tree on
   exact unit-level costs (hashed 3-mer L1 above 150 strings), optimal profile-profile DP per merge, then
   refinement (every tree split, then every single row realigned to the rest, kept if the cross cost
   falls; up to 3 rounds). `--unit-go` adds a gap-open per unit-gap run (0 shipped; 5 gives fewer SV pieces
   at VNTRs but hurts STRs). mafft `--anysymbol` / `--textmatrix` were tried and lose (1.18 vs 1.064),
   because mafft normalises user matrices.
5. **Expansion to bases:** each unit column expanded through the unit positions; same-slot insertions
   star-aligned; flanks aligned separately (mafft L-INS-i; FFT-NS-i > 100 pieces; FFT-NS-2 > 20 kb).
   Nearly lossless (+0.4% at L012184).
6. **Base-level polish:** each row realigned against all others under the exact base-level SP cost +
   gap-open 0.5 per pair, banded max(100, 3p), 2 rounds; fixes boundary shifts (TR499276 6.1x -> 1.00x).
7. **Adequacy guard:** MSA SP cost / sum of pairwise optima (truth-free); if > 1.2 the fallback
   (mafft L-INS-i, or abPOA with `--fallback-engine abpoa`) is also run and the cheaper kept.

Limitations stated: one unit per region (compound arrays modelled badly), unit order differences cannot
be stacked acyclically (L000034 up to 15% above optima), linear objective scatters unit gaps (more SV
pieces; worst truvari representation, raw F1 at matched controls 0.958 -> 0.323), identical unit
variants in different columns (k-mer redundancy 0.518, worst of the good methods), mafft fallback breaks
non-repeat SVs.

### 3.3 Stages 3-4 and pgrealign: what won and why (`results/stage3_final.md`, `results/stage4_chr20.md`, `results/stage4_chr6.md`, `results/full_graph_chr20.md`)

- **Stage 3, local vg call, 149 regions, 3 replicates.** Summed called-haplotype ED at hotspot VNTRs: MC
  120,918; unit_aware 80,048 (21 better / 10 worse); unit_aware__all 82,638; abPOA 85,029; L-INS-i 139,210
  (on its regions); truth graph 15,840. **No method passes the controls**: L-INS-i / unit-aware's mafft
  fallback break the 3 non-repeat SVs (L006556, L002292, L003883); POAs and projections break matched
  controls.
- **Stage 4, whole chr20, patch + re-map + re-call.** abPOA on hap32 over all 624 eligible repeats: SV F1
  0.5325 -> **0.5835** (FP 445 -> 302); held-out chr6: 0.5894 -> 0.6425 (FP 697 -> 453); refined +0.024 /
  +0.041. **Full-panel abPOA projected *hurts*** (0.4798; FN +65; small-variant FPs inside patched
  regions x3); the unit-aware aligner on hap32 or full panel is net negative (0.5128 / 0.4756), despite
  the best Stage-0 cost/opt.
- **Why the full panel hurt (4e, 4m):** a hap32 row aligned inside a 460-way POA matches each unit
  through whichever panel allele fits best, so its deletions break at other alleles' junctions (99.7% of
  hap32 gap runs coincide with a run of some non-hap32 allele); projection keeps the junctions, drops the
  alleles. Not a gap-placement artefact (equal-cost placement explains <= 4.4%); per-row or joint
  normalisation cannot fix it. More cut points = more vg sites = more recombinant (stitched) called
  haplotypes, which carry 71-92% of the call error (4i, 4z).
- **Variants tried on a 22-region test set / whole chr20:** abPOA `-G`, `-X`, `-O`, unbanded, `-J`,
  `-R`, frequency-first, seeded subsets, node merging (all <= 55% of the k-mer duplication removed,
  calls barely move); coarsened sites (much worse); linkage weight 8 (helps repeats, costs small
  variants genome-wide); FAMSA (best full-panel cost/opt 1.066 but fragments deletions; no DNA mode),
  Kalign, mafft FFT-NS-2, MUSCLE5 (too slow); **MST threading** (`mst`, each allele aligned to its
  nearest tree neighbour: fewest sites, best recall, but cost/opt 1.24 -> precision loss), `mst3`,
  **k-medoid abPOA backbone + threading** (`bbt*`: K is a dial between MST and full POA, never both
  good); centrolign pairs in the star (default settings write VNTRs as one DEL + one INS, 2.7-8.4x the
  optimum; `constraint_method 0` 1.11-1.24x pairwise but loses to abPOA in the star, 40x slower).
- **Centre-star wins (4u-4x).** Aligning every distinct full-panel allele pairwise to one centre
  (abPOA unbanded, default convex gaps), left-normalising indels, merging on the centre's columns and
  aligning same-slot insertions to each other with abPOA:
  - `st_chm13`: refined SV F1 **0.7034** (vs abPOA-on-hap32 0.6654, MC 0.6397), significant per region;
    first sample-independent full-panel arm to beat hap32 abPOA.
  - Pairwise consistency F1 vs direct optimal alignments (46,508 pairs, 630 regions, 4w): CHM13 / medoid
    star 0.872, abPOA hap32 0.839, MC 0.825, abPOA full panel 0.784. **The gain is all in pairs of
    different length**: 50-999 bp 0.836 vs 0.792; >= 1 kb 0.671 vs 0.524 (stacked left-normalised
    deletions keep each large indel in one piece). Stars lose a little at equal length (0.979 vs 0.994):
    two non-centre alleles meet only through the centre.
  - Cost: stars have poor unit cost/opt (st_chm13 1.259 on the test set), which is why they had been
    rejected earlier: insertions in different centre slots are never aligned to each other.
- **Why the medoid, not CHM13 (4x).** The CHM13-star comparison was flattered (the scored pairs include
  CHM13, which the star reproduces by construction); without them CHM13 star F1 0.844 < medoid 0.855. On
  24,034 HG002-free random panel pairs split by distance from CHM13, the CHM13 star degrades for alleles
  unlike CHM13 (farthest fifth 0.779 vs medoid 0.801): **reference bias**. The medoid (least
  haplotype-weighted multiset-15-mer Jaccard distance) is a real panel allele, follows the panel's
  composition, is the same for every sample, and calls as well as CHM13 within noise (refined 0.6957).
  An alignment-based medoid moved the centre in 60/636 regions and changed nothing (4y).
- **Whole-contig pgrealign build (chr20, then chr6 held out).** Every tandem repeat region (14,511 on
  chr20) realigned by the medoid star in the full 458-haplotype graph: chr20 short-read refined SV F1
  0.6342 -> 0.6741, ONT 0.6771 -> 0.7611; after adopting **uncapped medoid + scores x4 (long-gap
  extension 1)**: chr6 held out short 0.6968 -> **0.7905**, ONT 0.8074 -> **0.8699**, small variants up
  too; qc pair F1 0.9871 -> 0.9932 over 39,024 chr6 regions.
- **Repeat sites in vg call** (parked, `docs/repeat-sites.md`): genotyping each region as one site of
  whole panel walks cuts stitched haplotypes 182 -> 29 and haplotype ED by ~18-21%; the remaining error is
  mostly copy number / alleles absent from the panel, and MAPQ asymmetry in vg's read term.

### 3.4 The hard-set strategy benchmark (`hardset/README.md`, `results/hardset_chr20*.tsv`)

45 chr20 regions with pair F1 < 0.8; each strategy in its own process, 9 GB / 1,800 s caps.

| strategy | aligned | pair F1 median | LQ | better / worse vs star | cost ratio | time (s) |
|---|---|---|---|---|---|---|
| original graph | 45 | 0.647 | 0.509 | | | |
| **star** (medoid star, current) | 45 | 0.703 | 0.628 | - | 1.34 | 1,773 |
| star_longest | 44 | 0.699 | 0.580 | 26 / 18 | 1.72 | 1,807 |
| star_gapnorm | 45 | 0.717 | 0.642 | 26 / 11 | 2.08 | +282 |
| select (best of 3 stars on other pairs) | 45 | 0.720 | 0.633 | 24 / 2 | 1.50 | 3,862 |
| abpoa_mc (progressive, MC BAR) | 42 | 0.621 | 0.487 | 3 / 36 | 1.42 | 2,382 |
| abpoa_prog (progressive, pgrealign scores) | 41 | 0.621 | 0.501 | 2 / 36 | 1.41 | 2,327 |
| famsa | 44 | 0.667 | 0.506 | 10 / 33 | 1.37 | 6,396 |
| mafft FFT-NS-2 | 43 | 0.613 | 0.486 | 10 / 33 | 1.50 | 8,502 |
| cluster_par | 45 | 0.598 | 0.478 | 4 / 29 | 2.21 | 1,525 |
| cluster_merge | 44 | 0.700 | 0.566 | 12 / 16 | 1.36 | 3,648 |
| units (unit-aware) | 44 | 0.630 | 0.467 | 6 / 35 | 1.43 | 3,207 |
| star_centrolign_c0 | 42 | 0.659 | 0.549 | 8 / 32 | 1.56 | 21,604 |

"No aligner rescues the hard set": best-of-ten per region gives only 0.750. Progressive aligners build
~1/3 of star's sequence but lose 0.05-0.09 F1. Failures concentrate in the 4 regions > 27 kb (abPOA OOM
6x). No strategy repairs the 8 phase-flip regions.

### 3.5 Remaining failure classes

From `docs/results.md` section 5, `docs/sv-fixes-plan.md`, `results/full_graph_chr20.md`:
- **Very long regions / outlier alleles** (memory/time): L012272 (55 kb span, 115 kb full-panel allele),
  L015347, L001909, L007040, L011138, L012184 (one 39.6 kb allele vs CHM13 4.8 kb). spoa/abPOA memory is
  quadratic in the longest allele; all-pairs mafft predicted at 244 GB. pgrealign bands alignments over
  4.5e8 cells and still leaves `too_big` regions (19 on chr20 before banding was extended).
- **Compound / differently ordered arrays** (several motifs, unit-type order differs between haplotypes):
  L011463 (best 1.241), L009658 (1.175), L004145 (1.142), L000034 (1.127). 10/64 hotspots stay above 1.1
  with every method on hap32, 21/64 projected. Cannot be stacked acyclically.
- **Non-repeat SVs:** L006556, L002292, L003883 broken by L-INS-i/G-INS-i/unit-aware's fallback (POAs and
  E-INS-i reproduce MC exactly). In the star, non-homologous or inverted alleles get threaded through the
  centre by chance matches: hence pgrealign's `--min-share 0.5` homology guard (15-mer share), the x4
  scores (long gaps cheap), and the proposed forced-homology guard.
- **Short-period controls:** L011224 (10 bp CA repeat) misplaced by abPOA/spoa.
- **Star expansion:** cross-slot scatter (89.6% of the extra 3.3 Mb on chr20): expansions in different
  centre slots are never merged (TR_chr20_1893502: largest blocks in 187 slots); TR_chr20_64970081 grew
  12.6 kb -> 145 kb and is a real regression (call 616 bp vs truth 4,669 bp).
- **Representation**: correct haplotype, records cut differently from the truth VCF (46% of ONT errors,
  16% of short-read errors after realignment).
- **Not fixable by realignment:** HG002 allele > 1% from every panel allele, or an unsampled haplotype
  (about half of remaining real errors); vg call stitching across sites of one repeat.

---

## 4. pgrealign

### 4.1 Architecture (`README.md`, `docs/method.md`, `docs/reference.md`)

Pure Python stdlib + vg + abPOA (+ a C compiler for the catalogue finders). One module per stage:

| stage | does |
|---|---|
| `catalog` | CHM13 TR catalogue from RepeatMasker + `trfind` (own seed-based TR finder, `c/trfind.c`), optional TRF or `strf` (own TRF reimplementation, `c/strf.c`); classes STR (<= 6 bp), VNTR (>= 7), SAT, LC; shipped `catalogs/chm13v2.0_tr.bed.gz` |
| `prepare` | `vg chunk --contig` the contig component from the whole-genome GBZ, drop `gref_CHM13`, `vg snarls -T -P`, reference node index |
| `regions` | anchor each target (padded `--pad`) on two snarl boundary nodes **of the same chain** that the reference visits once, forward; the narrowest pair within `--max-span` wins; overlapping regions merge; else `too_long` |
| `extract` | one streaming pass over GFA W lines: per region, **runs** (spanning fwd/rev, prefix, suffix, internal, complex), oriented to the region's forward strand; distinct spanning sequences = **alleles** (`a1..`, weight = number of runs), partial ones = **fragments** (`f1..`, side left/right/internal); statuses `complex` (U-turns, nodes visited both ways = inversions), `no_runs`, `invariant`, `ok`. Writes `packages.jsonl.gz`, `occurrences.tsv.gz`, `interiors.tsv.gz` |
| `realign` | medoid star per `ok` package (`star.py`), parallel within a memory budget, resumable; writes `msas.jsonl.gz` |
| `replace` | induce a graph from each MSA (`induce.py`: one node per distinct base per column, unchop, topological numbering, chop at 1,024 bp), remove the interior nodes, add the induced graph + anchor edges, rewrite every recorded run (reversed for `-` runs); dense id renumbering by default; checks every rewritten run spells its member |
| `verify` | every path's sequence (upper-cased md5) and metadata identical to before |
| `qc` | pair/ref precision, recall, F1 vs optimal pairwise alignments, growth, forced homology, before (`--graph`) and after |
| `union` / `project` | realign a recipient graph (e.g. hap32) with a panel's alleles (panel weights choose the centre; recipient-only alleles weight 0), then keep only the recipient's rows |
| `merge` | join per-contig GBZs (non-overlapping id ranges) |
| `run` | prepare -> qc for one contig, resumable |

### 4.2 The medoid star (`pgrealign/star.py`, `docs/method.md` "Realign")

1. **Centre:** allele of positive weight minimising sum_j w_j (1 - J(i,j)), J = Jaccard of 15-mer
   **multisets** (copy number counts); exact up to 1e9 steps, else 1,000-hash bottom sketches; uncapped
   (the earlier 64-comparator cap made centres short and the star expand).
2. **Homology guard:** any allele >= 100 bp sharing < 0.5 of its distinct 15-mers with the centre ->
   region `not_homologous`, left as is (catches non-homologous SVs and inversions).
3. **Pairwise alignment to the centre:** abPOA `-m 0 -r 1`, global, **unbanded** unless centre x allele >
   `--max-cells` (450,000,000; then abPOA's adaptive band). Scores **M 8, X 16, O 16,96, E 8,1** (abPOA
   defaults x4 except long-gap extension 1), so long gaps cost 1/bp vs a match's 8 and the long piece
   wins from 12 bp: a replacement becomes one DEL + one INS, not many chance-match bubbles.
4. **Left normalisation** of every pairwise alignment (indels shifted left while score unchanged), so
   equivalent repeat-unit indels stack in the same centre columns.
5. **Merge:** n centre columns and n+1 insertion slots; distinct insertions in one slot are aligned to
   each other by abPOA (banded, longest first); members inserting the same sequence share the slot row.
   Insertions in different slots are **never** aligned to each other (the main known weakness).
6. **Fragments:** a fragment found in an allele (left = prefix, right = suffix, internal = leftmost
   occurrence) takes that allele's row; otherwise it is aligned to the centre like an allele.
7. Limits/statuses: `too_big` (banded cells of the longest sequence > max-cells, > 40 Mbp of alleles,
   memory estimate > budget), `timeout` (900 s), `aligner_failed`, `no_alleles`, `not_homologous`,
   `error`. Memory: abPOA reserves 20 B/cell of the whole matrix even banded.

### 4.3 How a region is cut out and patched back

- Cut: anchors are snarl boundary nodes of one chain, so every crossing walk visits both; the interior
  is the set of nodes visited by the region's runs (internal anchor-free walks included). Sequences are
  spelled in the region's forward orientation.
- Patch: the induced graph replaces the interior; edges L -> first node of each allele/left-fragment
  walk, last node -> R; an empty allele gets L -> R. Every run in `occurrences.tsv.gz` (walk key, region,
  type, orient, step range, member id) is replaced by the member's new walk. Anchors stay.

### 4.4 The interface plegma must satisfy

Three possible insertion points, from least to most integration:

1. **As the pairwise aligner inside the star** (`star.star_msa(alleles, fragments, align=...)`). Contract
   (`docs/reference.md` "The aligner contract"): `align(seqs, banded=False)` gets a list of >= 2
   distinct, non-empty, upper-case strings (A/C/G/T/N only from star_msa); calls are either (centre,
   allele/fragment) with the centre first, or the distinct insertions of one slot (>= 2, longest first);
   it must return one gapped row per input, in order, all the same length, `-` gaps, each row spelling
   its input; all-gap columns allowed. `banded` is a hint. Raise `star.TimedOut` past a deadline,
   `star.AlignerFailed` on failure. `qc.region_qc` uses the same contract for its optimal pairs.
2. **As a whole-region MSA** written straight into `msas.jsonl.gz`: one JSON object per region (`id`,
   `status: "ok"`, `rows: {member_id: gapped_row}`; all rows the same length, each spelling its member's
   upper-cased sequence; fragments included as rows). `replace` reads only `id`, `status`, `rows`; any
   other status leaves the region's old subgraph. Inputs come from `packages.jsonl.gz` (fields `id`,
   `status`, `alleles[{id, weight, seq}]`, `fragments[{id, weight, side, seq}]`, `ref`). Use
   `realign.write_msas(path, records)`.
3. **In vntr-eval**: an MSA FASTA with one row per input record (hap32 names, or `u0001..` for the panel
   union), then `msa_graph.py` (section 2.4).

Implications of the contract **[inferred]**: plegma must be global end to end (the anchors are shared
bases), must keep rows exactly spelling inputs (N included), should treat weights as haplotype counts
(the centre choice uses them), and must handle empty alleles and fragments (partial sequences touching
one end, or internal).

---

## 5. Implications for plegma

### 5.1 What an MSA must have to score well

Ranked by how strongly the vntr-eval/pgrealign evidence ties each property to calling, not only to
cost/opt:

1. **Global, exact rows.** End-to-end alignment (shared anchor bases); every row spells its input exactly,
   N runs included (as their own columns; `msa_graph.py` makes them N nodes). Fragments and empty alleles
   must be representable for pgrealign.
2. **Each large indel as one contiguous block, stacked in the same columns across all haplotypes that
   carry it.** This is the medoid star's whole advantage: its consistency gain is entirely in pairs of
   different length (F1 50-999 bp 0.836 vs 0.792; >= 1 kb 0.671 vs 0.524 for abPOA-on-hap32;
   `results/stage4_chr20.md` 4w). Fragmented deletions (gap runs per row) predict extra errors
   monotonically (4t). Unit-cost cost/opt cannot see this; affine cost and gap-run counts can.
3. **Pairwise consistency with direct optimal pairwise alignments,** especially **precision** (few
   off-optimum aligned pairs): extra pairs track refined SV FPs better than missed pairs (4w). Report
   qc-style pair F1 / precision, not just cost/opt.
4. **No junctions that the genotyped rows do not need.** Every cut point becomes a vg site and a place to
   stitch recombinant haplotypes (71-92% of called-haplotype error). Aligning a row through a graph of
   hundreds of alleles cuts its indels at other alleles' junctions (4m); projected full-panel POA added
   4-5k sites on chr20 and lost SV F1.
5. **Identical sequence on shared nodes** (`kmer_frac_extra`): tracks read placement (Spearman 0.75)
   and should approach the controls' ~0.05 (MC 0.60, good methods 0.40-0.52).
6. **Repeat units in phase, but gapped as whole contiguous blocks.** The unit-aware aligner stacked units
   in phase (`phase_check` 1.0) and got the best cost/opt (1.048), yet called worse on whole chr20, because
   its linear unit-gap objective scatters unit deletions (more SV pieces, records unlike the truth VCF's)
   and puts identical unit variants in different columns. In-phase units are necessary, not sufficient
   [inferred from 3.2/3.3].
7. **No forced homology**: unrelated or inverted alleles must not be threaded through each other by chance
   matches (pgrealign: `--min-share 0.5` guard; scores x4 so long gaps cost 1/bp vs match 8; candidate
   40-column > 30% mismatch window guard). Non-repeat SV controls must stay exactly 1.000.
8. **Sample- and reference-independence**: a CHM13-centred alignment is biased against alleles unlike
   CHM13 (farthest fifth 0.779 vs medoid 0.801); the panel medoid is the adopted centre.
9. **Bounded growth**: the star's main defect is cross-slot scatter (89.6% of 3.3 Mb excess on chr20;
   TR_chr20_64970081 12.6 kb -> 145 kb, a real calling regression). A POA that merges expansions placed in
   different centre slots would attack exactly this [inferred].

Stage-0 targets to quote: hotspot VNTR cost/opt <= 1.1 (MC 1.566, best 1.048 unit-aware, 1.083 abPOA,
1.067 spoa); matched controls not above MC (1.102) region by region; non-repeat SVs exactly 1.000;
`kmer_frac_extra` below abPOA-hap32's 0.41-0.45 at hotspots. Whole-contig targets: pgrealign medoid
star (qc pair F1 0.9932 on chr6; hard-set pair F1 median 0.703, `select` 0.720, best-of-ten 0.750).

### 5.2 Sizes and runtimes to beat

- **Input sizes:** hap32 hotspot median 30 distinct sequences of median 4.7 kb (max allele median 8.5 kb,
  max 55 kb); full panel median 375 distinct (up to 450), longest up to 115 kb, distinct bp per region up
  to 21.6 Mb (L012272). Whole chr20: 14,511 regions, 13.1 Mb replaced; chr6: 39,026 regions.
- **Per-region aligner time** (shared, loaded 10-core Mac): hap32 abPOA 0.1 s median (1.2 s hotspot,
  max 10 s), spoa 0.9 s, unit-aware 2.2 s (max 421 s), L-INS-i 4.8 s (hotspot 47 s, max 737 s). Full panel:
  abPOA 1 s (hotspot 25 s, max 269 s, but 11 regions refused for memory), unit-aware 13 s (hotspot 132 s,
  max 1,722 s, 9.4 GB; the only method covering 149/149), Cactus-settings abPOA 21 s.
- **pgrealign medoid star**, whole chr20 (458 haplotypes): 2,137 s wall at 8 jobs, 12,702 CPU s, peak 8.5
  GB for 14,511 regions (19 `too_big` before banding). Hard set: 1,773 s for 45 regions, <= 601 s each.
- **Memory walls:** abPOA ~20 B per DP cell (reserved even when banded), spoa 80 B x L^2; pgrealign bands
  above 4.5e8 cells (~21 kb x 21 kb). L012272 (55 kb hap32, 115 kb full panel) and the 25-103 kb chr20
  regions are the memory tests. **Anchoring that removes the quadratic term is where plegma can win
  outright** (coverage), even at equal accuracy.

### 5.3 The hardest regions

- Size/outliers: L012272, L015347, L011138, L001084, L001909, L007009, L007040, L008976, L012184
  (39.6 kb panel allele vs CHM13 4.8 kb), L000289, L016870, L002013.
- Compound / reordered arrays (above 1.1 with every method): L011463, L009658, L004145, L000034.
- Non-repeat SVs that mafft-style methods break: L006556, L002292, L003883; STR control L011224.
- Representation/truvari shift: L001894 (truth writes a 76 bp DEL + 72 bp INS pair).
- chr20 hard set (needs the big data): TR_chr20_1893502 (pair F1 0.19, 452 alleles, 18 kb),
  TR_chr20_64970081 (419 alleles; star blow-up), TR_chr20_65952628 (34 kb allele), TR756034 (14.4 kb,
  450 alleles, dominates chr20 error totals), TR_chr20_18259720, TR_chr20_62032077 (forced homology).

### 5.4 The pairwise MUM experiment [computed]

Following the design note (plegma adds sequences one at a time and anchors each **pairwise** against the
current graph, a path or a consensus, with MUMs unique in both, falling back to the largest SMEMs,
recursing), I measured pairwise MUMs on the 149 hap32 packages.

**Method** (`scratchpad/scripts/mumcheck.c`, driver `mum_driver.py`, summary `mum_summary.py`,
`mum_extra.py`; outputs `scratchpad/out/mum*.tsv`). Brute-force-tested on 40 synthetic repeat pairs
(`test_mumcheck.py`, 0 mismatches).
- Pairs over **distinct** hap32 sequences: `ref` = CHM13 vs every other (1,090 pairs at hotspot VNTRs);
  `med` = the hap32 medoid (pgrealign's rule: least haplotype-weighted multiset-15-mer Jaccard distance)
  vs every other; `rnd` = up to 15 random non-CHM13 pairs per region.
- MUM = maximal exact match >= L bp occurring exactly once in each sequence, left-maximal (MUMmer's
  definition; suffix array + LCP). L = 20, 50, 100.
- **Array zone**: the annotated core on CHM13; on other sequences [left flank, length - right flank)
  **[inferred zone]**. Caveat: the annotated core is not always the repeat's full extent (TR831034's 2-bp
  repeat runs into the "flank").
- **Consistency**: a MUM (p, q, l) is on *some* optimal global unit-cost alignment iff
  ED(A[:p], B[:q]) + ED(A[p+l:], B[q+l:]) = ED(A, B) (exact, tie-independent). Unit cost is evaluate.py's
  headline cost. Not checked against truth (truth.fa is held out, and no true MSA exists).
- **Partition loss**: chain the MUMs, force them as matches, align each gap optimally; report pooled
  sum(cost)/sum(opt). Two chainers: (a) **heaviest chain** (MUMmer-style, maximise total anchored bp);
  (b) **gap-penalised chain**: score = sum of lengths - 1 x sum |diagonal shift| (including from the start
  diagonal 0 and to the end diagonal m - n). **Recursive** = chain, then recompute MUMs unique within each
  gap's substrings and recurse until none; report the largest leaf window that still needs DP.
- 6,543 pairs; L012272's 30 CHM13 pairs had > 1.5e9 cells and were left out of the DP numbers. Runtime 2
  min at 7 processes for all 149 regions.

**Results, hotspot VNTRs** (pooled cost/opt of the forced-anchor pairwise alignment; 1.000 = no loss):

| pairs | MUMs in array per pair (median, L=20) | array MUMs on an optimal path (L=20 / 50 / 100) | flank MUMs on an optimal path (L=20) | heaviest chain, top level (L=20 / 50) | heaviest chain, recursive (L=20 / 50) | gap-penalised chain, top level (L=20 / 50 / 100) | gap-penalised, recursive (L=20 / 50 / 100) | largest DP window left, recursive gap-penalised L=20: median / p90 / max |
|---|---|---|---|---|---|---|---|---|
| CHM13 vs each | 18 | 17% / 19% / 21% | 89% | 1.650 / 1.612 | 1.833 / 1.772 | 1.034 / 1.031 / 1.026 | 1.063 / 1.054 / 1.040 | 1,326 / 5,137 / 17,879 bp |
| medoid vs each | 13 | 14% / 16% / 20% | 88% | 1.770 / 1.716 | 1.960 / 1.878 | 1.054 / 1.051 / 1.045 | 1.096 / 1.086 / 1.067 | 839 / 3,502 / 14,418 bp |
| random pairs | 12 | 18% / 20% / 25% | 81% | 1.546 / 1.522 | 1.708 / 1.660 | 1.040 / 1.037 / 1.034 | 1.068 / 1.061 / 1.050 | 1,001 / 3,747 / 18,446 bp |

Other strata (CHM13 pairs, L = 20; heaviest -> gap-penalised recursive): hotspot_other 7 array MUMs per
pair, 25% on an optimal path, 1.534 -> 1.057; matched controls 6 per pair, 62% on an optimal path,
1.230 -> 1.046; correct controls 1.092 -> 1.021; non-repeat SVs 1.001 -> 1.001.

What this says:
- **Unique MUMs inside VNTR arrays are common, not rare.** 99% of hotspot CHM13 pairs have at least one
  >= 20 bp MUM inside the array (median 18; their lengths sum to a median 51% of the array's length). Hotspot
  arrays are degenerate enough (unit variants) that long unique matches exist. Even at L = 100, 93% of
  pairs have one (median 7).
- **But most of them are false homology.** Only 14-21% of array MUMs lie on any optimal alignment, against
  81-89% in the flanks (and 55-62% in the arrays of matched controls). Mean excess forced by one array
  MUM is ~1,400-1,600 edits: they pair unit copies at the wrong offset (a unit variant shared by two different
  copy positions). Longer MUMs are only slightly better (17% -> 21% from 20 to 100 bp).
- **MUMmer-style heaviest chaining is disastrous in arrays**: 1.65 top level and 1.83 recursive, worse
  than MC's own 1.566. Recursing inside partitions makes it worse at every threshold, because MUMs
  unique *within a gap* are even less reliable.
- **A chain that pays for diagonal shifts fixes most of it**: 1.026-1.034 top level and 1.040-1.063
  recursive at hotspots (1.007-1.014 median per pair). The anchors are usable if the chainer charges the
  indel cost implied between anchors; the residual 3-6% pairwise loss is still large next to a 1.1
  multiple-alignment target, because each forced anchor also constrains every later pair through the
  graph [inferred].
- **The windows left for DP stay big.** Even with recursion, the largest unanchored window per pair is a
  median 0.8-1.3 kb, p90 3.5-5 kb, max 14-18 kb at hotspots (it grows with L: median 2.2 kb at L = 100).
  The top-level chained anchors cover only a median 11-15% of the array (pairs with no chained array anchor: 7-13%).
  So "small POA windows inside VNTRs" will not come from MUMs alone; plan for kb-scale banded or
  bit-parallel DP inside arrays, or anchors that are not exact-unique (unit-phase anchors).
- **Medoid as the anchor partner gives smaller windows** (median 0.8 vs 1.3 kb) but slightly higher loss
  than CHM13 (1.096 vs 1.063 recursive). The medoid is closer to most alleles, so more of them chain, but
  the medoid's own array is a less clean template [inferred].
- **Flanks are reliable anchors** (81-89% on an optimal path at L = 20; 92-97% at L = 50-100). Combined with
  the package design (identical anchor bases), the flank-to-array transition is where anchoring is safe.

### 5.5 Addition order: what the evidence says

- **Progressive POA on hap32: longest first** was best in 5/6 region-tool pairs (abPOA, spoa); Prim
  nearest-first from the medoid ("guide") was close (wins 1/6, by 0.003); CHM13-first ("given") was
  worst; random orders spread cost/opt by 0.11-0.33 and can be worse than MC (`tools/POA.md`).
- **Full panel:** longest first still best in 4/6, despite outliers (45/149 regions have a full-panel
  longest allele > 1.5x the hap32 longest); frequency-first, shortest-first, seeded subsets and guide-tree
  `-p` all worse (`results/stage4_chr20.md` 4e, 4h).
- **But the winning whole-contig method is not progressive at all:** the medoid star. Medoid beat CHM13
  (bias), the longest allele (pair F1 -0.021 on expanded regions; ties abPOA in calling), the abPOA
  consensus and the majority consensus (4u-4x). One-parent threading in Prim order from CHM13 (`mst`)
  removes junctions but costs cross-branch optimality (cost/opt 1.24); a k-medoid abPOA backbone
  restores optimality but brings the junctions back (4q-4s).
- **Suggested for plegma [inferred]:** start the graph from the **medoid** (sample-independent, closest
  to the panel), add the rest **nearest-to-the-graph first** (Prim on the multiset 15-mer distance, as
  pgrealign's medoid metric) and compare against **longest first**; anchor each new sequence against
  the path it is closest to (star-like anchoring keeps indels stacked), but do the in-window DP against
  the graph so that expansions in different "slots" can merge (the star's main defect). Always
  left-normalise or otherwise canonicalise indel placement in repeats.

### 5.6 Experiments to run first

1. **Baseline harness.** Wire plegma's MSA into `msa_graph.py` + `evaluate.py --no-truvari` on all 149
   regions (hap32 arm), and also compute qc-style pair F1 / precision against direct pairwise optima
   (pgrealign `qc.region_qc` with an `Abpoa` aligner, or the C edit distance here). Report per stratum
   against MC, abPOA-longest-first and the medoid star (`pgrealign.star.star_msa` with `Abpoa`,
   runnable on the packages' distinct sequences with hap32 counts as weights).
2. **Anchor chaining ablation** on the pairwise level (this report's harness): heaviest vs
   gap-penalised vs a proper affine chaining cost (minimap2-like), vs flank-only anchors, at L = 20-100;
   score with exact forced-anchor cost and with affine cost. Gate: recursive loss <= 1.01 at hotspots.
3. **Anchor against the graph vs a path vs the consensus**: re-run (2) with the partner being (a) the
   medoid, (b) the nearest already-added allele, (c) the current consensus; measure loss and window
   sizes as the panel grows (34 -> 375 sequences).
4. **SMEM fallback**: when a gap has no MUM, how often does the largest SMEM land at the right unit
   offset? Same exact consistency test (prefix + suffix ED). Expect the same false-homology problem as
   array MUMs; test with diagonal-consistent selection.
5. **Unit-phase anchors** for arrays: since exact-unique anchors cover only ~11-15% of an array, test
   anchors from a unit model (CHM13 period/unit as in `realign_units.py`) or from k-mers that are unique
   *after* removing unit periodicity; or no in-array anchors plus a banded/bit-parallel DP.
6. **Order ablation**: medoid-first + Prim vs longest-first vs CHM13-first vs random (3 seeds), on hap32
   and full panel; measure cost/opt, pair F1, gap runs per row and sites (top-level bubbles).
7. **Fragmentation check on the full panel**: align the union (`panel.py union`), project
   (`panel.py project`), and compare gap runs per hap32 row and sites with the hap32 arm; the full-panel
   arm is where every POA lost so far.
8. **Size frontier**: L012272, L015347, L001084, L012184 (full panel), and memory at 115 kb.
9. **Hard cases**: the compound arrays (L011463, L009658, L004145, L000034) and the non-repeat SV controls
   (must stay at 1.000).
10. **Then** hand pgrealign an `msas.jsonl.gz` (`id`, `status`, `rows`) for chr20 and run its
    replace/verify/qc chain, the stage that decides calling.

### 5.7 Caveats on this report

- Section 5.4 is pairwise, under unit edit cost, over distinct hap32 sequences; MSA-level effects
  (consistency across many pairs, gap placement for calling) are not measured. The gap-penalised chainer
  is a crude 1 x |diagonal shift| rule.
- The array zone of non-CHM13 sequences is inferred from CHM13's flank lengths.
- Everything about calling (Stages 3-4) is quoted from the repositories; nothing was re-run.
- The hard-set alleles and the chr20/chr6 graphs are not in the clone.
