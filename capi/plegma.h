/* plegma's C API: the consensus of a read set (a local assembly's, between two anchors), as `plegma poa -r 0`
 * writes it, from a C or C++ caller that assembles many. Built as build/libplegma.so (`pixi run build-capi`); link
 * with -lplegma and an rpath to it and to the Mojo runtime libraries beside it (see README).
 *
 * Threads: every call is re-entrant, and a model is read-only once made, so threads may share one.
 * Errors: nothing aborts the caller. A model that fails to load is NULL, its message in the caller's buffer; a
 * consensus that fails, or whose arguments are wrong, is a result holding the message (plegma_result_error).
 * Compatibility: plegma_options starts with its own size, so a caller built against an older plegma.h keeps working
 * with a newer library (the fields it does not know take their defaults).
 */
#ifndef PLEGMA_H
#define PLEGMA_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PLEGMA_VERSION_MAJOR 0
#define PLEGMA_VERSION_MINOR 3
#define PLEGMA_VERSION_PATCH 0

/* The library's version, "major.minor.patch" (the header's, PLEGMA_VERSION_*, is what the caller was built with). */
const char* plegma_version(void);

/* Where a read is fixed in a local assembly: on both anchors (spanning it), on the left or the right anchor alone
 * (starting there, ending anywhere; or ending there), or on neither (inside the region: it makes no consensus, and
 * plegma_clustered_consensus aligns it inside one). */
enum { PLEGMA_ANCHOR_BOTH = 0, PLEGMA_ANCHOR_LEFT = 1, PLEGMA_ANCHOR_RIGHT = 2, PLEGMA_ANCHOR_INSIDE = 3 };

/* A read's strand for the homopolymer model: 0 as the consensus runs, 1 reverse-complemented to it, or unknown
 * (both strands' probabilities weigh half each). */
enum { PLEGMA_STRAND_FORWARD = 0, PLEGMA_STRAND_REVERSE = 1, PLEGMA_STRAND_UNKNOWN = 2 };

/* How plegma_consensus runs: plegma poa's options. Set them with plegma_options_init first (plegma poa's defaults,
 * and `size`), then change what you need.
 *
 * - polish_rounds, likely_rounds: rounds of polishing by edit distance, then by likelihood (-P 50,20; >= 0).
 * - pinned: each read pinned to the consensus so far as minimap2 pins an alignment, the alignment only between the
 *   pins (-S, 1), or aligned to the whole graph (-N, 0). Pins pay on long, divergent sequences (a pangenome's
 *   alleles); on a local assembly's reads of one place they only cost (shasta2's chr21: 2.4 times the time, the same
 *   consensus), so a local assembler turns them off.
 * - min_cover: a one-directional assembly (no read spanning) stops where fewer reads go on (-c 5; >= 1).
 * - match, mismatch, gap_open1, gap_extend1, gap_open2, gap_extend2: abPOA's scoring (-M 2 -X 4 -O 4,24 -E 2,1); a
 *   g-long gap costs min(gap_open1 + g gap_extend1, gap_open2 + g gap_extend2), gap_open2 0: affine.
 * - band_b, band_f, band_max: the adaptive band, band_b + band_f L letters either side for a read of L (-b 10
 *   -f 0.01), at most band_max (-B 64; negative: no cap); band_b negative: no band. */
typedef struct plegma_options {
    uint32_t size; /* sizeof(plegma_options), set by plegma_options_init */
    int32_t polish_rounds;
    int32_t likely_rounds;
    int32_t pinned;
    int32_t min_cover;
    int32_t match;
    int32_t mismatch;
    int32_t gap_open1;
    int32_t gap_extend1;
    int32_t gap_open2;
    int32_t gap_extend2;
    int32_t band_b;
    float band_f;
    int32_t band_max;
    int32_t noise_run; /* plegma_clustered_consensus: clusters apart only by single letters in runs this long are one
                          (3; 0 never) */
    float two_lengths_db; /* clustered and phased consensus, with a model: how many dB better two run lengths must
                             explain the reads than one to keep two (15) */
} plegma_options;

void plegma_options_init(plegma_options* options);

typedef struct plegma_model plegma_model;
typedef struct plegma_result plegma_result;

/* A homopolymer model (shasta2's csv: base,strand,left,right,n,m,count,probability,logProbabilityDb), from a file
 * or from its text in memory. NULL on failure, the message (cut to error_size, NUL-ended) in error when given. */
plegma_model* plegma_model_load(const char* path, char* error, size_t error_size);
plegma_model* plegma_model_parse(const char* text, size_t length, char* error, size_t error_size);
void plegma_model_free(plegma_model* model);

/* The consensus of `count` reads: read i is sequences[i], lengths[i] letters long (A, C, G, T, lower case taken as
 * upper; any other byte is read as N, which matches nothing), fixed as anchors[i] says (NULL: all
 * PLEGMA_ANCHOR_BOTH), on strand strands[i] (NULL: all unknown). With a model every homopolymer run takes its most
 * likely length (plegma poa -H). options NULL: plegma poa's defaults. Never NULL; free it with plegma_result_free.
 * One part, polished; two (the left and the right side, unpolished) when no read spans both anchors and the two
 * sides could not be joined. A read may be empty (a read between adjacent anchors); reads far from the consensus
 * (more than 30% edits, an empty one among real ones) shape the graph but are left out of polishing. */
plegma_result* plegma_consensus(size_t count, const uint8_t* const* sequences, const size_t* lengths,
                                const int32_t* anchors, const int32_t* strands, const plegma_model* model,
                                const plegma_options* options);

/* Several consensus sequences of one read set (haplotypes: abPOA's -d, as a diploid caller asks for): at most
 * max_clusters (>= 1) parts, each the consensus of the reads nearest it, made and polished as plegma_consensus makes
 * one; a cluster stays apart only with at least min_freq (0 to 1) of the reads that make consensus sequences (not
 * PLEGMA_ANCHOR_INSIDE). The reads are clustered by the MSA columns where they disagree, then moved, k-means like, to
 * the consensus each aligns to more cheaply until none moves. Each read's alignment to its cluster's consensus
 * (plegma_result_placement) is at the options' scoring as costs (twice the scores: a mismatch 2 (mismatch + match),
 * a gap letter 2 gap_extend + match over twice the openings), gaps leftmost. Arguments as plegma_consensus's. */
plegma_result* plegma_clustered_consensus(size_t count, const uint8_t* const* sequences, const size_t* lengths,
                                          const int32_t* anchors, const int32_t* strands, const plegma_model* model,
                                          const plegma_options* options, int32_t max_clusters, float min_freq);

/* The two haplotypes' consensus sequences of reads already phased: haplotypes[i] 0 or 1 (-1: neither, placed where
 * it fits best). Each made as plegma_consensus makes one from its reads; where the two differ only by homopolymer
 * lengths and a model is given, each such run decided on both haplotypes' reads (their own lengths only when that
 * explains the reads clearly better than one length for both). Two parts always; members and placements as
 * plegma_clustered_consensus's. */
plegma_result* plegma_phased_consensus(size_t count, const uint8_t* const* sequences, const size_t* lengths,
                                       const int32_t* anchors, const int32_t* strands, const int32_t* haplotypes,
                                       const plegma_model* model, const plegma_options* options);

/* A read's alignment to its cluster's consensus: the cluster (a part), the alignment's cost, and the spans it
 * covers, consensus letters consensus_start ..< consensus_end and read letters read_start ..< read_end. */
typedef struct plegma_placement {
    int64_t cluster;
    int64_t cost;
    int64_t consensus_start;
    int64_t consensus_end;
    int64_t read_start;
    int64_t read_end;
} plegma_placement;

/* The result's parts (0 on failure); part i's letters (valid until the result is freed) and their count; the
 * failure's message, or NULL. */
size_t plegma_result_parts(const plegma_result* result);
const uint8_t* plegma_result_sequence(const plegma_result* result, size_t part, size_t* length);
const char* plegma_result_error(const plegma_result* result);

/* Of a plegma_clustered_consensus result: part i's reads (input indices, ascending; valid until the result is
 * freed) and their count; read r's alignment to its cluster's consensus as a CIGAR (=, X, I a letter the read has
 * and the consensus not, D one the consensus has and the read not; NUL-ended), its placement in *placement when
 * given. 0 reads and NULL for a plegma_consensus result, or a part or read out of range. */
size_t plegma_result_members(const plegma_result* result, size_t part, const int32_t** reads);
const char* plegma_result_cigar(const plegma_result* result, size_t read, plegma_placement* placement);
void plegma_result_free(plegma_result* result);

#ifdef __cplusplus
}
#endif

#endif
