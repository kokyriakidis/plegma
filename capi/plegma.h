/* plegma's C API: the consensus of a read set (a local assembly's, between two anchors), as `plegma poa -r 0`
 * writes it, from a C or C++ caller that assembles many. Built as build/libplegma.so (`pixi run build-capi`); link
 * with -lplegma and an rpath to it and to the Mojo runtime libraries beside it (see README).
 *
 * Threads: every call is re-entrant, and a model is read-only once made, so threads may share one.
 * Errors: nothing aborts the caller. A model that fails to load is NULL, its message in the caller's buffer; a
 * consensus that fails is a result holding the message (plegma_result_error).
 */
#ifndef PLEGMA_H
#define PLEGMA_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Where a read is fixed in a local assembly: on both anchors (spanning it), or on the left or the right anchor
 * alone (starting there, ending anywhere; or ending there). */
enum { PLEGMA_ANCHOR_BOTH = 0, PLEGMA_ANCHOR_LEFT = 1, PLEGMA_ANCHOR_RIGHT = 2 };

/* A read's strand for the homopolymer model: 0 as the consensus runs, 1 reverse-complemented to it, or unknown
 * (both strands' probabilities weigh half each). */
enum { PLEGMA_STRAND_FORWARD = 0, PLEGMA_STRAND_REVERSE = 1, PLEGMA_STRAND_UNKNOWN = 2 };

/* plegma poa's default polishing rounds (-P 50,20). */
enum { PLEGMA_POLISH_ROUNDS = 50, PLEGMA_LIKELY_ROUNDS = 20 };

typedef struct plegma_model plegma_model;
typedef struct plegma_result plegma_result;

/* A homopolymer model (shasta2's csv: base,strand,left,right,n,m,count,probability,logProbabilityDb), from a file
 * or from its text in memory. NULL on failure, the message (cut to error_size, NUL-ended) in error when given. */
plegma_model* plegma_model_load(const char* path, char* error, size_t error_size);
plegma_model* plegma_model_parse(const char* text, size_t length, char* error, size_t error_size);
void plegma_model_free(plegma_model* model);

/* The consensus of `count` reads: read i is sequences[i] (A, C, G, T; lower case taken as upper), lengths[i] long,
 * fixed as anchors[i] says (NULL: all PLEGMA_ANCHOR_BOTH), on strand strands[i] (NULL: all unknown). With a model
 * every homopolymer run takes its most likely length (plegma poa -H). Never NULL; free it with plegma_result_free.
 * One part, polished; two (the left and the right side, unpolished) when no read spans both anchors and the two
 * sides could not be joined. */
plegma_result* plegma_consensus(size_t count, const uint8_t* const* sequences, const size_t* lengths,
                                const int32_t* anchors, const int32_t* strands, const plegma_model* model,
                                int32_t polish_rounds, int32_t likely_rounds);

/* The result's parts (0 on failure); part i's letters (valid until the result is freed) and their count; the
 * failure's message, or NULL. */
size_t plegma_result_parts(const plegma_result* result);
const uint8_t* plegma_result_sequence(const plegma_result* result, size_t part, size_t* length);
const char* plegma_result_error(const plegma_result* result);
void plegma_result_free(plegma_result* result);

#ifdef __cplusplus
}
#endif

#endif
