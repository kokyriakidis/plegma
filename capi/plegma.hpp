// plegma's C++ API: a header-only wrapper over plegma.h. A Model owns a loaded homopolymer model (load it once,
// share it across threads); consensus() takes the reads and returns the consensus parts, throwing
// std::runtime_error on failure.
//
//     plegma::Model model = plegma::Model::load("HomopolymerModel-ONT-R10.4.1-sup5.0.0-2026-09-28.csv");
//     std::vector<plegma::Read> reads = {{"ACGT...", plegma::Anchor::Both, plegma::Strand::Forward}, ...};
//     std::vector<std::string> parts = plegma::consensus(reads, &model, plegma::Options().pinned(false));
#pragma once

#include "plegma.h"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace plegma {

static_assert(sizeof(plegma_options) == 68, "plegma_options laid out as the library reads it");

enum class Anchor : std::int32_t {
    Both = PLEGMA_ANCHOR_BOTH,
    Left = PLEGMA_ANCHOR_LEFT,
    Right = PLEGMA_ANCHOR_RIGHT,
    Inside = PLEGMA_ANCHOR_INSIDE
};
enum class Strand : std::int32_t {
    Forward = PLEGMA_STRAND_FORWARD,
    Reverse = PLEGMA_STRAND_REVERSE,
    Unknown = PLEGMA_STRAND_UNKNOWN
};

// A read of a local assembly: its letters (A, C, G, T), where it is fixed, its strand (for the homopolymer model).
struct Read {
    std::string_view sequence;
    Anchor anchor = Anchor::Both;
    Strand strand = Strand::Unknown;
};

// The library's version ("major.minor.patch").
inline std::string version() { return plegma_version(); }

// plegma poa's options (plegma_options in plegma.h, which says what each means): polishing rounds by edit distance,
// then by likelihood (-P); each read pinned to the consensus so far (-S) or aligned to the whole graph (-N: pinned
// false, what a local assembler wants); where a one-directional assembly stops (-c); the scoring (-M -X -O -E) and
// the adaptive band (-b -f -B). Default: plegma poa's. A value out of range fails the consensus call.
struct Options {
    plegma_options values;
    Options() { plegma_options_init(&values); }
    Options& polishRounds(std::int32_t rounds) { values.polish_rounds = rounds; return *this; }
    Options& likelyRounds(std::int32_t rounds) { values.likely_rounds = rounds; return *this; }
    Options& pinned(bool pinned) { values.pinned = pinned ? 1 : 0; return *this; }
    Options& minCover(std::int32_t cover) { values.min_cover = cover; return *this; }
    Options& scoring(std::int32_t match, std::int32_t mismatch)
    {
        values.match = match;
        values.mismatch = mismatch;
        return *this;
    }
    Options& gaps(std::int32_t open1, std::int32_t extend1, std::int32_t open2 = 0, std::int32_t extend2 = 0)
    {
        values.gap_open1 = open1;
        values.gap_extend1 = extend1;
        values.gap_open2 = open2;
        values.gap_extend2 = extend2;
        return *this;
    }
    Options& noiseRun(std::int32_t run) { values.noise_run = run; return *this; }
    Options& twoLengthsDb(float db) { values.two_lengths_db = db; return *this; }
    Options& threads(std::int32_t count) { values.threads = count; return *this; }
    Options& band(std::int32_t b, float f, std::int32_t max = 64)
    {
        values.band_b = b;
        values.band_f = f;
        values.band_max = max;
        return *this;
    }
};

// A homopolymer model (shasta2's csv), read-only once made.
class Model {
public:
    static Model load(const std::string& path)
    {
        char error[512] = {0};
        plegma_model* model = plegma_model_load(path.c_str(), error, sizeof error);
        if (not model) {
            throw std::runtime_error(std::string("plegma: ") + error);
        }
        return Model(model);
    }

    static Model parse(std::string_view text)
    {
        char error[512] = {0};
        plegma_model* model = plegma_model_parse(text.data(), text.size(), error, sizeof error);
        if (not model) {
            throw std::runtime_error(std::string("plegma: ") + error);
        }
        return Model(model);
    }

    const plegma_model* get() const { return model.get(); }

private:
    struct Free {
        void operator()(plegma_model* model) const { plegma_model_free(model); }
    };
    std::unique_ptr<plegma_model, Free> model;
    explicit Model(plegma_model* model) : model(model) {}
};

// Reusable scratch buffers for consensus calls: create one per worker thread and pass to consensus(),
// clusteredConsensus(), or phasedConsensus() to avoid reallocating buffers across windows.
class Scratch {
public:
    Scratch() : scratch(plegma_scratch_create()) {}
    plegma_scratch* get() const { return scratch.get(); }

private:
    struct Free {
        void operator()(plegma_scratch* s) const { plegma_scratch_free(s); }
    };
    std::unique_ptr<plegma_scratch, Free> scratch;
};

// The consensus of the reads, as plegma poa -r 0 writes it (with -H when a model is given): one part, or the left and
// the right side when no read spans both anchors and the two could not be joined.
inline std::vector<std::string> consensus(const std::vector<Read>& reads, const Model* model = nullptr,
                                          const Options& options = Options(), Scratch* scratch = nullptr)
{
    std::vector<const std::uint8_t*> sequences;
    std::vector<std::size_t> lengths;
    std::vector<std::int32_t> anchors;
    std::vector<std::int32_t> strands;
    sequences.reserve(reads.size());
    lengths.reserve(reads.size());
    anchors.reserve(reads.size());
    strands.reserve(reads.size());
    for (const Read& read : reads) {
        sequences.push_back(reinterpret_cast<const std::uint8_t*>(read.sequence.data()));
        lengths.push_back(read.sequence.size());
        anchors.push_back(static_cast<std::int32_t>(read.anchor));
        strands.push_back(static_cast<std::int32_t>(read.strand));
    }
    std::unique_ptr<plegma_result, void (*)(plegma_result*)> result(
        plegma_consensus_with_scratch(scratch ? scratch->get() : nullptr, reads.size(), sequences.data(),
                                      lengths.data(), anchors.data(), strands.data(), model ? model->get() : nullptr,
                                      &options.values),
        plegma_result_free);
    if (const char* error = plegma_result_error(result.get())) {
        throw std::runtime_error(std::string("plegma: ") + error);
    }
    std::vector<std::string> parts;
    for (std::size_t part = 0; part < plegma_result_parts(result.get()); part++) {
        std::size_t length = 0;
        const std::uint8_t* letters = plegma_result_sequence(result.get(), part, &length);
        parts.emplace_back(reinterpret_cast<const char*>(letters), length);
    }
    return parts;
}

// A read's alignment to its cluster's consensus (plegma_placement): the cluster, the CIGAR (=, X, I, D), its cost and
// the spans it covers, [consensusStart, consensusEnd) and [readStart, readEnd).
struct Placement {
    int cluster = -1;
    std::string cigar;
    std::int64_t cost = 0;
    std::size_t consensusStart = 0, consensusEnd = 0, readStart = 0, readEnd = 0;
};

// plegma_clustered_consensus's outcome: per cluster its consensus and its reads (input indices); per read its
// placement.
struct Clustering {
    std::vector<std::string> consensus;
    std::vector<std::vector<int>> members;
    std::vector<Placement> placements;
};

namespace detail {
inline Clustering clustering(plegma_result* raw, std::size_t count)
{
    std::unique_ptr<plegma_result, void (*)(plegma_result*)> result(raw, plegma_result_free);
    if (const char* error = plegma_result_error(result.get())) {
        throw std::runtime_error(std::string("plegma: ") + error);
    }
    Clustering out;
    for (std::size_t part = 0; part < plegma_result_parts(result.get()); part++) {
        std::size_t length = 0;
        const std::uint8_t* letters = plegma_result_sequence(result.get(), part, &length);
        out.consensus.emplace_back(reinterpret_cast<const char*>(letters), length);
        const std::int32_t* ids = nullptr;
        std::size_t n = plegma_result_members(result.get(), part, &ids);
        out.members.emplace_back(ids, ids + n);
    }
    for (std::size_t r = 0; r < count; r++) {
        plegma_placement held{};
        const char* cigar = plegma_result_cigar(result.get(), r, &held);
        Placement placement;
        if (!cigar) {
            // No placement for this read: the defaults (cluster -1), not the zeroed struct's cluster 0.
            out.placements.push_back(std::move(placement));
            continue;
        }
        placement.cluster = static_cast<int>(held.cluster);
        placement.cigar = cigar;
        placement.cost = held.cost;
        placement.consensusStart = static_cast<std::size_t>(held.consensus_start);
        placement.consensusEnd = static_cast<std::size_t>(held.consensus_end);
        placement.readStart = static_cast<std::size_t>(held.read_start);
        placement.readEnd = static_cast<std::size_t>(held.read_end);
        out.placements.push_back(std::move(placement));
    }
    return out;
}
}  // namespace detail

// At most maxClusters consensus sequences of the reads (haplotypes, abPOA's -d), each from the reads nearest it, a
// cluster kept apart only with at least minFreq of the reads (see plegma_clustered_consensus).
inline Clustering clusteredConsensus(const std::vector<Read>& reads, const Model* model = nullptr,
                                     const Options& options = Options(), int maxClusters = 2, float minFreq = 0.25f,
                                     Scratch* scratch = nullptr)
{
    std::vector<const std::uint8_t*> sequences;
    std::vector<std::size_t> lengths;
    std::vector<std::int32_t> anchors;
    std::vector<std::int32_t> strands;
    for (const Read& read : reads) {
        sequences.push_back(reinterpret_cast<const std::uint8_t*>(read.sequence.data()));
        lengths.push_back(read.sequence.size());
        anchors.push_back(static_cast<std::int32_t>(read.anchor));
        strands.push_back(static_cast<std::int32_t>(read.strand));
    }
    return detail::clustering(
        plegma_clustered_consensus_with_scratch(scratch ? scratch->get() : nullptr, reads.size(), sequences.data(),
                                                lengths.data(), anchors.data(), strands.data(),
                                                model ? model->get() : nullptr, &options.values, maxClusters, minFreq),
        reads.size());
}

// The two haplotypes' consensus sequences of reads already phased (haplotypes[i] 0, 1 or -1; see
// plegma_phased_consensus).
inline Clustering phasedConsensus(const std::vector<Read>& reads, const std::vector<int>& haplotypes,
                                  const Model* model = nullptr, const Options& options = Options(),
                                  Scratch* scratch = nullptr)
{
    if (haplotypes.size() != reads.size()) throw std::runtime_error("plegma: one haplotype per read");
    std::vector<const std::uint8_t*> sequences;
    std::vector<std::size_t> lengths;
    std::vector<std::int32_t> anchors, strands, haps;
    for (std::size_t r = 0; r < reads.size(); r++) {
        sequences.push_back(reinterpret_cast<const std::uint8_t*>(reads[r].sequence.data()));
        lengths.push_back(reads[r].sequence.size());
        anchors.push_back(static_cast<std::int32_t>(reads[r].anchor));
        strands.push_back(static_cast<std::int32_t>(reads[r].strand));
        haps.push_back(static_cast<std::int32_t>(haplotypes[r]));
    }
    return detail::clustering(
        plegma_phased_consensus_with_scratch(scratch ? scratch->get() : nullptr, reads.size(), sequences.data(),
                                             lengths.data(), anchors.data(), strands.data(), haps.data(),
                                             model ? model->get() : nullptr, &options.values),
        reads.size());
}

}  // namespace plegma
