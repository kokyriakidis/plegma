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

enum class Anchor : std::int32_t { Both = PLEGMA_ANCHOR_BOTH, Left = PLEGMA_ANCHOR_LEFT, Right = PLEGMA_ANCHOR_RIGHT };
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

// plegma poa's options (plegma_options in plegma.h): polishing rounds by edit distance, then by likelihood (-P);
// each read pinned to the consensus so far (-S) or aligned to the whole graph (-N: pinned false, what a local
// assembler wants); where a one-directional assembly stops (-c). Default: plegma poa's.
struct Options {
    plegma_options values;
    Options() { plegma_options_init(&values); }
    Options& polishRounds(std::int32_t rounds) { values.polish_rounds = rounds; return *this; }
    Options& likelyRounds(std::int32_t rounds) { values.likely_rounds = rounds; return *this; }
    Options& pinned(bool pinned) { values.pinned = pinned ? 1 : 0; return *this; }
    Options& minCover(std::int32_t cover) { values.min_cover = cover; return *this; }
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

// The consensus of the reads, as plegma poa -r 0 writes it (with -H when a model is given): one part, or the left and
// the right side when no read spans both anchors and the two could not be joined.
inline std::vector<std::string> consensus(const std::vector<Read>& reads, const Model* model = nullptr,
                                          const Options& options = Options())
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
        plegma_consensus(reads.size(), sequences.data(), lengths.data(), anchors.data(), strands.data(),
                         model ? model->get() : nullptr, &options.values),
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

}  // namespace plegma
