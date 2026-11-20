// Tests of plegma's C/C++ API (pixi run test-capi). With arguments, `test_capi [-H model.csv] [-N] file.fa...` instead
// prints each file's consensus as `plegma poa -r 0 [-H model.csv] file.fa` writes it (one FASTA record per part, the
// reads' anchor= and strand= header words read), to check the library against the command.
#include "plegma.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <thread>

namespace {

int failures = 0;

void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
    failures += ok ? 0 : 1;
}

std::size_t distance(const std::string& a, const std::string& b)
{
    std::vector<std::size_t> row(b.size() + 1);
    for (std::size_t j = 0; j <= b.size(); j++) row[j] = j;
    for (std::size_t i = 1; i <= a.size(); i++) {
        std::size_t diagonal = row[0];
        row[0] = i;
        for (std::size_t j = 1; j <= b.size(); j++) {
            std::size_t up = row[j];
            row[j] = std::min({row[j] + 1, row[j - 1] + 1, diagonal + (a[i - 1] == b[j - 1] ? 0 : 1)});
            diagonal = up;
        }
    }
    return row[b.size()];
}

std::string randomSequence(std::mt19937& rng, std::size_t length)
{
    std::string s(length, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}

std::string mutated(std::mt19937& rng, const std::string& s, double rate)
{
    std::uniform_real_distribution<double> u(0., 1.);
    std::string out;
    for (char c : s) {
        double r = u(rng);
        if (r < rate / 3) {
            out += "ACGT"[rng() % 4];
        } else if (r < 2 * rate / 3) {
            continue;
        } else if (r < rate) {
            out += c;
            out += "ACGT"[rng() % 4];
        } else {
            out += c;
        }
    }
    return out;
}

// shasta2's own test model: true lengths 2 and 3, observed 0 to 4; strand 0 reads a run right 90% of the time,
// strand 1 one letter short.
std::string toyModel()
{
    std::ostringstream csv;
    csv << "base,strand,left,right,n,m,count,probability,logProbabilityDb\n";
    for (char base : std::string("ACGT")) {
        for (int strand = 0; strand < 2; strand++) {
            for (int n = 2; n <= 3; n++) {
                int peak = strand == 0 ? n : n - 1;
                for (int m = 0; m <= 4; m++) {
                    double p = m == peak ? 0.9 : 0.025;
                    csv << base << "," << strand << ",*,*," << n << "," << m << "," << (m == peak ? 100 : 0) << ","
                        << p << "," << 10. * std::log10(p) << "\n";
                }
            }
        }
    }
    return csv.str();
}

struct Case {
    std::string truth;
    std::vector<std::string> sequences;
    std::vector<plegma::Read> reads() const
    {
        std::vector<plegma::Read> out;
        for (std::size_t k = 0; k < sequences.size(); k++) {
            plegma::Anchor anchor = k < 6 ? plegma::Anchor::Both : (k % 2 ? plegma::Anchor::Left : plegma::Anchor::Right);
            out.push_back({sequences[k], anchor, k % 2 ? plegma::Strand::Reverse : plegma::Strand::Forward});
        }
        return out;
    }
};

Case makeCase(std::mt19937& rng, std::size_t length)
{
    Case c;
    c.truth = randomSequence(rng, length);
    for (int k = 0; k < 10; k++) {
        std::string piece = c.truth;
        if (k >= 6 and k % 2) piece = c.truth.substr(0, length / 2 + rng() % (length / 2));  // left: a prefix
        if (k >= 6 and not(k % 2)) piece = c.truth.substr(rng() % (length / 2));  // right: a suffix
        c.sequences.push_back(mutated(rng, piece, 0.05));
    }
    return c;
}

std::vector<plegma::Read> readFasta(const std::string& path, std::vector<std::string>& storage)
{
    std::ifstream in(path);
    std::string line;
    std::vector<std::string> headers;
    storage.clear();
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            headers.push_back(line);
            storage.emplace_back();
        } else if (not storage.empty()) {
            storage.back() += line;
        }
    }
    std::vector<plegma::Read> reads;
    for (std::size_t k = 0; k < storage.size(); k++) {
        plegma::Read read{storage[k]};
        if (headers[k].find(" anchor=left") != std::string::npos) read.anchor = plegma::Anchor::Left;
        if (headers[k].find(" anchor=right") != std::string::npos) read.anchor = plegma::Anchor::Right;
        if (headers[k].find(" strand=0") != std::string::npos) read.strand = plegma::Strand::Forward;
        if (headers[k].find(" strand=1") != std::string::npos) read.strand = plegma::Strand::Reverse;
        reads.push_back(read);
    }
    return reads;
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc > 1) {
        // The library's consensus of each file, as plegma poa -r 0 writes it.
        std::unique_ptr<plegma::Model> model;
        plegma::Options options;
        int first = 1;
        while (first < argc and argv[first][0] == '-') {
            std::string flag = argv[first];
            if (flag == "-H" and first + 1 < argc) {
                model = std::make_unique<plegma::Model>(plegma::Model::load(argv[first + 1]));
                first += 2;
            } else if (flag == "-N") {
                options.pinned(false);
                first += 1;
            } else {
                std::cerr << "test_capi: unknown option " << flag << std::endl;
                return 2;
            }
        }
        for (int a = first; a < argc; a++) {
            std::vector<std::string> storage;
            std::vector<plegma::Read> reads = readFasta(argv[a], storage);
            std::vector<std::string> parts = plegma::consensus(reads, model.get(), options);
            if (parts.size() == 1) {
                std::cout << ">Consensus_sequence\n" << parts[0] << "\n";
            } else {
                std::cout << ">Consensus_sequence_left\n" << parts[0] << "\n";
                std::cout << ">Consensus_sequence_right\n" << parts[1] << "\n";
            }
        }
        return 0;
    }

    std::mt19937 rng(17);

    // Spanning reads alone, then with reads fixed at one anchor.
    {
        std::string truth = randomSequence(rng, 2000);
        std::vector<std::string> storage;
        std::vector<plegma::Read> reads;
        for (int k = 0; k < 10; k++) storage.push_back(mutated(rng, truth, 0.05));
        for (const std::string& s : storage) reads.push_back({s});
        std::vector<std::string> parts = plegma::consensus(reads);
        check(parts.size() == 1 and distance(parts[0], truth) <= 2, "spanning reads give the truth");
    }
    {
        Case c = makeCase(rng, 2000);
        std::vector<std::string> parts = plegma::consensus(c.reads());
        check(parts.size() == 1 and distance(parts[0], c.truth) <= 3, "anchored reads give the truth");
    }

    // Models: a missing file fails with a message; a model's text parses; a model changes a run's length.
    try {
        plegma::Model::load("/nonexistent/model.csv");
        check(false, "a missing model file fails");
    } catch (const std::runtime_error& e) {
        check(std::string(e.what()).size() > 8, std::string("a missing model file fails: ") + e.what());
    }
    try {
        plegma::Model::parse("not,a,model\n");
        check(false, "a malformed model fails");
    } catch (const std::runtime_error& e) {
        check(true, std::string("a malformed model fails: ") + e.what());
    }
    plegma::Model toy = plegma::Model::parse(toyModel());
    {
        // A truth with one run (TTT) and no other: three strand 0 reads report it whole, five strand 1 reads one
        // short, as the toy model says strand 1 reads it.
        std::string truth;
        for (int i = 0; truth.size() < 120; i++) {
            char c = "ACG"[rng() % 3];
            if (not truth.empty() and truth.back() == c) continue;
            truth += c;
            if (truth.size() == 60) truth += "TTT";
        }
        std::string shorter = truth;
        shorter.erase(shorter.find("TTT"), 1);
        std::vector<plegma::Read> reads;
        for (int k = 0; k < 8; k++) {
            reads.push_back({k < 3 ? std::string_view(truth) : std::string_view(shorter), plegma::Anchor::Both,
                             k < 3 ? plegma::Strand::Forward : plegma::Strand::Reverse});
        }
        check(plegma::consensus(reads)[0] == shorter, "without the model the reads' vote holds");
        check(plegma::consensus(reads, &toy)[0] == truth, "with the model the run takes its likely length");
    }

    // Threads sharing one model give what one thread gives.
    {
        std::vector<Case> cases;
        for (int k = 0; k < 32; k++) cases.push_back(makeCase(rng, 800 + 50 * k));
        std::vector<std::string> alone(cases.size());
        for (std::size_t k = 0; k < cases.size(); k++) alone[k] = plegma::consensus(cases[k].reads(), &toy)[0];
        std::vector<std::string> together(cases.size());
        std::vector<std::thread> threads;
        for (int t = 0; t < 8; t++) {
            threads.emplace_back([&, t] {
                for (std::size_t k = t; k < cases.size(); k += 8) {
                    together[k] = plegma::consensus(cases[k].reads(), &toy)[0];
                }
            });
        }
        for (std::thread& thread : threads) thread.join();
        check(alone == together, "8 threads sharing a model agree with one");
    }

    // Options: unpinned gives the same consensus here; no polishing gives the bundle as it is.
    {
        Case c = makeCase(rng, 3000);
        std::vector<std::string> pinned = plegma::consensus(c.reads(), &toy);
        std::vector<std::string> unpinned = plegma::consensus(c.reads(), &toy, plegma::Options().pinned(false));
        check(pinned == unpinned, "unpinned agrees with pinned on reads of one place");
        plegma_options defaults;
        plegma_options_init(&defaults);
        check(defaults.polish_rounds == 50 and defaults.likely_rounds == 20 and defaults.pinned == 1 and
                  defaults.min_cover == 5,
              "plegma_options_init sets plegma poa's defaults");
    }

    // Edge cases: no reads; one read.
    {
        std::vector<plegma::Read> none;
        std::vector<std::string> parts = plegma::consensus(none);
        check(parts.size() == 1 and parts[0].empty(), "no reads, an empty consensus");
        std::string one = randomSequence(rng, 300);
        check(plegma::consensus({{one}})[0] == one, "one read is its own consensus");
    }

    // The version, the options' size and defaults; a caller built against an older, shorter plegma_options.
    {
        std::string want = std::to_string(PLEGMA_VERSION_MAJOR) + "." + std::to_string(PLEGMA_VERSION_MINOR) + "." +
                           std::to_string(PLEGMA_VERSION_PATCH);
        check(plegma::version() == want, "plegma_version agrees with the header: " + plegma::version());
        plegma_options o;
        plegma_options_init(&o);
        check(o.size == sizeof(plegma_options) and o.match == 2 and o.mismatch == 4 and o.gap_open1 == 4 and
                  o.gap_extend1 == 2 and o.gap_open2 == 24 and o.gap_extend2 == 1 and o.band_b == 10 and
                  std::fabs(o.band_f - 0.01f) < 1e-6 and o.band_max == 64 and o.noise_run == 3,
              "plegma_options_init sets the size, scoring and band");
        Case c = makeCase(rng, 1500);
        std::vector<std::string> full = plegma::consensus(c.reads(), nullptr, plegma::Options().pinned(false));
        plegma::Options older = plegma::Options().pinned(false);
        older.values.size = 20;  // through min_cover: the fields after it take their defaults
        older.values.match = 99;
        check(plegma::consensus(c.reads(), nullptr, older) == full, "fields past an older size take their defaults");
    }

    // Wrong arguments fail with a message, never a crash.
    {
        auto fails = [](plegma_result* result) {
            bool failed = plegma_result_error(result) != nullptr and plegma_result_parts(result) == 0;
            plegma_result_free(result);
            return failed;
        };
        std::string read = "ACGTACGTAC";
        const std::uint8_t* letters[] = {reinterpret_cast<const std::uint8_t*>(read.data())};
        std::size_t lengths[] = {read.size()};
        std::int32_t badAnchor[] = {7};
        std::int32_t badStrand[] = {9};
        check(fails(plegma_consensus(1, nullptr, lengths, nullptr, nullptr, nullptr, nullptr)), "no sequences fails");
        check(fails(plegma_consensus(1, letters, nullptr, nullptr, nullptr, nullptr, nullptr)), "no lengths fails");
        check(fails(plegma_consensus(1, letters, lengths, badAnchor, nullptr, nullptr, nullptr)), "a bad anchor fails");
        check(fails(plegma_consensus(1, letters, lengths, nullptr, badStrand, nullptr, nullptr)), "a bad strand fails");
        plegma_options o;
        plegma_options_init(&o);
        o.polish_rounds = -1;
        check(fails(plegma_consensus(1, letters, lengths, nullptr, nullptr, nullptr, &o)), "negative rounds fail");
        plegma_options_init(&o);
        o.size = 0;
        check(fails(plegma_consensus(1, letters, lengths, nullptr, nullptr, nullptr, &o)),
              "options not set by plegma_options_init fail");
        // A count or length no caller means (a size_t of -1, say) fails before anything is allocated for it.
        check(fails(plegma_consensus(SIZE_MAX, letters, lengths, nullptr, nullptr, nullptr, nullptr)),
              "a count of SIZE_MAX fails");
        check(fails(plegma_clustered_consensus(SIZE_MAX, letters, lengths, nullptr, nullptr, nullptr, nullptr, 2, 0.25)),
              "a clustered count of SIZE_MAX fails");
        std::int32_t hap[] = {0};
        check(fails(plegma_phased_consensus(SIZE_MAX, letters, lengths, nullptr, nullptr, hap, nullptr, nullptr)),
              "a phased count of SIZE_MAX fails");
        std::size_t huge[] = {std::size_t(1) << 62};
        check(fails(plegma_consensus(1, letters, huge, nullptr, nullptr, nullptr, nullptr)), "a length of 2^62 fails");
        plegma_options_init(nullptr);
        check(true, "plegma_options_init(NULL) does nothing");
        const std::uint8_t* none[] = {nullptr};
        std::size_t zero[] = {0};
        plegma_result* empty = plegma_consensus(1, none, zero, nullptr, nullptr, nullptr, nullptr);
        check(plegma_result_error(empty) == nullptr and plegma_result_parts(empty) == 1,
              "an empty read may have no letters pointer");
        plegma_result_free(empty);
        check(plegma_result_parts(nullptr) == 0 and plegma_result_error(nullptr) == nullptr,
              "a null result has no parts");
    }

    // An empty read among real ones (a read between adjacent anchors) leaves the consensus as the others make it;
    // bytes other than A, C, G, T are read as N and never reach the consensus.
    {
        std::string truth = randomSequence(rng, 300);
        std::vector<std::string> storage;
        for (int k = 0; k < 3; k++) storage.push_back(mutated(rng, truth, 0.06));
        std::vector<plegma::Read> reads;
        for (const std::string& s : storage) reads.push_back({s});
        plegma::Options options = plegma::Options().pinned(false);
        std::string without = plegma::consensus(reads, nullptr, options)[0];
        reads.push_back({std::string_view()});
        check(plegma::consensus(reads, nullptr, options)[0] == without, "an empty read leaves the consensus alone");
        std::vector<std::string> odd = storage;
        for (std::string& s : odd) s[100] = '-';
        std::vector<plegma::Read> oddReads;
        for (const std::string& s : odd) oddReads.push_back({s});
        std::string got = plegma::consensus(oddReads, nullptr, options)[0];
        check(got.find_first_not_of("ACGTN") == std::string::npos, "other bytes never reach the consensus");
    }

    // Two haplotypes: two clusters, each its haplotype's consensus, every read in its own and its CIGAR spanning it;
    // one haplotype: one cluster.
    {
        std::string one = randomSequence(rng, 700);
        std::string two = one;
        for (std::size_t i = 35; i < two.size(); i += 70) two[i] = two[i] == 'A' ? 'C' : 'A';
        two.erase(300, 4);
        std::vector<std::string> storage;
        for (int k = 0; k < 16; k++) storage.push_back(mutated(rng, k % 2 ? two : one, 0.04));
        storage.push_back(one.substr(200, 250));  // inside the region: no anchor
        std::vector<plegma::Read> reads;
        for (std::size_t k = 0; k < storage.size(); k++) {
            reads.push_back({storage[k], k + 1 == storage.size() ? plegma::Anchor::Inside : plegma::Anchor::Both});
        }
        plegma::Clustering got =
            plegma::clusteredConsensus(reads, nullptr, plegma::Options().pinned(false), 2, 0.2f);
        bool split = got.consensus.size() == 2;
        check(split, "two haplotypes, two clusters");
        if (split) {
            bool truths = true, placed = true;
            for (int c = 0; c < 2; c++) {
                std::size_t a = distance(got.consensus[c], one), b = distance(got.consensus[c], two);
                truths = truths and std::min(a, b) <= 1;
            }
            for (std::size_t r = 0; r + 1 < reads.size(); r++) {
                const plegma::Placement& p = got.placements[r];
                const std::string& truth = r % 2 ? two : one;
                placed = placed and distance(got.consensus[p.cluster], truth) <
                                        distance(got.consensus[1 - p.cluster], truth);
                placed = placed and p.readStart == 0 and p.readEnd == storage[r].size() and not p.cigar.empty();
            }
            const plegma::Placement& inside = got.placements.back();
            check(truths, "each cluster's consensus is its haplotype");
            check(placed, "every read in its haplotype's cluster, aligned whole");
            check(inside.cost == 0 and inside.consensusEnd - inside.consensusStart == 250,
                  "a read inside aligned within its haplotype's consensus");
            check(got.members[0].size() + got.members[1].size() == reads.size(), "members cover every read");
        }
        std::vector<plegma::Read> same;
        for (int k = 0; k < 16; k += 2) same.push_back({storage[k]});
        check(plegma::clusteredConsensus(same, nullptr, plegma::Options().pinned(false), 2, 0.2f).consensus.size() == 1,
              "one haplotype, one cluster");
        bool threw = false;
        try {
            plegma::clusteredConsensus(same, nullptr, plegma::Options(), 0, 0.2f);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw, "max_clusters 0 fails");
        plegma_result* plain = plegma_consensus(0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
        const std::int32_t* ids = nullptr;
        check(plegma_result_members(plain, 0, &ids) == 0 and plegma_result_cigar(plain, 0, nullptr) == nullptr,
              "a plain consensus has no members or placements");
        plegma_result_free(plain);
    }

    std::cout << (failures ? "FAILED " : "all passed ") << failures << std::endl;
    return failures ? 1 : 0;
}
