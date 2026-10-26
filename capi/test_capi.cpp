// Tests of plegma's C/C++ API (pixi run test-capi). With arguments, `test_capi [-H model.csv] file.fa...` instead
// prints each file's consensus as `plegma poa -r 0 [-H model.csv] file.fa` writes it (one FASTA record per part, the
// reads' anchor= and strand= header words read), to check the library against the command.
#include "plegma.hpp"

#include <algorithm>
#include <cmath>
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
        int first = 1;
        if (std::string(argv[1]) == "-H" and argc > 2) {
            model = std::make_unique<plegma::Model>(plegma::Model::load(argv[2]));
            first = 3;
        }
        for (int a = first; a < argc; a++) {
            std::vector<std::string> storage;
            std::vector<plegma::Read> reads = readFasta(argv[a], storage);
            std::vector<std::string> parts = plegma::consensus(reads, model.get());
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

    // Edge cases: no reads; one read.
    {
        std::vector<plegma::Read> none;
        std::vector<std::string> parts = plegma::consensus(none);
        check(parts.size() == 1 and parts[0].empty(), "no reads, an empty consensus");
        std::string one = randomSequence(rng, 300);
        check(plegma::consensus({{one}})[0] == one, "one read is its own consensus");
    }

    std::cout << (failures ? "FAILED " : "all passed ") << failures << std::endl;
    return failures ? 1 : 0;
}
