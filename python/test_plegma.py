"""Tests of plegma's Python binding (pixi run test-python)."""

import os
import random
import sys
import threading
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import plegma  # noqa: E402


def _random(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutated(rng, s, rate):
    out = []
    for c in s:
        r = rng.random()
        if r < rate / 3:
            out.append(rng.choice("ACGT"))
        elif r < 2 * rate / 3:
            continue
        elif r < rate:
            out.append(c + rng.choice("ACGT"))
        else:
            out.append(c)
    return "".join(out)


def _distance(a, b):
    row = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        diagonal, row[0] = row[0], i
        for j in range(1, len(b) + 1):
            up = row[j]
            row[j] = min(row[j] + 1, row[j - 1] + 1, diagonal + (a[i - 1] != b[j - 1]))
            diagonal = up
    return row[len(b)]


class ConsensusTest(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(7)
        self.truth = _random(self.rng, 800)
        self.reads = [_mutated(self.rng, self.truth, 0.05) for _ in range(8)]

    def test_version(self):
        self.assertRegex(plegma.version(), r"^\d+\.\d+\.\d+$")

    def test_spanning_reads_give_the_truth(self):
        parts = plegma.consensus(self.reads, pinned=False)
        self.assertEqual(len(parts), 1)
        self.assertLessEqual(_distance(parts[0], self.truth), 2)

    def test_bytes_and_str_agree(self):
        self.assertEqual(
            plegma.consensus(self.reads, pinned=False), plegma.consensus([r.encode() for r in self.reads], pinned=False)
        )

    def test_anchored_reads(self):
        reads = self.reads[:4] + [self.reads[4][:500], self.reads[5][-500:]]
        parts = plegma.consensus(reads, anchors=["both"] * 4 + ["left", "right"], strands=[0, 1, None, 0, 1, 0])
        self.assertEqual(len(parts), 1)
        self.assertLessEqual(_distance(parts[0], self.truth), 3)

    def test_options_and_errors(self):
        self.assertEqual(plegma.consensus([]), [""])
        with self.assertRaises(plegma.Error):
            plegma.consensus(self.reads, polish_rounds=-1)
        with self.assertRaises(TypeError):
            plegma.consensus(self.reads, no_such_option=1)
        with self.assertRaises(ValueError):
            plegma.consensus(self.reads, anchors=["both"])
        with self.assertRaises(ValueError):
            plegma.consensus(self.reads[:1], anchors=["sideways"])
        with self.assertRaises(plegma.Error):
            plegma.Model.load("/nonexistent/model.csv")
        # Integers a C int32_t cannot hold fail here, not wrapped into another value.
        with self.assertRaises(ValueError):
            plegma.consensus(self.reads[:1], anchors=[2**32])
        with self.assertRaises(ValueError):
            plegma.consensus(self.reads, polish_rounds=2**32 + 1)
        with self.assertRaises(plegma.Error):
            plegma.consensus(self.reads, two_lengths_db=float("nan"))

    def test_threads_agree(self):
        alone = plegma.consensus(self.reads, pinned=False)
        got = [None] * 4

        def run(k):
            got[k] = plegma.consensus(self.reads, pinned=False)

        threads = [threading.Thread(target=run, args=(k,)) for k in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self.assertTrue(all(g == alone for g in got))



def _other(s):
    """Another allele of `s`: a SNP every 60 letters, a 3-letter deletion at 200."""
    out = []
    for i, c in enumerate(s):
        if 200 <= i < 203:
            continue
        out.append("ACGT"[("ACGT".index(c) + 1) % 4] if i % 60 == 30 else c)
    return "".join(out)


def _cigar_spells(cigar, consensus, read, c0, r0):
    i, j, count = c0, r0, ""
    for ch in cigar:
        if ch.isdigit():
            count += ch
            continue
        for _ in range(int(count)):
            if ch in "=X":
                if (consensus[i] == read[j]) != (ch == "="):
                    return False
                i, j = i + 1, j + 1
            elif ch == "I":
                j += 1
            else:
                i += 1
        count = ""
    return True


class ClusteredTest(unittest.TestCase):
    def setUp(self):
        rng = random.Random(11)
        self.one = _random(rng, 700)
        self.two = _other(self.one)
        self.reads = [_mutated(rng, self.one if k % 2 == 0 else self.two, 0.04) for k in range(16)]

    def test_two_haplotypes(self):
        got = plegma.clustered_consensus(self.reads, max_clusters=2, min_freq=0.2, pinned=False)
        self.assertEqual(len(got.consensus), 2)
        self.assertEqual(sorted(len(m) for m in got.members), [8, 8])
        for k, cons in enumerate(got.consensus):
            self.assertLessEqual(min(_distance(cons, self.one), _distance(cons, self.two)), 1)
        for r, placed in enumerate(got.placements):
            self.assertIn(r, got.members[placed.cluster])
            mine = got.consensus[placed.cluster]
            self.assertTrue(_cigar_spells(placed.cigar, mine, self.reads[r], placed.consensus_start, placed.read_start))

    def test_phased(self):
        haplotypes = [k % 2 for k in range(len(self.reads))]
        got = plegma.phased_consensus(self.reads, haplotypes, pinned=False)
        self.assertEqual(len(got.consensus), 2)
        self.assertLessEqual(_distance(got.consensus[0], self.one), 1)
        self.assertLessEqual(_distance(got.consensus[1], self.two), 1)
        self.assertEqual(got.members[0], list(range(0, 16, 2)))

    def test_inside_and_errors(self):
        reads = self.reads + [self.one[200:450]]
        anchors = ["both"] * len(self.reads) + ["inside"]
        got = plegma.clustered_consensus(reads, anchors=anchors, max_clusters=2, min_freq=0.2, pinned=False)
        self.assertEqual(got.placements[-1].cost, 0)
        with self.assertRaises(plegma.Error):
            plegma.clustered_consensus(self.reads, max_clusters=0)
        with self.assertRaises(ValueError):
            plegma.phased_consensus(self.reads, [0, 1])
        with self.assertRaises(TypeError):
            plegma.clustered_consensus(self.reads, no_such=1)


if __name__ == "__main__":
    unittest.main()
