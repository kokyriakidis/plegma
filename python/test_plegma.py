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


if __name__ == "__main__":
    unittest.main()
