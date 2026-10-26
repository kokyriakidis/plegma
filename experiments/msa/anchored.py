#!/usr/bin/env python3
"""Local assemblies between two anchors: reads spanning both, reads fixed at one anchor only, and none spanning, as
a local assembler (Shasta's LocalAssembly) hands them to a consensus tool; plegma poa against pericles (theseus-lib's
consensus branch).

usage: python3 -I experiments/msa/anchored.py <work dir> <plegma> <pericles> <karan> [--jobs n] [--only regex]

Each case is a truth of L letters (random, a tandem array of 30 copies of a 40-letter unit, 3% apart, in its middle)
and reads of it at an error rate (substitutions, insertions and deletions alike): `both` reads span it whole, `left`
reads start at its start and end at a random place past 30% of it, `right` reads end at its end and start before 70%
of it. Scenarios: many_both (10 both, 5 left, 5 right), few_both (1, 10, 10), no_both (0, 10, 10), left_only (0, 15,
0: a telomere, one-directional), right_only (0, 0, 15: the other end's) and overhang (5, 10, 10, the one-sided reads
cut as shasta2's LocalAssembly7 cuts them: 0.1 L + 30 letters past the far anchor, into the sequence beyond it). plegma reads each read's anchor from its header (`anchor=both|left|right`); pericles
from its own (`> reversed ends_free weight drop lag`: both 0 0 1 0 0, left 0 1 1 0 0, right 1 1 1 0 0), its first read
spanning and the spanning reads first, as it requires; it cannot run a case with none.

Each consensus is scored by its edit distance to the truth (karan); a one-directional one (left_only) against the
truth's best prefix (`--mode prefix`), right_only against its best suffix (`--mode suffix`), with its length beside the
truth's. Prints per case and tool the errors, the
time and the peak memory, then totals per scenario.
"""
import os
import random
import re
import subprocess
import sys
import zlib
from multiprocessing import Pool

SCENARIOS = {"many_both": (10, 5, 5), "few_both": (1, 10, 10), "no_both": (0, 10, 10), "left_only": (0, 15, 0),
             "right_only": (0, 0, 15), "overhang": (5, 10, 10)}
LENGTHS = [5000, 20000, 50000]
RATES = [0.05, 0.10]
SEEDS = [1, 2]
PLEGMA = PERICLES = KARAN = None


def mutated(rng, s, rate):
    out = []
    for c in s:
        r = rng.random()
        if r < rate / 3:
            out.append(rng.choice("ACGT"))
        elif r < 2 * rate / 3:
            continue
        elif r < rate:
            out.append(c)
            out.append(rng.choice("ACGT"))
        else:
            out.append(c)
    return "".join(out)


def truth_of(rng, length):
    unit = "".join(rng.choice("ACGT") for _ in range(40))
    array = "".join(mutated(rng, unit, 0.03) for _ in range(30))
    flank = (length - len(array)) // 2
    left = "".join(rng.choice("ACGT") for _ in range(flank))
    right = "".join(rng.choice("ACGT") for _ in range(length - len(array) - flank))
    return left + array + right


def make_case(path, scenario, length, rate, seed):
    # A checksum, not hash(): Python salts string hashes per process, and the cases must be the same every run.
    rng = random.Random(zlib.crc32(("%s %d %.2f %d" % (scenario, length, rate, seed)).encode()))
    truth = truth_of(rng, length)
    n_both, n_left, n_right = SCENARIOS[scenario]
    reads = []
    for _ in range(n_both):
        reads.append(("both", mutated(rng, truth, rate)))
    if scenario == "overhang":
        margin = int(0.1 * length) + 30
        before = "".join(rng.choice("ACGT") for _ in range(margin))
        after = "".join(rng.choice("ACGT") for _ in range(margin))
        for _ in range(n_left):
            reads.append(("left", mutated(rng, truth + after, rate)))
        for _ in range(n_right):
            reads.append(("right", mutated(rng, before + truth, rate)))
        n_left = n_right = 0
    for _ in range(n_left):
        end = rng.randint(int(0.3 * length), int(0.95 * length))
        reads.append(("left", mutated(rng, truth[:end], rate)))
    for _ in range(n_right):
        start = rng.randint(int(0.05 * length), int(0.7 * length))
        reads.append(("right", mutated(rng, truth[start:], rate)))
    with open(path + ".truth.fa", "w") as f:
        f.write(">truth\n%s\n" % truth)
    with open(path + ".fa", "w") as f:
        for k, (anchor, seq) in enumerate(reads):
            f.write(">r%d anchor=%s\n%s\n" % (k, anchor, seq))
    if n_both > 0:
        head = {"both": "0 0 1 0 0", "left": "0 1 1 0 0", "right": "1 1 1 0 0"}
        with open(path + ".pericles.fa", "w") as f:
            for anchor, seq in reads:
                f.write("> %s\n%s\n" % (head[anchor], seq))


TIMEOUT = 600
"""Seconds a tool gets per case."""


def timed(cmd, out):
    """Runs `cmd` with stdout to `out`; (seconds, peak KB) or None when it failed or ran past TIMEOUT."""
    stats = out + ".time"
    proc = subprocess.run(["/usr/bin/time", "-o", stats, "-f", "%e %M", "timeout", str(TIMEOUT)] + cmd,
                          stdout=open(out, "w"), stderr=subprocess.DEVNULL)
    if proc.returncode != 0:
        return None
    seconds, kb = open(stats).read().split()[-2:]
    return float(seconds), int(kb)


def consensus_of(path):
    """The first FASTA record's sequence of `path` (pericles writes its consensus as one record)."""
    seq, on = [], False
    for line in open(path):
        if line.startswith(">"):
            if on:
                break
            on = True
            continue
        if on:
            seq.append(line.strip())
    return "".join(seq)


def errors(truth_path, cons, mode, tmp):
    with open(tmp, "w") as f:
        f.write(">c\n%s\n" % cons)
    out = subprocess.run([KARAN, "--distance", "--mode", mode, truth_path, tmp], capture_output=True, text=True).stdout
    # --distance: query, reference, cost.
    last = out.strip().splitlines()[-1].split("\t")
    return int(last[2]) if len(last) > 2 and last[2].isdigit() else None


def run_case(task):
    path, scenario = task
    truth_path = path + ".truth.fa"
    mode = {"left_only": "prefix", "right_only": "suffix"}.get(scenario, "global")
    rows = []
    got = timed([PLEGMA, "poa", "-r", "0", path + ".fa"], path + ".plegma.out")
    if got:
        cons = consensus_of(path + ".plegma.out")
        rows.append(("plegma", errors(truth_path, cons, mode, path + ".plegma.c.fa"), len(cons)) + got)
    else:
        rows.append(("plegma", None, 0, 0.0, 0))
    if os.path.exists(path + ".pericles.fa"):
        for name, kind in (("pericles_hb", "2"), ("pericles_mv", "3")):
            out = path + "." + name + ".out"
            got = timed([PERICLES, "-m", "0", "-x", "2", "-o", "3", "-e", "1", "-t", kind, "-f", out, "-s",
                         path + ".pericles.fa"], out + ".log")
            if got and os.path.exists(out):
                cons = consensus_of(out)
                rows.append((name, errors(truth_path, cons, mode, out + ".c.fa"), len(cons)) + got)
            else:
                rows.append((name, None, 0, 0.0, 0))
    return path, scenario, rows


def main():
    global PLEGMA, PERICLES, KARAN
    args = sys.argv[1:]
    jobs = int(args[args.index("--jobs") + 1]) if "--jobs" in args else 4
    only = args[args.index("--only") + 1] if "--only" in args else ""
    work, PLEGMA, PERICLES, KARAN = args[:4]
    os.makedirs(work, exist_ok=True)
    tasks = []
    for scenario in SCENARIOS:
        for length in LENGTHS:
            for rate in RATES:
                for seed in SEEDS:
                    name = "%s_L%d_r%02d_s%d" % (scenario, length, int(rate * 100), seed)
                    if only and not re.search(only, name):
                        continue
                    path = os.path.join(work, name)
                    if not os.path.exists(path + ".fa"):
                        make_case(path, scenario, length, rate, seed)
                    tasks.append((path, scenario))
    totals = {}
    with Pool(jobs, initializer=init, initargs=(PLEGMA, PERICLES, KARAN)) as pool:
        for path, scenario, rows in pool.imap(run_case, tasks):
            name = os.path.basename(path)
            truth_len = len(consensus_of(path + ".truth.fa"))
            cells = []
            for tool, err, length, seconds, kb in rows:
                cells.append("%s %s err, %d/%d letters, %.2f s, %d MB" % (
                    tool, "FAIL" if err is None else err, length, truth_len, seconds, kb // 1024))
                t = totals.setdefault((scenario, tool), [0, 0, 0.0, 0])
                t[0] += 1
                t[1] += err if err is not None else 0
                t[2] += seconds
                t[3] = max(t[3], kb)
                if err is None:
                    t.append("fail")
            print("%-28s %s" % (name, " | ".join(cells)), flush=True)
    print("\ntotals per scenario")
    for (scenario, tool), t in sorted(totals.items()):
        fails = len(t) - 4
        print("%-10s %-12s cases %2d errors %6d time %8.2f s peak %6d MB%s" % (
            scenario, tool, t[0], t[1], t[2], t[3] // 1024, " failed %d" % fails if fails else ""))


def init(plegma, pericles, karan):
    global PLEGMA, PERICLES, KARAN
    PLEGMA, PERICLES, KARAN = plegma, pericles, karan


if __name__ == "__main__":
    main()
