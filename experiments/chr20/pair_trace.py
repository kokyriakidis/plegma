#!/usr/bin/env python3
"""Why an MSA's allele pairs disagree with their optimal alignments (pgrealign qc's pair F1): every letter the MSA
pairs otherwise than the pair's optimum, put in a class by where the MSA holds it, summed over a qc_loop subset.

usage: python3 -I experiments/chr20/pair_trace.py <qc_loop out dir> <label> [--jobs n] [--top n] [--region id]

For each region of <out dir>/<label>.msas.jsonl.gz and each of qc's random allele pairs (CHM13's allele excluded,
as qc draws them), each letter of either allele whose partner in the MSA's induced alignment (a letter of the other
allele, or none) is not its partner in the optimal alignment (asm5 scores, qc_loop's cached optima) is wrong, and is
put by its MSA column:
  ref, other deletes   on a CHM13 letter's column, the other allele gapped there: the other's deletion against
                       CHM13 placed elsewhere than the pair's optimum puts it
  ref, both letters    on a CHM13 letter's column, both alleles lettered: mismatches and shifts
  insertion, same slot in a column CHM13 gaps (an insertion), the other allele inserting at the same slot: the
                       runs inserted at one place aligned to each other otherwise than the pair's optimum
  insertion, elsewhere in a column CHM13 gaps, the other allele inserting nothing at that slot: the same inserted
                       copies placed at another slot (or another place altogether) by the other allele
Prints the classes' shares, the regions with the most wrong letters and, with --region, that region's worst pairs
spelled out by class and slot.
"""
import collections
import gzip
import json
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qc_loop  # noqa: E402  (pgrealign's qc, its align_pair cached)
from pgrealign import qc  # noqa: E402
from pgrealign.star import Abpoa  # noqa: E402

GAP = "-"
MAX_CELLS = 450_000_000
CLASSES = ("ref, other deletes", "ref, both letters", "insertion, same slot", "insertion, elsewhere")


def partners(ra, rb):
    """For each letter of a and of b, its partner in the pairwise rows (index of the other's letter, or -1)."""
    pa, pb = [], []
    for x, y in zip(ra, rb):
        if x == GAP and y == GAP:
            continue
        if x != GAP and y != GAP:
            pa.append(len(pb))
            pb.append(len(pa) - 1)
        elif x != GAP:
            pa.append(-1)
        else:
            pb.append(-1)
    return pa, pb


def columns(row):
    return [c for c, x in enumerate(row) if x != GAP]


def slots(ref_row):
    """For each column, the CHM13 letters before it (an insertion column's slot)."""
    out, n = [], 0
    for x in ref_row:
        out.append(n)
        if x != GAP:
            n += 1
    return out


def classify(rows, ref, a, b, opt, slot_of, inserted_slots):
    ra, rb = rows[a], rows[b]
    ia, ib = qc.induced(ra, rb)
    msa_a, msa_b = partners(ia, ib)
    opt_a, opt_b = partners(*opt)
    counts = collections.Counter()
    where = collections.Counter()
    for me, other, msa_p, opt_p in ((a, b, msa_a, opt_a), (b, a, msa_b, opt_b)):
        cols = columns(rows[me])
        for i, (m, o) in enumerate(zip(msa_p, opt_p)):
            if m == o:
                continue
            c = cols[i]
            if rows[ref][c] != GAP:
                kind = CLASSES[0] if rows[other][c] == GAP else CLASSES[1]
            else:
                kind = CLASSES[2] if slot_of[c] in inserted_slots[other] else CLASSES[3]
            counts[kind] += 1
            where[(kind, slot_of[c])] += 1
    return counts, where


def region(task):
    pkg, rows, abpoa = task
    seqs = {x["id"]: x["seq"].upper() for x in pkg["alleles"]}
    ref = pkg["ref"]
    slot_of = slots(rows[ref])
    inserted_slots = {k: {slot_of[c] for c, x in enumerate(r) if x != GAP and rows[ref][c] == GAP} for k, r in rows.items()}
    total = collections.Counter()
    pairs = []
    with Abpoa(abpoa, qc_loop.SCORES["asm5"]) as align:
        for kind, a, b in qc.compared_pairs(pkg, rows, 20, 20):
            if kind != "pair":
                continue
            # Banded past qc's max_cells, as qc does (and so its cached optimum).
            opt = qc.align_pair(seqs[a], seqs[b], align, banded=len(seqs[a]) * len(seqs[b]) > MAX_CELLS)
            counts, where = classify(rows, ref, a, b, opt, slot_of, inserted_slots)
            total.update(counts)
            pairs.append((sum(counts.values()), a, b, dict(counts), where.most_common(4)))
    pairs.sort(reverse=True)
    return pkg["id"], total, pairs


def main():
    args = sys.argv[1:]
    jobs = int(qc_loop.option(args, "--jobs", "8"))
    top = int(qc_loop.option(args, "--top", "15"))
    only = qc_loop.option(args, "--region", "")
    out, label = args[:2]
    os.environ["QC_OPTIMA_CACHE"] = qc_loop.CACHE = os.path.abspath(os.path.join(out, "optima"))
    msas = {}
    with gzip.open(os.path.join(out, label + ".msas.jsonl.gz"), "rt") as f:
        for line in f:
            m = json.loads(line)
            if m.get("status") == "ok" and len(m["rows"]) > 1 and (not only or m["id"] == only):
                msas[m["id"]] = m["rows"]
    tasks = []
    abpoa = os.environ.get("ABPOA", os.path.expanduser("~/plegma-bench/msa/bin/abpoa"))
    packages = sorted((x for x in os.listdir(out) if x.startswith("packages.")), key=len)[-1:]
    with gzip.open(os.path.join(out, packages[0]), "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] in msas and len(p["alleles"]) > 1 and p.get("ref") in msas[p["id"]]:
                tasks.append((p, msas[p["id"]], abpoa))
    with Pool(jobs) as pool:
        results = pool.map(region, tasks)
    grand = collections.Counter()
    for _, total, _ in results:
        grand.update(total)
    n = sum(grand.values())
    print("%s: %d regions, %d wrong letters in allele pairs" % (label, len(results), n))
    for k in CLASSES:
        print("  %-22s %8d  %5.1f%%" % (k, grand[k], 100 * grand[k] / max(n, 1)))
    print("regions with the most wrong letters:")
    for rid, total, pairs in sorted(results, key=lambda r: -sum(r[1].values()))[:top]:
        print("  %-20s %7d  %s" % (rid, sum(total.values()), ", ".join("%s %d" % (k, total[k]) for k in CLASSES
                                                                         if total[k])))
    if only:
        for rid, total, pairs in results:
            for wrong, a, b, counts, where in pairs[:6]:
                print("  %s %s-%s %d: %s | at (class, slot): %s" % (rid, a, b, wrong, counts, where))


if __name__ == "__main__":
    main()
