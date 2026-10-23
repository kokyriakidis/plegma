#!/usr/bin/env python3
"""Why an MSA loses pair F1 on a hard-set region: each scored pair's induced alignment against qc's optimum.

usage: python3 -I experiments/hardset/pairs.py <vntr-eval> <work dir> <region> <label> [<label> ...]

For each pair qc scores (cached optima in <work dir>/optima, MSAs in <work dir>/msa/<label>/, as hardset.py writes
them): the induced alignment's cost over the optimum's (pgrealign's scores, hardset_strategies.pair_cost), its F1,
and where each of the first allele's letters goes against the optimum: `agree`; aligned elsewhere by under 3 letters
(`near`: gap placement), by a whole number of periods within 2 letters (`copy`: the same letters placed in another
copy of the repeat) or otherwise (`far`); aligned in one and gapped in the other (`gapped`). Then the totals per
label, and how much of the lost F1 sits in pairs whose induced alignment costs no more than the optimum.
"""
import collections
import csv
import gzip
import json
import os
import sys


def partners(ra, rb):
    """Each letter of the first sequence: its partner in the second, or None."""
    out = []
    j = 0
    for x, y in zip(ra, rb):
        if x != "-":
            out.append(j if y != "-" else None)
        if y != "-":
            j += 1
    return out


def classify(opt_a, msa_a, period):
    counts = collections.Counter()
    for o, m in zip(opt_a, msa_a):
        if o == m:
            counts["agree"] += 1
        elif o is None or m is None:
            counts["gapped"] += 1
        else:
            off = m - o
            if abs(off) < 3:
                counts["near"] += 1
            elif period > 0 and abs(off) >= period - 2 and min(abs(off) % period, period - abs(off) % period) <= 2:
                counts["copy"] += 1
            else:
                counts["far"] += 1
    return counts


def main():
    vntr_eval, work, rid = sys.argv[1:4]
    labels = sys.argv[4:]
    sys.path.insert(0, os.path.join(vntr_eval, "tools"))
    import hardset_strategies as hs  # noqa: E402
    from pgrealign import qc  # noqa: E402

    with open(os.path.join(vntr_eval, "hardset", "chr20.tsv")) as f:
        region = next(r for r in csv.DictReader(f, delimiter="\t") if r["id"] == rid)
    period = int(region["period"] or 0)
    with gzip.open(os.path.join(work, "optima", rid + ".jsonl.gz"), "rt") as f:
        optima = [json.loads(line) for line in f]
    msas = {}
    for lab in labels:
        rows, name = {}, None
        for line in open(os.path.join(work, "msa", lab, rid + ".fa")):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                rows[name] = ""
            else:
                rows[name] += line
        msas[lab] = rows
    print("%s: period %d, %d pairs scored" % (rid, period, sum(1 for r in optima if not r.get("skipped"))))
    totals = {lab: collections.Counter() for lab in labels}
    lost = {lab: [0, 0] for lab in labels}  # letters off the optimum in pairs costing no more / more
    for rec in optima:
        if rec.get("skipped"):
            continue
        oa, ob = rec["rows"]
        opt_cost = hs.pair_cost(oa, ob)
        o_part = partners(oa, ob)
        line = "  %-4s %-6s %-6s opt cost %6d |" % (rec["kind"], rec["a"], rec["b"], opt_cost)
        for lab in labels:
            ra, rb = qc.induced(msas[lab][rec["a"]], msas[lab][rec["b"]])
            cost = hs.pair_cost(ra, rb)
            o = qc.aligned_pairs(oa, ob)
            m = qc.aligned_pairs(ra, rb)
            f1 = 2 * len(o & m) / (len(o) + len(m)) if o or m else 1.0
            c = classify(o_part, partners(ra, rb), period)
            totals[lab].update(c)
            off = sum(v for k, v in c.items() if k != "agree")
            lost[lab][0 if cost <= opt_cost else 1] += off
            line += " %s F1 %.3f excess %5d near %4d copy %4d far %4d gapped %4d |" % (
                lab, f1, cost - opt_cost, c["near"], c["copy"], c["far"], c["gapped"])
        print(line)
    for lab in labels:
        t = totals[lab]
        n = sum(t.values())
        print("%-12s letters %d: agree %.1f%%, near %.1f%%, copy %.1f%%, far %.1f%%, gapped %.1f%% | off letters in "
              "pairs at or under the optimum's cost %d, above it %d" % (
                  lab, n, 100 * t["agree"] / n, 100 * t["near"] / n, 100 * t["copy"] / n, 100 * t["far"] / n,
                  100 * t["gapped"] / n, lost[lab][0], lost[lab][1]))


if __name__ == "__main__":
    main()
