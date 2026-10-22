#!/usr/bin/env python3
"""The targets to extract: every catalogue repeat of 1 kb or more, vntr-eval's VNTRs of 1 kb or more and its hard
set, as one pgrealign catalogue.

usage: python3 -I experiments/regions/targets.py <suite dir>

Reads `<suite dir>/manifest.tsv` (gather.py) and writes `<suite dir>/targets_ge1kb.tsv` in pgrealign's catalogue
format (its header line, then chrom, start, end, id, class, period, motif, n_members, member_classes, sources):
catalogue entries as they are, vntr-eval's own entries under their own ids with `sources` naming their set. pgrealign
regions merges overlapping targets into one region, so a repeat in several lists becomes one region whose `targets`
names them all.
"""
import csv
import os
import sys

HEADER = "#chrom\tstart\tend\tid\tclass\tperiod\tmotif\tn_members\tmember_classes\tsources\n"
ORDER = {"chr%s" % c: k for k, c in enumerate([str(n) for n in range(1, 23)] + ["X", "Y", "M"])}


def main():
    suite = sys.argv[1]
    rows = []
    seen = set()
    with open(os.path.join(suite, "manifest.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["set"] == "packages" or int(r["end"]) - int(r["start"]) < 1000 and r["set"] == "catalog":
                continue
            if r["id"] in seen:
                continue
            seen.add(r["id"])
            kind = r["class"].split(",")[0] if r["class"] not in (".", "") else "VNTR"
            period = r["period"] if r["period"] not in ("", ".") else "0"
            period = str(int(float(period)))
            rows.append((r["contig"], int(r["start"]), int(r["end"]), r["id"], kind, period, r["motif"] or ".",
                         "1", "%s:1" % kind, r["set"] + ":1"))
    rows.sort(key=lambda x: (ORDER.get(x[0], 99), x[0], x[1], x[3]))
    with open(os.path.join(suite, "targets_ge1kb.tsv"), "w") as f:
        f.write(HEADER)
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(len(rows), "targets")


if __name__ == "__main__":
    main()
