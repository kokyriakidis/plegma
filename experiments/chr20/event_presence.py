#!/usr/bin/env python3
"""For each truth event of the local genotyping check (local_call.py run), and each spurious call: whether the panel
holds it, whether each MSA holds it, so a wrong call is a panel's, an MSA's or the genotyper's.

usage: python3 -I experiments/chr20/event_presence.py <packages.jsonl.gz> <prep dir> <run dir> <regions a,b|file>
           <label> [<label> ...]

A truth variant of MIN letters or more (the region's T2T-Q100 structural VCF, as prep sliced it) is `in the panel`
when some allele's optimal alignment to CHM13's allele (pgrealign's, fragmentation.py's cached gaps) has a gap on its
side (a deletion: CHM13 letters unaligned; an insertion: the allele's) within NEAR letters of it and SIZE of its size;
`in the MSA` when the same holds for the alignment of CHM13's row and the allele's that the label's MSA induces
(rows.json in the run dir). A spurious call (truvari's fp) is likewise `a path of the MSA` when some allele's induced
gap matches it. Prints per region, label and event the verdict, and the totals: truth events missed although the MSA
holds them (mapping or genotyping), missed and not in the MSA though in the panel (the MSA's), not in the panel; and
spurious calls the MSA's paths hold (a wrong allele chosen) or do not (a path the graph made). Needs VNTR_PGREALIGN,
VNTR_ABPOA and FRAG_CACHE.
"""
import collections
import gzip
import json
import os
import sys

MIN = 50
NEAR = 200
SIZE = 0.2
HERE = os.path.dirname(os.path.abspath(__file__))


def events(path):
    """(local position of the variant's first letter, side D or I, size) of a VCF's records of MIN letters or more."""
    out = []
    if not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            w = line.split("\t")
            ref, alts = w[3], w[4].split(",")
            for alt in alts:
                change = len(alt) - len(ref)
                if abs(change) >= MIN:
                    out.append((int(w[1]), "I" if change > 0 else "D", abs(change)))
    return out


def matches(event, gaps, offset):
    pos, side, size = event
    for s, start, length in gaps:
        if s == side and abs(length - size) <= SIZE * size and abs(start + offset - pos) <= NEAR + size:
            return True
    return False


def main():
    packages, prep, run, regions = sys.argv[1:5]
    labels = sys.argv[5:]
    sys.path.insert(0, HERE)
    from fragmentation import gaps, induced, optimal_gaps  # noqa: E402

    ids = open(regions).read().split() if os.path.exists(regions) else regions.split(",")
    wanted = set(ids)
    pkgs = {}
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] in wanted:
                pkgs[p["id"]] = p
    totals = {label: collections.Counter() for label in labels}
    for rid in ids:
        pkg = pkgs.get(rid)
        if pkg is None:
            continue
        meta = json.load(open(os.path.join(prep, rid, "region.json")))
        # CHM13's interior starts after the left flank: local position = interior letter + flank + 1 (VCF, 1-based).
        offset = len(meta["left"]) + 1
        ref = next(a for a in pkg["alleles"] if a["id"] == pkg["ref"])
        optimal = optimal_gaps(pkg, ref, os.path.join(os.environ["FRAG_CACHE"], rid + ".json"))
        panel = [[tuple(g) for g in v] for v in optimal.values()]
        truth = events(os.path.join(prep, rid, "truth.chr20.stvar.vcf.gz"))
        for label in labels:
            d = os.path.join(run, "%s.%s" % (label, rid))
            rows_path = os.path.join(d, "rows.json")
            if not os.path.exists(rows_path):
                continue
            rows = json.load(open(rows_path))
            induced_gaps = [gaps(*induced(rows[ref["id"]], rows[a["id"]])) for a in pkg["alleles"]
                            if a["id"] != ref["id"] and a["id"] in rows]
            tdir = os.path.join(d, "score", "chr20.truvari")
            missed = set(events(os.path.join(tdir, "fn.vcf.gz")))
            spurious = events(os.path.join(tdir, "fp.vcf.gz"))
            t = totals[label]
            for e in truth:
                in_panel = any(matches(e, g, offset) for g in panel)
                in_msa = any(matches(e, g, offset) for g in induced_gaps)
                state = "missed" if e in missed else "found"
                kind = "%s, %s" % ("in the MSA" if in_msa else "not in the MSA",
                                   "in the panel" if in_panel else "not in the panel")
                t["truth " + state + " | " + kind] += 1
                print("%-20s %-8s truth %s%d@%d %s: %s" % (rid, label, e[1].replace("D", "DEL").replace("I", "INS"),
                                                          e[2], e[0], state, kind))
            for e in spurious:
                on_path = any(matches(e, g, offset) for g in induced_gaps)
                in_panel = any(matches(e, g, offset) for g in panel)
                t["spurious | %s, %s" % ("an MSA path" if on_path else "no MSA path",
                                         "in the panel" if in_panel else "not in the panel")] += 1
                print("%-20s %-8s spurious %s%d@%d: %s, %s" % (
                    rid, label, e[1].replace("D", "DEL").replace("I", "INS"), e[2], e[0],
                    "an MSA path" if on_path else "no MSA path", "in the panel" if in_panel else "not in the panel"))
    for label, t in totals.items():
        print("%s:" % label)
        for k, v in sorted(t.items()):
            print("   %-60s %d" % (k, v))


if __name__ == "__main__":
    main()
