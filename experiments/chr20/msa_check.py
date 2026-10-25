#!/usr/bin/env python3
"""How good an MSA is for calling, with the genotyper taken out: the genotype chosen without vg call, then written
as calls by the MSA alone, and scored as vg-call-eval scores calls.

usage: python3 -I experiments/chr20/msa_check.py <batch_call out dir> <prep dir> <packages.jsonl.gz> <label>...
           [--regions file] [--jobs n] [--out dir]

For each label, per region its graph's rows (<out dir>/<label>.<region>/rows.json, what vg call had) and two
genotypes from those rows:
  truth   the rows nearest HG002's two truth haplotypes: the best genotype this graph holds
  oracle  the pair of rows the region's reads support best (oracle.py: K-mer counts, depth as well as presence)
each written two ways:
  msa      every difference between the row and CHM13's row in the MSA's columns, one record per run of differing
           columns: what a caller reports through this MSA's graph
  optimal  the row aligned exactly to CHM13's allele on its own (plegma ref on the pair, the binary in PLEGMA): how
           the allele's differences read at their best, whatever MSA holds it
and scored together over all regions (batch_call.score: truvari for >= 50 bp, aardvark for small variants). So
  truth x optimal   the ceiling of these alleles (what the panel and the sampling allow)
  truth x msa       the ceiling of this MSA: below truth x optimal by what its representation costs
  oracle x msa      what a genotyper using all the reads' information would score through this MSA
and vg call's own score (batch_calls.tsv) below oracle x msa by what the genotyper loses.
Per region, also: `split`, rows whose indels (>= 10 letters) through the MSA outnumber their optimal alignment's
(a change cut in pieces); `scattered`, changes of one kind and length (within a tenth) that different rows hold at
different places against CHM13 (one change as several bubbles); `oracle_ok`, whether the oracle's pair is the truth
pair. Writes <out>/<label>.tsv per region and prints the totals.
"""
import collections
import gzip
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import batch_call  # noqa: E402
import haps  # noqa: E402
import local_call  # noqa: E402
import oracle  # noqa: E402

EVENT = 10


def records_from(ref_row, row, offset, contig):
    """(pos, ref, alt) of each run of columns where `row` differs from `ref_row`, on the contig (1-based, the row's
    first letter at `offset`), an indel anchored on the reference letter before it."""
    out = []
    at = offset  # the next reference letter's 0-based position on the contig
    run_ref, run_alt, start = [], [], None
    for x, y in zip(ref_row, row):
        if x == "-" and y == "-":
            continue
        if x == y:
            if start is not None:
                out.append((start, "".join(run_ref), "".join(run_alt)))
                run_ref, run_alt, start = [], [], None
            at += 1
            continue
        if start is None:
            start = at
        if x != "-":
            run_ref.append(x)
            at += 1
        if y != "-":
            run_alt.append(y)
    if start is not None:
        out.append((start, "".join(run_ref), "".join(run_alt)))
    recs = []
    for start, r, a in out:
        if len(r) == len(a):
            recs.append((start + 1, r, a))
        else:
            anchor = contig[start - 1]
            recs.append((start, anchor + r, anchor + a))
    return recs


def indels(recs):
    return [(len(a) - len(r), p) for p, r, a in recs if abs(len(a) - len(r)) >= EVENT]


def vcf_lines(rid, hap_records):
    """Phased GT records of two haplotypes' record lists."""
    carriers = collections.defaultdict(set)
    for h, recs in enumerate(hap_records):
        for rec in recs:
            carriers[rec].add(h)
    lines = []
    for (pos, r, a), hs in sorted(carriers.items()):
        gt = "1|1" if hs == {0, 1} else ("1|0" if hs == {0} else "0|1")
        lines.append("%s\t%d\t.\t%s\t%s\t60\tPASS\t.\tGT\t%s\n" % (rid, pos, r, a, gt))
    return lines


def optimal_rows(ref_seq, seq):
    """The two rows of plegma ref on (reference, seq): the exact alignment to the reference."""
    if seq == ref_seq:
        return ref_seq, seq
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "in.fa")
        with open(path, "w") as f:
            f.write(">r\n%s\n>q\n%s\n" % (ref_seq, seq))
        out = subprocess.run([os.environ["PLEGMA"], "ref", "-r", "1", path], capture_output=True, text=True,
                             check=True).stdout
    rows, name = {}, None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            rows[name] = []
        elif name:
            rows[name].append(line.strip().upper())
    return "".join(rows["r"]), "".join(rows["q"])


def same_length(a, b):
    return abs(a - b) <= max(5, max(abs(a), abs(b)) // 10)


def scattered(rows, ref_id):
    """The changes that different rows hold at different places against CHM13: events of one kind and length (within
    a tenth) within reach of each other (twice the length, 100 letters at least), clustered, the clusters holding more
    than one place counted."""
    ref_row = rows[ref_id]
    events = set()
    for a, r in rows.items():
        if a == ref_id:
            continue
        for d, p in indels(records_from(ref_row, r, 0, "N" * (len(ref_row) + 1))):
            events.add((d > 0, d, p))
    events = sorted(events)
    parent = list(range(len(events)))

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, (k1, d1, p1) in enumerate(events):
        for j in range(i + 1, len(events)):
            k2, d2, p2 = events[j]
            if k1 == k2 and same_length(d1, d2) and abs(p1 - p2) <= max(2 * abs(d1), 100):
                parent[root(i)] = root(j)
    places = collections.defaultdict(set)
    for i, (_, _, p) in enumerate(events):
        places[root(i)].add(p)
    return sum(1 for ps in places.values() if len(ps) > 1)


def check(task):
    out_dir, prep_dir, pkg, label = task
    rid = pkg["id"]
    try:
        meta, contig, truth = haps.region(out_dir, prep_dir, pkg, rid)
        rows = json.load(open(os.path.join(out_dir, "%s.%s" % (label, rid), "rows.json")))
        ref_id = pkg["ref"]
        ref_seq = rows[ref_id].replace("-", "")
        offset = len(meta["left"])
        if contig[offset : offset + len(ref_seq)] != ref_seq:
            return rid, None, "contig is not left + reference + right"
        pool = {}
        for a, r in rows.items():
            pool.setdefault(meta["left"] + r.replace("-", "") + meta["right"], []).append(a)
        weight = {a["id"]: a["weight"] for a in pkg["alleles"]}
        truth_pair = [haps.nearest(t, pool, weight)[1] for t in truth]
        oa, ob, margin, sa, sb = oracle.best_pair(rows, meta["left"], meta["right"], prep_dir, rid)
        straight = haps.distance(sa, truth[0]) + haps.distance(sb, truth[1])
        crossed = haps.distance(sa, truth[1]) + haps.distance(sb, truth[0])
        oracle_pair = [oa, ob] if straight <= crossed else [ob, oa]
        best = sum(haps.distance(meta["left"] + rows[a].replace("-", "") + meta["right"], truth[h])
                   for h, a in enumerate(truth_pair))
        oracle_ok = min(straight, crossed) <= best + 20  # as oracle.py's evaluation counts it
        written = {}
        split = 0
        cache = {}
        for name, pair in (("truth", truth_pair), ("oracle", oracle_pair)):
            for how in ("msa", "optimal"):
                hap_records = []
                for a in pair:
                    if how == "msa":
                        recs = records_from(rows[ref_id], rows[a], offset, contig)
                    else:
                        if a not in cache:
                            r0, r1 = optimal_rows(ref_seq, rows[a].replace("-", ""))
                            cache[a] = records_from(r0, r1, offset, contig)
                        recs = cache[a]
                    hap_records.append(recs)
                written[name, how] = vcf_lines(rid, hap_records)
        for a in set(truth_pair) | set(oracle_pair):
            through = indels(records_from(rows[ref_id], rows[a], offset, contig))
            if len(through) > len(indels(cache[a])):
                split += 1
        note = {"oracle_ok": oracle_ok, "margin": margin, "split": split, "scattered": scattered(rows, ref_id),
                "truth_pair": truth_pair, "oracle_pair": oracle_pair}
        return rid, written, note
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        return rid, None, "%s: %s" % (type(error).__name__, error)


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)

    def option(name, default):
        if name in args:
            i = args.index(name)
            value = args[i + 1]
            del args[i : i + 2]
            return value
        return default

    jobs = int(option("--jobs", "8"))
    regions = option("--regions", "")
    report = option("--out", "")
    out_dir, prep_dir, packages = args[:3]
    labels = args[3:]
    report = report or os.path.join(out_dir, "msa_check")
    os.makedirs(report, exist_ok=True)
    ids = local_call.region_ids(regions) if regions else sorted(os.listdir(prep_dir))
    pkgs = local_call.packages_of(packages, set(ids))
    vce = os.environ["VCE"]
    for sub in ("scripts/wgs", "scripts/tier2", "scripts", "src"):
        sys.path.insert(0, os.path.join(vce, sub))
    os.environ["PATH"] = os.pathsep.join(x for x in (os.environ.get("BCFTOOLS_DIR"), os.environ.get("AARDVARK_DIR"),
                                                     os.environ["PATH"]) if x)
    truvari = os.environ.get("TRUVARI", os.path.join(vce, "work", "truvari-venv", "bin", "truvari"))
    header = open(os.path.join(out_dir, "cache", "header.vcf")).readlines()
    for label in labels:
        tasks = [(out_dir, prep_dir, pkgs[i], label) for i in ids
                 if i in pkgs and os.path.exists(os.path.join(out_dir, "%s.%s" % (label, i), "rows.json"))]
        results, notes, failed = {}, {}, {}
        with Pool(jobs) as pool:
            for rid, written, note in pool.imap_unordered(check, tasks):
                if written is None:
                    failed[rid] = note
                else:
                    results[rid], notes[rid] = written, note
        order = [i for i in ids if i in results]
        vg = {}
        path = os.path.join(out_dir, label, "batch_calls.tsv")
        if os.path.exists(path):
            for line in open(path).readlines()[1:]:
                w = line.rstrip("\n").split("\t")
                vg[w[0]] = (int(w[1]), int(w[2]), int(w[3]))
        print("%s: %d regions (%d failed); oracle's pair the truth's in %d; regions with split indels %d (%d rows), "
              "with scattered changes %d (%d changes)" % (
                  label, len(order), len(failed), sum(notes[r]["oracle_ok"] for r in order),
                  sum(1 for r in order if notes[r]["split"]), sum(notes[r]["split"] for r in order),
                  sum(1 for r in order if notes[r]["scattered"]), sum(notes[r]["scattered"] for r in order)))
        per_region = collections.defaultdict(dict)
        for name in ("truth", "oracle"):
            for how in ("msa", "optimal"):
                r, sv, small, per, _ = batch_call.score(order, {rid: results[rid][name, how] for rid in order}, header,
                                                        prep_dir, os.path.join(report, "%s.%s.%s" % (label, name, how)),
                                                        jobs, truvari)
                print("  %-7s x %-8s structural TP %4d FN %4d FP %4d F1 %.4f | small TP %5d FN %5d FP %5d F1 %.4f" % (
                    name, how, sv.truth_tp, sv.truth_fn, sv.query_fp, sv.f1, small.truth_tp, small.truth_fn,
                    small.query_fp, small.f1), flush=True)
                for rid in order:
                    p = per.get(rid, collections.Counter())
                    per_region[rid]["%s_%s" % (name, how)] = "%d/%d/%d" % (p["tp"], p["fn"], p["fp"])
        tp = sum(vg[r][0] for r in order if r in vg)
        fn = sum(vg[r][1] for r in order if r in vg)
        fp = sum(vg[r][2] for r in order if r in vg)
        prec, rec = tp / max(tp + fp, 1), tp / max(tp + fn, 1)
        print("  vg call            structural TP %4d FN %4d FP %4d F1 %.4f (its own run, these regions)" % (
            tp, fn, fp, 2 * prec * rec / max(prec + rec, 1e-9)))
        with open(os.path.join(report, label + ".tsv"), "w") as f:
            cols = ["truth_optimal", "truth_msa", "oracle_optimal", "oracle_msa"]
            f.write("region\t" + "\t".join(cols) + "\tvg_call\toracle_ok\tsplit\tscattered\ttruth_pair\toracle_pair\n")
            for rid in order:
                n = notes[rid]
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    rid, "\t".join(per_region[rid][c] for c in cols), "%d/%d/%d" % vg.get(rid, (0, 0, 0)),
                    int(n["oracle_ok"]), n["split"], n["scattered"], ",".join(n["truth_pair"]) + "\t" +
                    ",".join(n["oracle_pair"])))
        for rid, error in sorted(failed.items())[:5]:
            print("   failed %s: %s" % (rid, str(error)[:200]))


if __name__ == "__main__":
    main()
