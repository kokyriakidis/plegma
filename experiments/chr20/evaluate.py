#!/usr/bin/env python3
"""The local harness, one entry point: chr20 tandem-repeat regions genotyped one at a time, as the chr20 arms
genotype the whole chromosome, to compare MSAs in minutes and see why each errs.

  prep    <chr20 work dir> <packages.jsonl.gz> <regions file> <prep dir> [--arms a,b,...] [--genome fa]
          [--extra reads.fq.gz,...] [--jobs n]
      each region's window of CHM13 (FLANK each side), its truth and the reads the arms map into it or its flanks
      (local_call.prep); with --extra and --genome, the reads the chr20 FASTQs lack that share K-mers only these
      regions hold (extra_reads.py, in place); then the prep's checks (prep_check.py: windows, mates, truth inside the
      scored interior, read depth, starved interiors).
  run     <prep dir> <packages.jsonl.gz> <out dir> <label>=<msas.jsonl.gz | plegma flags>[@ENV=v;ENV=v] ...
          [--regions file] [--jobs n] [--minimal]
      each label's calls (batch_call.py: its MSA, the graph, vg haplotypes sampling on the region's own reads, giraffe,
      vg call, the scores), with the environment switches after `@`; calls cached by the rows, so an unchanged MSA costs
      nothing. Writes <out dir>/<label>/qc.tsv (every step's numbers per region) and batch_calls.tsv. With --minimal,
      each graph only CHM13's row and the rows nearest HG002's truth (a diagnostic: do the reads and vg call get the
      truth when nothing else is there?).
  report  <prep dir> <packages.jsonl.gz> <out dir> <label> ... [--regions file]
      per region and label: the MSA's own properties, which need no truth (its width over its longest allele, its
      scattered changes: one change of one length at several places, its indels per allele against CHM13); the called
      haplotypes' edit distance to HG002's truth against the least the graph's rows allow (right when within 1%,
      however the records split them); the records' structural TP/FN/FP; the sampling and mapping numbers; totals.
  strata  <out dir> <prep dir> <packages.jsonl.gz> <label> [--others a,b] [--regions file] [--jobs n] [--out tsv]
      each structural error put in a class: on the MSA's side or not (strata.py).
  check   <out dir> <prep dir> <packages.jsonl.gz> <label> ... [--regions file] [--jobs n] [--out dir]
      the MSA judged without the genotyper: the truth's and the reads' genotypes written through the MSA and optimally,
      scored; split and scattered changes (msa_check.py).

`run` needs VG, GAFBASE, BCFTOOLS_DIR, AARDVARK_DIR, VCE, VNTR_PGREALIGN and, for plegma flags, PLEGMA; `prep` VG and
TRUTH.
"""
import csv
import gzip
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def option(args, name, default):
    if name in args:
        i = args.index(name)
        value = args[i + 1]
        del args[i : i + 2]
        return value
    return default


def prep(args):
    import extra_reads
    import local_call
    import prep_check

    arms = option(args, "--arms", "original,star,plegma,plegmaQ").split(",")
    genome = option(args, "--genome", "")
    extra = [x for x in option(args, "--extra", "").split(",") if x]
    jobs = option(args, "--jobs", "6")
    work, packages, regions, prep_dir = args[:4]
    local_call.prep([work, packages, regions, prep_dir] + arms)
    if extra and genome:
        extra_reads.main([prep_dir, packages, prep_dir, genome] + extra)
    prep_check.main([prep_dir, packages, "--jobs", jobs])


def run(args):
    minimal = "--minimal" in args
    args = [a for a in args if a != "--minimal"]
    regions = option(args, "--regions", "")
    jobs = option(args, "--jobs", "8")
    prep_dir, packages, out_dir = args[:3]
    for spec in args[3:]:
        label, rest = spec.split("=", 1)
        flags, _, env_text = rest.partition("@")
        env = dict(os.environ)
        for pair in filter(None, env_text.split(";")):
            k, v = pair.split("=", 1)
            env[k] = v
        cmd = ["python3", "-I", os.path.join(HERE, "batch_call.py"), prep_dir, packages, out_dir, "%s=%s" % (label, flags),
               "--jobs", jobs, "--vg-haplotypes", "region"] + (["--regions", regions] if regions else []) + (
                   ["--minimal"] if minimal else [])
        proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
        lines = [x for x in proc.stdout.splitlines() if x.startswith(label) or x.startswith("   ")]
        print("\n".join(x[:240] for x in lines) if lines else "%s failed: %s" % (label, proc.stderr[-400:]), flush=True)


def report(args):
    import haps
    import msa_check

    regions_file = option(args, "--regions", "")
    prep_dir, packages, out_dir = args[:3]
    labels = args[3:]
    qcs = {}
    for label in labels:
        path = os.path.join(out_dir, label, "qc.tsv")
        qcs[label] = {r["region"]: r for r in csv.DictReader(open(path), delimiter="\t")} if os.path.exists(path) else {}
    regions = (open(regions_file).read().split() if regions_file else
               sorted(set().union(*[set(q) for q in qcs.values()])))
    pkgs = {}
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] in regions:
                pkgs[p["id"]] = p
    totals = {label: [0, 0, 0, 0, 0, 0, 0, 0, 0.0] for label in labels}
    print("%-20s %-6s %10s %6s %5s %6s %15s %7s %13s" % ("region", "label", "tp/fn/fp", "width", "scat", "ev/al",
                                                         "called/least", "kept", "mapped/reads"))
    for rid in regions:
        if rid not in pkgs:
            continue
        pkg = pkgs[rid]
        longest = max(len(a["seq"]) for a in pkg["alleles"])
        try:
            meta, contig, truth = haps.region(out_dir, prep_dir, pkg, rid)
        except OSError:
            continue
        weight = {a["id"]: a["weight"] for a in pkg["alleles"]}
        for label in labels:
            q = qcs[label].get(rid)
            if q is None or q.get("error"):
                print("%-20s %-6s failed %s" % (rid, label, (q or {}).get("error", "missing")[:80]))
                continue
            tp, fn, fp = int(q["sv_tp"]), int(q["sv_fn"]), int(q["sv_fp"])
            try:
                _, sampled, called, pair = haps.called(out_dir, rid, label, meta, contig, truth)
                got = sum(haps.distance(called[h], truth[pair[h]]) for h in (0, 1))
                least = sum(haps.nearest(t, sampled, weight)[0] for t in truth)
            except (OSError, KeyError, ValueError):
                got = least = -1
            rows = json.load(open(os.path.join(out_dir, "%s.%s" % (label, rid), "rows.json")))
            scat = msa_check.scattered(rows, pkg["ref"]) if pkg.get("ref") in rows else -1
            width = int(q["columns"]) / longest if q.get("columns") else 0
            # Parsimony, with no truth: each row's indels of 10 letters or more against CHM13's row.
            ref_row = rows.get(pkg["ref"])
            events = [len(msa_check.indels(msa_check.records_from(ref_row, r, 0, "N" * (len(ref_row) + 1))))
                      for a, r in rows.items() if a != pkg["ref"]] if ref_row else []
            per_allele = sum(events) / max(len(events), 1)
            # The haplotypes called as near HG002's as the graph's rows allow (a score of the haplotypes, however
            # the records split them): within 1% of the truth's length past the least.
            right = got >= 0 and got <= least + max(20, (len(truth[0]) + len(truth[1])) // 100)
            print("%-20s %-6s %10s %6.1f %5d %6.2f %7d/%-7d %7s %6s/%-6s" % (
                rid, label, "%d/%d/%d" % (tp, fn, fp), width, scat, per_allele, got, least, q.get("kept_alleles", ""),
                q.get("mapped", ""), q.get("reads", "")))
            t = totals[label]
            t[0] += tp
            t[1] += fn
            t[2] += fp
            t[3] += max(got - least, 0)
            t[4] += max(scat, 0)
            t[5] += 1
            t[6] += right
            t[7] += width
            t[8] += per_allele
    print("totals (the MSA's own properties need no truth: width over the longest allele, scattered changes, indels per")
    print("allele; the haplotypes' and records' scores use HG002's)")
    for label, (tp, fn, fp, excess, scat, n, right, width, per_allele) in totals.items():
        p, r = tp / max(tp + fp, 1), tp / max(tp + fn, 1)
        print("  %-6s width %.2f scattered %d indels/allele %.2f | haplotypes right %d/%d, excess %d | records: structural "
              "TP %d FN %d FP %d F1 %.3f" % (label, width / max(n, 1), scat, per_allele / max(n, 1), right, n, excess,
                                             tp, fn, fp, 2 * p * r / max(p + r, 1e-9)))


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ("prep", "run", "report", "strata", "check"):
        sys.exit(__doc__)
    command, args = sys.argv[1], sys.argv[2:]
    if command == "prep":
        prep(args)
    elif command == "run":
        run(args)
    elif command == "report":
        report(args)
    elif command == "strata":
        import strata

        strata.main(args)
    else:
        import msa_check

        msa_check.main(args)


if __name__ == "__main__":
    main()
