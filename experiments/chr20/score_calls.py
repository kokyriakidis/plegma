#!/usr/bin/env python3
"""The chr20 arms' calls against HG002's truth (T2T-Q100 v1.1), scored as vg-call-eval scores them: small variants by
aardvark (genotype match), variants of 50 bp or more by truvari (multiallelics split, --pick ac), each inside the
truth's benchmark regions, and the same restricted to the tandem-repeat regions the arms realign and to vntr-eval's
45 hard-set regions.

usage: python3 -I experiments/chr20/score_calls.py <work dir> <truth dir> <vg-call-eval> <packages.jsonl.gz>
           <hard set chr20.tsv> <arm> [<arm> ...] [--refined]

Reuses vg-call-eval's scripts/wgs/bench_wgs.score_contig on a layout of its own per (arm, set):
<work dir>/calls/<arm>.<set>/chr20/{chr20.vcf.gz (the arm's short.vcf.gz), chr20.fa, truth.chr20.{smvar,stvar}.vcf.gz,
truth.chr20.{smvar,stvar}.bed (the benchmark bed, intersected with the set's regions)}. Needs bcftools, aardvark and
truvari (BCFTOOLS_DIR, AARDVARK_DIR on PATH; TRUVARI) and REF, CHM13's chr20 FASTA (indexed). Prints per set and arm
the small-variant F1 with its true and false counts and the structural F1; with --refined, the structural F1 again
after truvari refine (refined.py: the truth's and the calls' haplotypes written alike where they disagree), which no
representation changes.
"""
import csv
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path


def read_bed(path):
    out = []
    with open(path) as f:
        for line in f:
            w = line.split("\t")
            if len(w) >= 3 and w[0] == "chr20":
                out.append((int(w[1]), int(w[2])))
    out.sort()
    return out


def merged(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def intersect(a, b):
    """Two sorted lists of disjoint intervals: their intersection."""
    out = []
    i = j = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def write_bed(path, intervals):
    with open(path, "w") as f:
        for s, e in intervals:
            f.write("chr20\t%d\t%d\n" % (s, e))


def main():
    args = sys.argv[1:]
    refine = "--refined" in args
    args = [a for a in args if a != "--refined"]
    work, truth, vce, packages, hardset = (Path(x) for x in args[:5])
    arms = args[5:]
    sys.path.insert(0, str(vce / "scripts" / "wgs"))
    sys.path.insert(0, str(vce / "scripts" / "tier2"))
    sys.path.insert(0, str(vce / "scripts"))
    sys.path.insert(0, str(vce / "src"))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    os.environ["PATH"] = os.pathsep.join(
        x for x in (os.environ.get("BCFTOOLS_DIR"), os.environ.get("AARDVARK_DIR"), os.environ["PATH"]) if x)
    import bench_wgs  # noqa: E402
    import bench_metrics as bm  # noqa: E402

    truvari = os.environ.get("TRUVARI", str(vce / "work" / "truvari-venv" / "bin" / "truvari"))
    base = work / "calls"
    shared = base / "truth"
    shared.mkdir(parents=True, exist_ok=True)
    # CHM13's chr20 (REF, indexed), the reference every arm's calls are on.
    for ext in ("", ".fai"):
        if not (base / ("chr20.fa" + ext)).exists():
            os.symlink(os.environ["REF"] + ext, base / ("chr20.fa" + ext))
    # Truth slices for chr20, once.
    for kind in ("smvar", "stvar"):
        vcf = shared / ("truth.chr20.%s.vcf.gz" % kind)
        if not vcf.exists():
            subprocess.run(["bcftools", "view", "-r", "chr20", "-Oz", "-o", str(vcf),
                            str(truth / ("CHM13v2.0_HG2-T2TQ100-V1.1_%s.vcf.gz" % kind))], check=True)
            subprocess.run(["bcftools", "index", "-t", "-f", str(vcf)], check=True)
    # The sets of regions: the realigned tandem repeats (every package's CHM13 span) and the hard set.
    repeats = []
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p.get("contig") == "chr20" and p.get("ref_end"):
                repeats.append((int(p["ref_start"]), int(p["ref_end"])))
    with open(hardset) as f:
        hard = [(int(r["ref_start"]), int(r["ref_end"])) for r in csv.DictReader(f, delimiter="\t")]
    sets = {"chr20": None, "repeats": merged(repeats), "hardset": merged(hard)}
    beds = {}
    for kind in ("smvar", "stvar"):
        bench = merged(read_bed(truth / ("CHM13v2.0_HG2-T2TQ100-V1.1_%s.benchmark.bed" % kind)))
        for name, regions in sets.items():
            beds[(name, kind)] = bench if regions is None else intersect(bench, regions)
    print("benchmark letters per set: %s" % ", ".join(
        "%s %d small / %d structural" % (name, sum(e - s for s, e in beds[(name, "smvar")]),
                                         sum(e - s for s, e in beds[(name, "stvar")])) for name in sets))
    rows = []
    for name in sets:
        for arm in arms:
            d = base / ("%s.%s" % (arm, name))
            c = d / "chr20"
            c.mkdir(parents=True, exist_ok=True)
            for src, dst in ((work / ("arm_%s" % arm) / "short.vcf.gz", "chr20.vcf.gz"),
                             (work / ("arm_%s" % arm) / "short.vcf.gz.tbi", "chr20.vcf.gz.tbi"),
                             (work / "calls" / "chr20.fa", "chr20.fa"), (work / "calls" / "chr20.fa.fai", "chr20.fa.fai")):
                if not (c / dst).exists():
                    os.symlink(src, c / dst)
            for kind in ("smvar", "stvar"):
                for ext in ("vcf.gz", "vcf.gz.tbi"):
                    dst = c / ("truth.chr20.%s.%s" % (kind, ext))
                    if not dst.exists():
                        os.symlink(shared / ("truth.chr20.%s.%s" % (kind, ext)), dst)
                write_bed(c / ("truth.chr20.%s.bed" % kind), beds[(name, kind)])
            r = bench_wgs.score_contig(d, "chr20", "HG002", 8, truvari)
            a = bench_wgs.pick(r.get("aardvark"), "GT", "ALL")
            small = bm.aardvark_counts(a) if a else None
            t = r.get("truvari")
            sv = bm.truvari_counts(t) if t else None
            rows.append((name, arm, small, sv, r))
            def counts(x):
                return "TP %d FN %d FP %d" % (x.truth_tp, x.truth_fn, x.query_fp) if x else ""

            print("%-8s %-9s small F1 %s (%s) | structural F1 %s (%s) %s" % (
                name, arm, "%.4f" % small.f1 if small else "-", counts(small),
                "%.4f" % sv.f1 if sv else "-", counts(sv),
                "; ".join("%s: %s" % (k, r[k]) for k in ("error", "aardvark_error", "truvari_error") if k in r)),
                flush=True)
            if refine and sv:
                import refined  # noqa: E402

                tp, fn, fp, f1 = refined.refined(str(d / "score" / "chr20.truvari"))["refined"]
                print("%-8s %-9s refined structural F1 %.4f (TP %d FN %d FP %d)" % (name, arm, f1, tp, fn, fp),
                      flush=True)


if __name__ == "__main__":
    main()
