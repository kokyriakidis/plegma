#!/usr/bin/env python3
"""How HG002's truth reads when written by an aligner: its two haplotypes per region written back as records against
CHM13 and scored against the truth itself, so what is lost is representation alone (the haplotypes are exact).

usage: python3 -I experiments/chr20/representation.py <batch_call out dir> <prep dir> <packages.jsonl.gz>
           [--regions file] [--jobs n] [--out dir]

Writers, per haplotype (the region's window with the truth applied, against the window's CHM13):
  optimal   the interior aligned exactly to CHM13's allele (plegma ref on the pair, the binary in PLEGMA): how ref writes
  minimap2  the whole window aligned by minimap2 at dipcall's settings (-x asm5 -z200000,10000 -r2k; MINIMAP2), the
            truth's own pipeline's aligner, its primary alignment
  chained   the interior aligned as minimap2 pins it (mm2like.py: its seeds, chain and pins, karan between pins)
  chained_window  the same on the whole window, as minimap2 sees it
  chained_context the interior aligned between the pins the whole window finds inside it (the flanks seed, but are
            not aligned: what ref can do with the flanks given)
  pinned    plegma ref -F (the Mojo port, plegma/pins.mojo) on the interior with the window's flanks
Prints each writer's structural and small TP/FN/FP; writes per region tsv to --out.
"""
import collections
import os
import re
import subprocess
import sys
import tempfile
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import batch_call  # noqa: E402
import haps  # noqa: E402
import local_call  # noqa: E402
import mm2like  # noqa: E402
import msa_check  # noqa: E402

MM2_OPTIONS = ["-c", "-x", "asm5", "-z200000,10000", "-r2k"]
WRITERS = ("optimal", "minimap2", "chained", "chained_window", "chained_context", "pinned")


def minimap2_rows(ref, seq):
    """(ref row, seq row, ref offset) of minimap2's primary alignment of `seq` to `ref`, None when it has none."""
    with tempfile.TemporaryDirectory() as tmp:
        r, q = os.path.join(tmp, "r.fa"), os.path.join(tmp, "q.fa")
        with open(r, "w") as f:
            f.write(">r\n%s\n" % ref)
        with open(q, "w") as f:
            f.write(">q\n%s\n" % seq)
        out = subprocess.run([os.environ["MINIMAP2"]] + MM2_OPTIONS + [r, q], capture_output=True, text=True,
                             check=True).stdout
    for line in out.splitlines():
        w = line.split("\t")
        if "tp:A:P" not in w or w[4] != "+":
            continue
        q_start, t_start = int(w[2]), int(w[7])
        cigar = next(x[5:] for x in w[12:] if x.startswith("cg:Z:"))
        rr, qq = [], []
        i, j = t_start, q_start
        for n, op in re.findall(r"(\d+)([MIDNX=])", cigar):
            n = int(n)
            if op in "M=X":
                rr.append(ref[i : i + n])
                qq.append(seq[j : j + n])
                i += n
                j += n
            elif op == "I":
                rr.append("-" * n)
                qq.append(seq[j : j + n])
                j += n
            else:
                rr.append(ref[i : i + n])
                qq.append("-" * n)
                i += n
        return "".join(rr), "".join(qq), t_start
    return None


def pinned_rows(ref, seq, left, right):
    """The two rows of plegma ref -F on (reference, seq), the flanks given."""
    if seq == ref:
        return ref, seq
    with tempfile.TemporaryDirectory() as tmp:
        path, flanks = os.path.join(tmp, "in.fa"), os.path.join(tmp, "flanks.fa")
        with open(path, "w") as f:
            f.write(">r\n%s\n>q\n%s\n" % (ref, seq))
        with open(flanks, "w") as f:
            f.write(">left\n%s\n>right\n%s\n" % (left, right))
        out = subprocess.run([os.environ["PLEGMA"], "ref", "-F", flanks, "-r", "1", path], capture_output=True,
                             text=True, check=True).stdout
    rows, name = {}, None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            rows[name] = []
        elif name:
            rows[name].append(line.strip().upper())
    return "".join(rows["r"]), "".join(rows["q"])


def write(task):
    out_dir, prep_dir, pkg = task
    rid = pkg["id"]
    try:
        meta, contig, truth = haps.region(out_dir, prep_dir, pkg, rid)
        lo, hi = len(meta["left"]), len(contig) - len(meta["right"])
        written = {}
        chained = mm2like.chained_rows({h: (contig[lo:hi], truth[h][lo : len(truth[h]) - len(meta["right"])])
                                        for h in (0, 1)})
        chained_window = mm2like.chained_rows({h: (contig, truth[h]) for h in (0, 1)})
        chained_context = mm2like.chained_rows(
            {h: (contig[lo:hi], truth[h][lo : len(truth[h]) - len(meta["right"])]) for h in (0, 1)},
            context={h: (meta["left"], meta["right"]) for h in (0, 1)})
        for how in WRITERS:
            hap_records = []
            for n, h in enumerate(truth):
                if how == "chained":
                    r0, r1 = chained[n]
                    hap_records.append(msa_check.records_from(r0, r1, lo, contig))
                elif how == "chained_context":
                    r0, r1 = chained_context[n]
                    hap_records.append(msa_check.records_from(r0, r1, lo, contig))
                elif how == "pinned":
                    interior = h[lo : len(h) - len(meta["right"])]
                    r0, r1 = pinned_rows(contig[lo:hi], interior, meta["left"], meta["right"])
                    hap_records.append(msa_check.records_from(r0, r1, lo, contig))
                elif how == "chained_window":
                    r0, r1 = chained_window[n]
                    hap_records.append(msa_check.records_from(r0, r1, 0, contig))
                elif how == "optimal":
                    interior = h[lo : len(h) - len(meta["right"])]
                    r0, r1 = msa_check.optimal_rows(contig[lo:hi], interior)
                    hap_records.append(msa_check.records_from(r0, r1, lo, contig))
                else:
                    got = minimap2_rows(contig, h)
                    if got is None:
                        return rid, None, "no minimap2 alignment"
                    r0, r1, start = got
                    hap_records.append(msa_check.records_from(r0, r1, start, contig))
            written[how] = msa_check.vcf_lines(rid, hap_records)
        return rid, written, ""
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
    out_dir, prep_dir, packages = args[:3]
    report = option("--out", os.path.join(out_dir, "representation"))
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
    results, failed = {}, {}
    with Pool(jobs) as pool:
        for rid, written, error in pool.imap_unordered(write, [(out_dir, prep_dir, pkgs[i]) for i in ids if i in pkgs]):
            if written is None:
                failed[rid] = error
            else:
                results[rid] = written
    order = [i for i in ids if i in results]
    print("%d regions (%d failed)" % (len(order), len(failed)))
    per_region = collections.defaultdict(dict)
    for how in WRITERS:
        _, sv, small, per, _ = batch_call.score(order, {rid: results[rid][how] for rid in order}, header, prep_dir,
                                                os.path.join(report, how), jobs, truvari)
        print("  %-8s structural TP %4d FN %4d FP %4d F1 %.4f | small TP %5d FN %5d FP %5d F1 %.4f" % (
            how, sv.truth_tp, sv.truth_fn, sv.query_fp, sv.f1, small.truth_tp, small.truth_fn, small.query_fp,
            small.f1), flush=True)
        for rid in order:
            p = per.get(rid, collections.Counter())
            per_region[rid][how] = "%d/%d/%d" % (p["tp"], p["fn"], p["fp"])
    with open(os.path.join(report, "representation.tsv"), "w") as f:
        f.write("region\t" + "\t".join(WRITERS) + "\n")
        for rid in order:
            f.write(rid + "\t" + "\t".join(per_region[rid][how] for how in WRITERS) + "\n")
    for rid, error in sorted(failed.items())[:5]:
        print("   failed %s: %s" % (rid, error[:200]))


if __name__ == "__main__":
    main()
