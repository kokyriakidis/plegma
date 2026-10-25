#!/usr/bin/env python3
"""A structural score that does not hang on representation: truvari refine (phab) over a truvari bench.

truvari bench matches records by position, size and sequence, so a haplotype written as one 132-letter insertion
where the truth writes two of 66 scores 2 FN and 1 FP though the haplotype is the same. truvari refine realigns the
truth's and the calls' haplotypes in each region where they disagree (MAFFT) and writes both again the same way
before matching: what is left is the genotype's error, not the representation's.

usage: python3 -I experiments/chr20/refined.py <truvari bench dir> [...]
  refined(bench_dir) -> {"bench": (tp, fn, fp, f1), "refined": (tp, fn, fp, f1)}

The bench is rerun from the bench dir's params.json with the truth's sample renamed TRUTH and the calls' QUERY (both
are HG002, which phab's merged VCF cannot hold twice) and every contig's `_` made `.` (phab names a haplotype
<sample>_<1|2>_<contig>:<start>-<end> and msatovcf cuts it at its last `_`: a contig holding one loses every
genotype), then refined with --use-original-vcfs. Needs bcftools
(BCFTOOLS_DIR), truvari (TRUVARI or vg-call-eval's venv under VCE) and MAFFT (MAFFT_DIR, by default the vntr
micromamba environment's bin).
"""
import json
import os
import shutil
import subprocess
import sys


def counts(path):
    d = json.load(open(path))
    return d["TP-base"], d["FN"], d["FP"], d["f1"] or 0.0


def refined(bench_dir, threads=8):
    params = json.load(open(os.path.join(bench_dir, "params.json")))
    truvari = os.environ.get("TRUVARI") or os.path.join(os.environ["VCE"], "work", "truvari-venv", "bin", "truvari")
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join(x for x in (os.environ.get("MAFFT_DIR", os.path.expanduser(
        "~/micromamba/envs/vntr/bin")), os.environ.get("BCFTOOLS_DIR"), env["PATH"]) if x)
    out = bench_dir.rstrip("/") + ".refined"
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    contigs = [line.split("\t")[0] for line in open(params["reference"] + ".fai")]
    chrs = os.path.join(out, "contigs.map")
    with open(chrs, "w") as f:
        f.writelines("%s\t%s\n" % (c, c.replace("_", ".")) for c in contigs)
    # The reference again under the new names, each contig on one line, and its index.
    reference = os.path.join(out, "reference.fa")
    with open(params["reference"]) as src, open(reference, "w") as dst, open(reference + ".fai", "w") as fai:
        at = 0
        name, parts = None, []

        def flush():
            nonlocal at
            if name is None:
                return
            seq = "".join(parts)
            dst.write(">%s\n%s\n" % (name, seq))
            fai.write("%s\t%d\t%d\t%d\t%d\n" % (name, len(seq), at + len(name) + 2, len(seq), len(seq) + 1))
            at += len(name) + 2 + len(seq) + 1

        for line in src:
            if line.startswith(">"):
                flush()
                name, parts = line[1:].split()[0].replace("_", "."), []
            else:
                parts.append(line.strip())
        flush()
    bed = os.path.join(out, "include.bed")
    with open(params["includebed"]) as src, open(bed, "w") as dst:
        for line in src:
            w = line.split("\t", 1)
            dst.write(w[0].replace("_", ".") + "\t" + w[1])
    renamed = {}
    for key, sample in (("base", "TRUTH"), ("comp", "QUERY")):
        names = os.path.join(out, key + ".names")
        with open(names, "w") as f:
            f.write("%s\n" % sample)
        path = os.path.join(out, key + ".vcf.gz")
        subprocess.run("bcftools reheader -s %s %s | bcftools annotate --rename-chrs %s -Oz -o %s" % (
            names, params[key], chrs, path), shell=True, env=env, check=True)
        subprocess.run(["bcftools", "index", "-t", "-f", path], env=env, check=True)
        renamed[key] = path
    bench = os.path.join(out, "bench")
    cmd = [truvari, "bench", "-b", renamed["base"], "-c", renamed["comp"], "-f", reference, "-o", bench,
           "--bSample", "TRUTH", "--cSample", "QUERY", "--includebed", bed, "--pick", params["pick"],
           "--refdist", str(params["refdist"]),
           "--pctseq", str(params["pctseq"]), "--pctsize", str(params["pctsize"]), "--sizemin", str(params["sizemin"]),
           "--sizefilt", str(params["sizefilt"]), "--sizemax", str(params["sizemax"])]
    subprocess.run(cmd, env=env, check=True, capture_output=True)
    subprocess.run([truvari, "refine", "-u", "-a", "mafft", "-t", str(threads), bench], env=env, check=True,
                   capture_output=True)
    return {"bench": counts(os.path.join(bench, "summary.json")),
            "refined": counts(os.path.join(bench, "refine.variant_summary.json"))}


def main(argv=None):
    for d in (sys.argv[1:] if argv is None else argv):
        r = refined(d)
        print("%s\n  bench   TP %d FN %d FP %d F1 %.4f\n  refined TP %d FN %d FP %d F1 %.4f" % (
            (d,) + r["bench"] + r["refined"]))


if __name__ == "__main__":
    main()
