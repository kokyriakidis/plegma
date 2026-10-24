#!/usr/bin/env python3
"""Genotyping one region at a time: an MSA's graph against HG002's truth in seconds per region, instead of the whole
chr20 pipeline's hours, for comparing MSAs (its numbers are relative: the haplotypes are the region's alleles, not 32
sampled ones, and the flanks are CHM13's).

  prep:  python3 -I experiments/chr20/local_call.py prep <chr20 work dir> <packages.jsonl.gz> <regions a,b|file>
             <prep dir> <arm> [<arm> ...]
  run:   python3 -I experiments/chr20/local_call.py run <prep dir> <packages.jsonl.gz> <out dir>
             <label>=<msas.jsonl.gz | plegma flags> [...] [--jobs n] [--regions a,b|file]

`prep`, once for a set of regions: the reads that the chr20 arms (`arm_<arm>`) align into each region or its flanks
(its nodes and CHM13's within FLANK letters of them, as region_reads.py finds them), both mates, written per region from the arms' FASTQs; CHM13's sequence over the region
and FLANK letters each side; and the truth over that interval (the T2T-Q100 small and structural VCFs and benchmark
beds that score_calls.py sliced into <work dir>/calls/truth), shifted to the interval's coordinates.

`run`, per label and region (in parallel): the MSA's rows (a pgrealign msas.jsonl.gz's, or plegma's on the distinct
alleles), every allele flanked by CHM13's flanks, induced into a graph as pgrealign replace does (induce.py; each allele
a haplotype walk per haplotype spelling it, up to 8, CHM13's the reference); `vg autoindex` for giraffe, `vg giraffe`
(paired as the arms map, then the singles), GAF-Base and `vg call` with the arms' flags, and the calls scored by
vg-call-eval's score_contig (aardvark for small variants, truvari for 50 bp or more). Prints per region and label the
structural and small true, missed and spurious counts, and the totals; writes <out dir>/local_calls.tsv.

Needs VG, GAFBASE, BCFTOOLS_DIR, AARDVARK_DIR, VCE (vg-call-eval), VNTR_PGREALIGN and, for flags, PLEGMA.
"""
import collections
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

FLANK = 1000
MAX_COPIES = 8
HERE = os.path.dirname(os.path.abspath(__file__))


def region_ids(text):
    return open(text).read().split() if os.path.exists(text) else text.split(",")


def packages_of(path, wanted):
    out = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] in wanted:
                out[p["id"]] = p
    return out


def fastq_records(path):
    with gzip.open(path, "rt") as f:
        while True:
            head = f.readline()
            if not head:
                return
            seq, plus, qual = f.readline(), f.readline(), f.readline()
            yield head[1:].split()[0].split("/")[0], head + seq + plus + qual


def fastq_records_plain(path):
    """(name, record text) of an uncompressed FASTQ."""
    with open(path) as f:
        while True:
            head = f.readline()
            if not head:
                return
            seq, plus, qual = f.readline(), f.readline(), f.readline()
            yield head[1:].split()[0].split("/")[0], head + seq + plus + qual


def vcf_slice(src, dst, start, end, length):
    """The records of `src` with start <= POS - 1 < end, shifted to the interval, as a bgzipped indexed VCF whose
    contig chr20 is `length` long."""
    with gzip.open(src, "rt") as f, open(dst[:-3], "w") as out:
        for line in f:
            if line.startswith("##contig"):
                out.write("##contig=<ID=chr20,length=%d>\n" % length)
            elif line.startswith("#"):
                out.write(line)
            else:
                w = line.split("\t", 2)
                pos = int(w[1]) - 1
                if start <= pos and pos + 1 < end:
                    out.write("%s\t%d\t%s" % (w[0], pos - start + 1, w[2]))
    subprocess.run(["bgzip", "-f", dst[:-3]], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", dst], check=True)


def bed_slice(src, dst, start, end):
    with open(src) as f, open(dst, "w") as out:
        for line in f:
            w = line.split("\t")
            s, e = max(int(w[1]), start), min(int(w[2]), end)
            if w[0] == "chr20" and s < e:
                out.write("chr20\t%d\t%d\n" % (s - start, e - start))


def prep(args):
    work, packages, regions, out = args[:4]
    arms = args[4:]
    sys.path.insert(0, HERE)
    import region_reads  # noqa: E402

    ids = region_ids(regions)
    pkgs = packages_of(packages, set(ids))
    os.makedirs(out, exist_ok=True)
    saved = os.path.join(out, "reads_by_region.json")
    reads = collections.defaultdict(set)
    if os.path.exists(saved):
        # The arms' graphs and alignments read once: the read-to-region map is kept.
        for read, rids in json.load(open(saved)).items():
            reads[read] = set(rids)
    else:
        vg = os.environ["VG"]
        walks = [region_reads.node_lengths_and_chm13(vg, os.path.join(work, "arm_%s" % a, "sampled.gbz"))
                 for a in arms]
        intervals = region_reads.region_intervals(packages, walks[0], walks)
        wanted_ids = set(pkgs)
        # The reads the arms map within FLANK letters of a region, as well as into it: the local graph has those flanks.
        flanks = region_reads.flank_nodes(packages, walks[0], walks, FLANK, wanted_ids)
        for a, arm in enumerate(arms):
            mine = [x for x in intervals[a] if x[2] in wanted_ids]
            for read, rids in region_reads.read_regions(os.path.join(work, "arm_%s" % arm, "short.gaf.gz"), mine,
                                                        every=True, extra=flanks[a]).items():
                reads[read.rsplit("/", 1)[0]].update(rids)
        with open(saved, "w") as f:
            json.dump({r: sorted(x) for r, x in reads.items()}, f)
    print("reads in the %d regions: %d" % (len(pkgs), len(reads)), flush=True)
    chrom = "".join(line.strip() for line in open(os.path.join(work, "calls", "chr20.fa")) if not line.startswith(">"))
    for rid, p in pkgs.items():
        d = os.path.join(out, rid)
        os.makedirs(d, exist_ok=True)
        start, end = max(0, int(p["ref_start"]) - FLANK), min(len(chrom), int(p["ref_end"]) + FLANK)
        meta = {"id": rid, "start": start, "end": end, "left": chrom[start : int(p["ref_start"])].upper(),
                "right": chrom[int(p["ref_end"]) : end].upper()}
        with open(os.path.join(d, "region.json"), "w") as f:
            json.dump(meta, f)
        with open(os.path.join(d, "chr20.fa"), "w") as f:
            f.write(">chr20\n%s\n" % chrom[start:end].upper())
        with open(os.path.join(d, "chr20.fa.fai"), "w") as f:
            f.write("chr20\t%d\t7\t%d\t%d\n" % (end - start, end - start, end - start + 1))
        truth = os.path.join(work, "calls", "truth")
        for kind in ("smvar", "stvar"):
            vcf_slice(os.path.join(truth, "truth.chr20.%s.vcf.gz" % kind),
                      os.path.join(d, "truth.chr20.%s.vcf.gz" % kind), start, end, end - start)
            bed_slice(os.path.join(os.environ["TRUTH"], "CHM13v2.0_HG2-T2TQ100-V1.1_%s.benchmark.bed" % kind),
                      os.path.join(d, "truth.chr20.%s.bed" % kind), start, end)
    # Each region's reads gathered in memory and written once (a file per region and FASTQ open at a time would
    # pass the open-file limit).
    for name, fq in (("R1", "hg002.chr20.R1.fq.gz"), ("R2", "hg002.chr20.R2.fq.gz"), ("single", "hg002.chr20.single.fq.gz")):
        kept = collections.defaultdict(list)
        for read, text in fastq_records(os.path.join(work, fq)):
            for rid in reads.get(read, ()):
                kept[rid].append(text)
        for rid in pkgs:
            with open(os.path.join(out, rid, name + ".fq"), "w") as f:
                f.writelines(kept.get(rid, ()))
    print("prepared %d regions in %s" % (len(pkgs), out))


def msa_rows(source, pkg, msas):
    """{allele id: row} of the region's alleles, from a msas.jsonl.gz record or plegma's run on the distinct alleles."""
    alleles = pkg["alleles"]
    if source.endswith(".jsonl.gz"):
        rec = msas.get(pkg["id"])
        if rec is None or rec.get("status") != "ok":
            return None
        rows = {a["id"]: rec["rows"][a["id"]] for a in alleles if a["id"] in rec["rows"]}
    else:
        # CHM13's allele first (plegma -A -Q anchors the MSA on the first sequence).
        ref_seq = next((a["seq"].upper() for a in alleles if a["id"] == pkg.get("ref")), None)
        distinct = list(dict.fromkeys(([ref_seq] if ref_seq else []) + [a["seq"].upper() for a in alleles]))
        # WEIGHTED=1: each distinct sequence written once per haplotype spelling it, so -A weighs it so.
        copies = collections.Counter()
        for a in alleles:
            copies[a["seq"].upper()] += a["weight"] if os.environ.get("WEIGHTED") == "1" else 0
        if len(distinct) < 2:
            row_of = {distinct[0]: distinct[0]} if distinct else {}
        else:
            with tempfile.TemporaryDirectory() as tmp:
                path = os.path.join(tmp, "in.fa")
                with open(path, "w") as f:
                    for k, s in enumerate(distinct):
                        f.write(">s%d\n%s\n" % (k, s))
                        for c in range(copies[s] - 1):
                            f.write(">s%d_%d\n%s\n" % (k, c, s))
                text = subprocess.run([os.environ["PLEGMA"]] + source.split() + ["-r", "1", path], capture_output=True,
                                      text=True, check=True).stdout
            got, name = {}, None
            for line in text.splitlines():
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    got[name] = []
                elif name is not None:
                    got[name].append(line.strip().upper())
            row_of = {s: "".join(got["s%d" % k]) for k, s in enumerate(distinct)}
        rows = {a["id"]: row_of[a["seq"].upper()] for a in alleles}
    # Columns no allele uses go (a fragment's).
    width = len(next(iter(rows.values()))) if rows else 0
    keep = [c for c in range(width) if any(r[c] != "-" for r in rows.values())]
    return {a: "".join(r[c] for c in keep) for a, r in rows.items()}


def run_one(task):
    label, source, pkg, msas, prep_dir, out_dir = task
    sys.path.insert(0, os.environ["VNTR_PGREALIGN"])
    from pgrealign import induce  # noqa: E402

    rid = pkg["id"]
    src = os.path.join(prep_dir, rid)
    meta = json.load(open(os.path.join(src, "region.json")))
    try:
        rows = msa_rows(source, pkg, msas)
    except subprocess.CalledProcessError as error:
        return label, rid, {"error": "plegma failed: %s" % (error.stderr or "")[-120:]}
    if rows is None:
        return label, rid, {"error": "no MSA"}
    ref = pkg.get("ref")
    left, right = meta["left"], meta["right"]
    flanked = {a: left + r + right for a, r in rows.items()}
    g = induce.induce(flanked)
    d = os.path.join(out_dir, "%s.%s" % (label, rid))
    c = os.path.join(d, "chr20")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(c)
    # The MSA's rows (unflanked), for event_presence.py.
    with open(os.path.join(d, "rows.json"), "w") as f:
        json.dump(rows, f)
    gfa = os.path.join(d, "local.gfa")
    weight = {a["id"]: a["weight"] for a in pkg["alleles"]}
    with open(gfa, "w") as f:
        f.write("H\tVN:Z:1.1\tRS:Z:CHM13\n")
        for n in sorted(g.seq):
            f.write("S\t%d\t%s\n" % (n, g.seq[n]))
        for u, v in sorted(g.edges):
            f.write("L\t%d\t+\t%d\t+\t0M\n" % (u, v))
        for a, walk in g.walks.items():
            if not walk:
                continue
            steps = "".join(">%d" % n for n in walk)
            length = sum(len(g.seq[n]) for n in walk)
            if a == ref:
                f.write("W\tCHM13\t0\tchr20\t0\t%d\t%s\n" % (length, steps))
            for h in range(min(weight.get(a, 1), MAX_COPIES)):
                f.write("W\t%s\t%d\tchr20\t0\t%d\t%s\n" % (a, h + 1, length, steps))
    vg = os.environ["VG"]
    prefix = os.path.join(d, "map")

    def run(cmd, out=None):
        with open(out, "w") if out else open(os.devnull, "w") as o:
            proc = subprocess.run(cmd, stdout=o, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            raise RuntimeError(" ".join(cmd[:3]) + ": " + proc.stderr[-300:])

    try:
        run([vg, "autoindex", "-w", "sr-giraffe", "-g", gfa, "-p", prefix, "-t", "1", "-T", d])
        gbz = prefix + ".giraffe.gbz"
        idx = ["-Z", gbz, "-d", prefix + ".dist", "-m", prefix + ".shortread.withzip.min", "-z",
               prefix + ".shortread.zipcodes", "-o", "gaf", "-t", "1"]
        paired = os.path.join(d, "paired.gaf")
        single = os.path.join(d, "single.gaf")
        if os.path.getsize(os.path.join(src, "R1.fq")) > 0:
            run([vg, "giraffe"] + idx + ["-f", os.path.join(src, "R1.fq"), "-f", os.path.join(src, "R2.fq"),
                                         "--fragment-mean", "401.7", "--fragment-stdev", "166.4"], paired)
        else:
            open(paired, "w").close()
        if os.path.getsize(os.path.join(src, "single.fq")) > 0:
            run([vg, "giraffe"] + idx + ["-f", os.path.join(src, "single.fq")], single)
        else:
            open(single, "w").close()
        gaf = os.path.join(d, "reads.gaf")
        with open(gaf, "w") as f:
            for part in (paired, single):
                for line in open(part):
                    if not line.startswith("@"):
                        f.write(line)
        gafbase = os.environ["GAFBASE"]
        sorted_gaf = os.path.join(d, "reads.sorted.gaf")
        run([gafbase, "sort", gaf, "-o", sorted_gaf])
        db = os.path.join(d, "reads.gaf.db")
        run([gafbase, "construct", sorted_gaf, "-r", gbz, "-o", db, "--overwrite"])
        vcf = os.path.join(c, "chr20.vcf")
        run([vg, "call", gbz, "-p", "CHM13#0#chr20", "-d", "2", "-t", "1", "-s", "HG002", "--mismap-max", "0.95",
             "--read-likelihood", "--phased", "--gaf-base", db, "--gbz-base", gbz], vcf)
        run(["bgzip", "-f", vcf])
        run(["tabix", "-f", "-p", "vcf", vcf + ".gz"])
    except RuntimeError as error:
        return label, rid, {"error": str(error)}
    for name in ("chr20.fa", "chr20.fa.fai", "truth.chr20.smvar.vcf.gz", "truth.chr20.smvar.vcf.gz.tbi",
                 "truth.chr20.stvar.vcf.gz", "truth.chr20.stvar.vcf.gz.tbi", "truth.chr20.smvar.bed",
                 "truth.chr20.stvar.bed"):
        os.symlink(os.path.join(src, name), os.path.join(c, name))
    vce = os.environ["VCE"]
    for sub in ("scripts/wgs", "scripts/tier2", "scripts", "src"):
        sys.path.insert(0, os.path.join(vce, sub))
    os.environ["PATH"] = os.pathsep.join(x for x in (os.environ.get("BCFTOOLS_DIR"), os.environ.get("AARDVARK_DIR"),
                                                     os.environ["PATH"]) if x)
    import bench_wgs  # noqa: E402
    import bench_metrics as bm  # noqa: E402

    truvari = os.environ.get("TRUVARI", os.path.join(vce, "work", "truvari-venv", "bin", "truvari"))
    r = bench_wgs.score_contig(Path(d), "chr20", "HG002", 1, truvari)
    a = bench_wgs.pick(r.get("aardvark"), "GT", "ALL")
    small = bm.aardvark_counts(a)
    sv = bm.truvari_counts(r.get("truvari"))
    out = {"sv": (sv.truth_tp, sv.truth_fn, sv.query_fp), "small": (small.truth_tp, small.truth_fn, small.query_fp),
           "nodes": len(g.seq)}
    for key in ("aardvark_error", "truvari_error", "error"):
        if key in r:
            out[key] = r[key]
    return label, rid, out


def run_all(args):
    jobs = 8
    if "--jobs" in args:
        jobs = int(args[args.index("--jobs") + 1])
        del args[args.index("--jobs") : args.index("--jobs") + 2]
    chosen = None
    if "--regions" in args:
        chosen = region_ids(args[args.index("--regions") + 1])
        del args[args.index("--regions") : args.index("--regions") + 2]
    prep_dir, packages, out_dir = args[:3]
    labels = [tuple(x.split("=", 1)) for x in args[3:]]
    ids = chosen or sorted(x for x in os.listdir(prep_dir) if os.path.exists(os.path.join(prep_dir, x, "region.json")))
    pkgs = packages_of(packages, set(ids))
    msas = {}
    for label, source in labels:
        if source.endswith(".jsonl.gz"):
            with gzip.open(source, "rt") as f:
                for line in f:
                    r = json.loads(line)
                    if r["id"] in pkgs:
                        msas.setdefault(label, {})[r["id"]] = r
    os.makedirs(out_dir, exist_ok=True)
    # vg call runs gbz-base, which sits beside gaf-base.
    os.environ["PATH"] = os.path.dirname(os.environ["GAFBASE"]) + os.pathsep + os.environ["PATH"]
    # Each task carries its own region's MSA record only (all of them per task cost minutes in pickling).
    tasks = [(label, source, pkgs[i], {i: msas[label][i]} if i in msas.get(label, {}) else {}, prep_dir, out_dir)
             for label, source in labels for i in ids if i in pkgs]
    totals = {label: collections.Counter() for label, _ in labels}
    found = {}
    with Pool(jobs) as pool:
        for label, rid, res in pool.imap_unordered(run_one, tasks):
            found[(label, rid)] = res
            if "sv" in res:
                for k, v in zip(("sv_tp", "sv_fn", "sv_fp"), res["sv"]):
                    totals[label][k] += v
                for k, v in zip(("sm_tp", "sm_fn", "sm_fp"), res["small"]):
                    totals[label][k] += v
            else:
                totals[label]["failed"] += 1
    with open(os.path.join(out_dir, "local_calls.tsv"), "w") as f:
        f.write("region\tlabel\tsv_tp\tsv_fn\tsv_fp\tsmall_tp\tsmall_fn\tsmall_fp\tnodes\terror\n")
        for rid in ids:
            line = []
            for label, _ in labels:
                res = found.get((label, rid), {})
                sv, sm = res.get("sv", ("", "", "")), res.get("small", ("", "", ""))
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ((rid, label) + tuple(sv) + tuple(sm) + (
                    res.get("nodes", ""), res.get("error", res.get("truvari_error", "")))))
                line.append("%s SV %s/%s/%s small %s/%s/%s" % ((label,) + tuple(sv) + tuple(sm)) if "sv" in res
                            else "%s %s" % (label, res.get("error", "?")[:60]))
            print("%-22s %s" % (rid, " | ".join(line)))
    for label, t in totals.items():
        p = t["sv_tp"] / max(t["sv_tp"] + t["sv_fp"], 1)
        r = t["sv_tp"] / max(t["sv_tp"] + t["sv_fn"], 1)
        sp = t["sm_tp"] / max(t["sm_tp"] + t["sm_fp"], 1)
        sr = t["sm_tp"] / max(t["sm_tp"] + t["sm_fn"], 1)
        print("%-12s structural TP %d FN %d FP %d (F1 %.4f) | small TP %d FN %d FP %d (F1 %.4f) | failed %d" % (
            label, t["sv_tp"], t["sv_fn"], t["sv_fp"], 2 * p * r / max(p + r, 1e-9), t["sm_tp"], t["sm_fn"],
            t["sm_fp"], 2 * sp * sr / max(sp + sr, 1e-9), t["failed"]))


def main():
    if sys.argv[1] == "prep":
        prep(sys.argv[2:])
    else:
        run_all(sys.argv[2:])


if __name__ == "__main__":
    main()
