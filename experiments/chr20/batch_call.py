#!/usr/bin/env python3
"""local_call.py's per-region genotyping, made fast: every region's graph indexed, mapped and called on its own in
parallel (seconds each), its calls cached by its MSA, so a rerun maps and calls only the regions whose MSA changed;
then all of them scored at once.

usage: python3 -I experiments/chr20/batch_call.py <prep dir> <packages.jsonl.gz> <out dir>
           <label>=<msas.jsonl.gz | plegma flags> [...] [--regions a,b|file] [--jobs n] [--seed-rows] [--sample n]

Prep is local_call.py prep's (each region's flanks, reads, CHM13 slice and truth). Per label:
1. each region's MSA (in parallel; as local_call.py makes it; plegma's cached by the binary, the flags and the alleles,
   and with `--seed-rows` taken from a previous run's <out dir>/<label>.<region>/rows.json when not cached yet),
   restricted (`--sample`, 32; 0 keeps all; `--vg-haplotypes <kff | region>`: vg haplotypes itself samples, by those
   counts or the region's own reads', see call_region)
   to CHM13's row and the alleles the region's reads support best, as the
   production pipeline's vg haplotypes samples 32 haplotypes by the sample's k-mers before mapping (a graph of every
   allele, hundreds of haplotypes, takes minutes to index and map and genotypes a harder graph than production's), flanked by CHM13's flanks and induced into a graph
   (pgrealign's induce.py), and a key: the hash of its rows;
2. each region whose key has no cached calls (<out dir>/cache/<key>.vcf), in parallel: its graph as a GFA (its
   CHM13 walk the contig `CHM13#0#<region>`, each allele a haplotype walk per haplotype spelling it, up to
   MAX_COPIES), `vg autoindex`, `vg giraffe` of its reads, GAF-Base and `vg call` with the arms' flags, the calls
   into the cache (one graph of all the regions made each call load and decompose all of them: minutes a call);
3. every region's calls (cached or new) scored together by vg-call-eval's score_contig over one multi-contig
   reference and truth (aardvark for small variants, truvari for 50 bp or more), and the structural true, missed and
   spurious counts split by region.
Prints per label the totals and the time of each step; writes <out dir>/<label>/batch_calls.tsv (per region), keeps
each label's batch under <out dir>/<label>/, and per region <out dir>/<label>.<region>/ with its rows and its missed
and spurious calls (what strata.py, msa_check.py and evaluate.py report read). Needs VG, GAFBASE, BCFTOOLS_DIR,
AARDVARK_DIR, VCE, VNTR_PGREALIGN and, for plegma flags, PLEGMA.
"""
import collections
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

MAX_COPIES = 8
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import local_call  # noqa: E402


def msa_key(source, pkg):
    """The key of a region's MSA made by plegma flags: the binary, the flags, the environment switches plegma reads and
    the region's alleles."""
    binary = os.environ.get("PLEGMA", "")
    stamp = ""
    if binary and os.path.exists(binary):
        st = os.stat(binary)
        stamp = "%d.%d" % (st.st_size, int(st.st_mtime))
    alleles = sorted((a["id"], a["seq"].upper(), a["weight"]) for a in pkg["alleles"])
    # The environment switches plegma reads (PLEGMA_QCOSTS, PLEGMA_VOTERS, ...) and WEIGHTED change the MSA too.
    switches = sorted((k, v) for k, v in os.environ.items() if (k.startswith("PLEGMA_") or k == "WEIGHTED"))
    text = json.dumps([binary, stamp, source, pkg.get("ref"), alleles, switches])
    return hashlib.sha1(text.encode()).hexdigest()[:20]


K = 29


def read_kmer_counts(prep_dir, rid):
    """The 29-mers of a region's reads, both strands, with their counts."""
    comp = str.maketrans("ACGT", "TGCA")
    out = collections.Counter()
    for name in ("R1", "R2", "single"):
        path = os.path.join(prep_dir, rid, name + ".fq")
        with open(path) as f:
            for k, line in enumerate(f):
                if k % 4 != 1:
                    continue
                seq = line.strip().upper()
                for s in (seq, seq.translate(comp)[::-1]):
                    for i in range(len(s) - K + 1):
                        out[s[i : i + K]] += 1
    return out


PRESENT_DISCOUNT = 0.9
HET_ADJUSTMENT = 0.05
ABSENT_SCORE = 0.8


def sampled(rows, pkg, prep_dir, count):
    """The MSA restricted to CHM13's row and `count` alleles chosen as vg haplotypes chooses the haplotypes the
    production graph keeps: each allele's distinct 29-mers classed by the reads' counts against one haplotype's depth
    (the median count of the K-mers every allele holds, halved) as absent, heterozygous or homozygous; greedily, the
    allele whose K-mers score most (present +1, absent -ABSENT_SCORE) joins, and the K-mers it holds then count less
    (a homozygous one PRESENT_DISCOUNT times, a heterozygous one HET_ADJUSTMENT: the second haplotype is the one that
    holds what the first does not). So the two alleles the reads support come first, and an allele whose only mark is
    a deletion's junction is not crowded out by near copies of the reference. The columns no kept row uses dropped."""
    import numpy

    counts = read_kmer_counts(prep_dir, pkg["id"])
    weight = {a["id"]: a["weight"] for a in pkg["alleles"]}
    by_seq = {}
    for a in sorted(rows, key=lambda a: -weight.get(a, 0)):
        by_seq.setdefault(rows[a].replace("-", ""), a)
    seqs = list(by_seq)
    index = {}
    holds = []
    for seq in seqs:
        mine = {seq[i : i + K] for i in range(len(seq) - K + 1)}
        holds.append([index.setdefault(km, len(index)) for km in mine])
    m = numpy.zeros((len(seqs), len(index)), dtype=numpy.float32)
    for i, cols in enumerate(holds):
        m[i, cols] = 1.0
    observed = numpy.zeros(len(index))
    for km, j in index.items():
        observed[j] = counts.get(km, 0)
    shared = m.min(axis=0) > 0
    pool = observed[shared] if shared.any() else observed[observed > 0]
    haploid = max(1.0, float(numpy.median(pool)) / 2) if len(pool) else 1.0
    absent = observed < 0.3 * haploid
    het = ~absent & (observed < 1.5 * haploid)
    score = numpy.where(absent, -ABSENT_SCORE, 1.0)
    chosen = []
    left = list(range(len(seqs)))
    while left and len(chosen) < count:
        totals = m[left] @ score
        best = left[int(numpy.argmax(totals))]
        chosen.append(best)
        left.remove(best)
        mine = m[best] > 0
        score = numpy.where(mine & ~absent & ~het, score * PRESENT_DISCOUNT, score)
        score = numpy.where(mine & het, numpy.minimum(score, HET_ADJUSTMENT), score)
    keep = {by_seq[seqs[i]] for i in chosen} | {pkg.get("ref")}
    kept = {a: r for a, r in rows.items() if a in keep}
    width = len(next(iter(kept.values())))
    cols = [c for c in range(width) if any(r[c] != "-" for r in kept.values())]
    return {a: "".join(r[c] for c in cols) for a, r in kept.items()}


def graph_of(task):
    """A region's MSA rows (cached by msa_key when plegma makes them; with `seed`, a previous run's rows.json of the
    label taken as the MSA when the cache has none), its flanked graph and the rows' key."""
    label, source, pkg, msas, prep_dir, cache, seed, sample = task
    sys.path.insert(0, os.environ["VNTR_PGREALIGN"])
    from pgrealign import induce  # noqa: E402

    rid = pkg["id"]
    meta = json.load(open(os.path.join(prep_dir, rid, "region.json")))
    rows = None
    cached = None
    if not source.endswith(".jsonl.gz"):
        cached = os.path.join(cache, "msa." + msa_key(source, pkg) + ".json")
        if os.path.exists(cached):
            rows = json.load(open(cached))
        elif seed and os.path.exists(seed):
            rows = json.load(open(seed))
    try:
        if rows is None:
            rows = local_call.msa_rows(source, pkg, msas)
    except subprocess.CalledProcessError as error:
        return rid, None, "plegma failed: %s" % (error.stderr or "")[-120:]
    if cached and rows is not None and not os.path.exists(cached):
        with open(cached, "w") as f:
            json.dump(rows, f)
    if rows is None:
        return rid, None, "no MSA"
    if sample > 0:
        rows = sampled(rows, pkg, prep_dir, sample)
    flanked = {a: meta["left"] + r + meta["right"] for a, r in rows.items()}
    g = induce.induce(flanked)
    # The steps' invariants, checked here and reported (qc.tsv): every row spells its allele, every walk its flanked
    # allele, CHM13's its window of chr20.
    alleles = {a["id"]: a["seq"].upper() for a in pkg["alleles"]}
    window = "".join(x.strip() for x in open(os.path.join(prep_dir, rid, "chr20.fa")) if not x.startswith(">")).upper()
    ref = pkg.get("ref")
    qc = {"alleles": len(alleles), "rows": len(rows),
          "rows_spell": sum(1 for a, r in rows.items() if r.replace("-", "").upper() == alleles.get(a)),
          "walks_spell": sum(1 for a, w in g.walks.items()
                             if "".join(g.seq[n] for n in w).upper() == flanked[a].replace("-", "").upper()),
          "chm13_window": int(ref in rows and (meta["left"] + alleles[ref] + meta["right"]).upper() == window),
          "nodes": len(g.seq), "edges": len(g.edges), "columns": len(next(iter(rows.values())))}
    # vg haplotypes samples in call_region: the key says so (the rows are then all of them).
    key = hashlib.sha1(json.dumps([rid, sorted(rows.items())] + (["vg haplotypes"] if sample < 0 else []))
                       .encode()).hexdigest()[:20]
    weight = {a["id"]: a["weight"] for a in pkg["alleles"]}
    return rid, {"key": key, "rows": rows, "seq": g.seq, "edges": sorted(g.edges), "walks": g.walks,
                 "ref": pkg.get("ref"), "weight": weight, "qc": qc}, None


def run(cmd, out=None, log=None):
    with open(out, "w") if out else open(os.devnull, "w") as o:
        proc = subprocess.run([str(c) for c in cmd], stdout=o, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        # The cause first (a C++ exception's what(), or vg's error line), then the tail: a crash's stack trace
        # otherwise buries it.
        cause = [x.strip() for x in proc.stderr.splitlines() if "what():" in x or x.lstrip().startswith("error")]
        raise RuntimeError(" ".join(str(c) for c in cmd[:3]) + ": " + " | ".join(cause[:2]) + " | " +
                           proc.stderr[-300:])


def call_region(task):
    """One region's graph indexed, its reads mapped and its CHM13 path called, on a graph of its own (a call over a
    graph of many regions loads and decomposes all of them, however few it calls): (region, VCF header, records, the
    alleles vg haplotypes kept or None). With `kff`, every allele a haplotype walk once per haplotype spelling it, as
    the panel holds them, and the graph mapped to is vg haplotypes' sample of 32 by the K-mer counts, as the chr20
    arms make theirs (test_arm.sh: autoindex -w sampling, then --include-reference --set-reference CHM13)."""
    rid, g, prep_dir, d, kff = task
    if kff == "region":
        kff = region_kff(prep_dir, rid)
    os.makedirs(d, exist_ok=True)
    gfa = os.path.join(d, "local.gfa")
    with open(gfa, "w") as f:
        f.write("H\tVN:Z:1.1\tRS:Z:CHM13\n")
        for n in sorted(g["seq"]):
            f.write("S\t%d\t%s\n" % (n, g["seq"][n]))
        for u, v in g["edges"]:
            f.write("L\t%d\t+\t%d\t+\t0M\n" % (u, v))
        for a, walk in g["walks"].items():
            if not walk:
                continue
            steps = "".join(">%d" % n for n in walk)
            length = sum(len(g["seq"][n]) for n in walk)
            if a == g["ref"]:
                f.write("W\tCHM13\t0\t%s\t0\t%d\t%s\n" % (rid, length, steps))
            copies = g["weight"].get(a, 1) if kff else min(g["weight"].get(a, 1), MAX_COPIES)
            for h in range(copies):
                f.write("W\t%s\t%d\t%s\t0\t%d\t%s\n" % (a, h + 1, rid, length, steps))
    src = os.path.join(prep_dir, rid)
    vg = os.environ["VG"]
    prefix = os.path.join(d, "map")
    kept = None
    qc = {"sampling": "vg haplotypes" if kff else "none"}
    started = time.time()
    try:
        if kff:
            full = os.path.join(d, "full")
            run([vg, "autoindex", "-n", "-w", "sampling", "-g", gfa, "-p", full, "-t", "1", "-T", d])
            sampled_gbz = os.path.join(d, "sampled.gbz")
            run([vg, "haplotypes", "-t", "1", "-i", full + ".hapl", "-k", kff, "--include-reference", "--set-reference",
                 "CHM13", "--num-haplotypes", "32", "-g", sampled_gbz, full + ".gbz"])
            kept = sampled_alleles(vg, sampled_gbz, g)
            qc["kept_alleles"] = len(kept)
            qc["sampled_walks_unmatched"], qc["chm13_path_kept"] = sampled_unmatched(vg, sampled_gbz, g)
            # From a GBZ, autoindex writes the indexes only: the graph is sampled.gbz itself.
            run([vg, "autoindex", "-w", "sr-giraffe", "-G", sampled_gbz, "-p", prefix, "-t", "1", "-T", d])
            gbz = sampled_gbz
        else:
            run([vg, "autoindex", "-w", "sr-giraffe", "-g", gfa, "-p", prefix, "-t", "1", "-T", d])
            gbz = prefix + ".giraffe.gbz"
        idx = ["-Z", gbz, "-d", prefix + ".dist", "-m", prefix + ".shortread.withzip.min", "-z",
               prefix + ".shortread.zipcodes", "-o", "gaf", "-t", "1"]
        gaf = os.path.join(d, "reads.gaf")
        with open(gaf, "w") as out:
            for args in ([os.path.join(src, "R1.fq"), os.path.join(src, "R2.fq")], [os.path.join(src, "single.fq")]):
                if os.path.getsize(args[0]) == 0:
                    continue
                cmd = [vg, "giraffe"] + idx + sum((["-f", a] for a in args), [])
                if len(args) == 2:
                    cmd += ["--fragment-mean", "401.7", "--fragment-stdev", "166.4"]
                proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if proc.returncode != 0:
                    raise RuntimeError("giraffe: " + proc.stderr[-400:])
                out.writelines(x + "\n" for x in proc.stdout.splitlines() if x and not x.startswith("@"))
        gafbase = os.environ["GAFBASE"]
        run([gafbase, "sort", gaf, "-o", gaf + ".sorted"])
        db = os.path.join(d, "reads.gaf.db")
        run([gafbase, "construct", gaf + ".sorted", "-r", gbz, "-o", db, "--overwrite"])
        vcf = os.path.join(d, "calls.vcf")
        run([vg, "call", gbz, "-p", "CHM13#0#%s" % rid, "-d", "2", "-t", "1", "-s", "HG002", "--mismap-max", "0.95",
             "--read-likelihood", "--phased", "--gaf-base", db, "--gbz-base", gbz], vcf)
    except RuntimeError as error:
        return rid, None, str(error), None, None
    header, body = [], []
    for line in open(vcf):
        if line.startswith("#"):
            if not line.startswith("##contig"):
                header.append(line)
        else:
            body.append(line)
    qc.update(read_qc(src, gaf))
    qc.update(call_qc(body))
    qc["seconds"] = round(time.time() - started, 1)
    shutil.rmtree(d, ignore_errors=True)
    return rid, header, body, kept, qc


def read_qc(src, gaf):
    """The mapping's numbers: reads given (each mate one), alignments, mapped, at MAPQ 0, at MAPQ >= 30, mean
    identity of the mapped (GAF's matches over its block), and those under 0.9 (forced: a read from elsewhere, or
    reaching past the window)."""
    given = 0
    for name in ("R1", "R2", "single"):
        with open(os.path.join(src, name + ".fq")) as f:
            given += sum(1 for k, _ in enumerate(f) if k % 4 == 0)
    aligned = mapped = q0 = q30 = low = 0
    matches = block = 0
    for line in open(gaf):
        w = line.split("\t", 12)
        aligned += 1
        if w[5] == "*":
            continue
        mapped += 1
        q = int(w[11])
        q0 += q == 0
        q30 += q >= 30
        matches += int(w[9])
        block += int(w[10])
        low += int(w[9]) < 0.9 * int(w[10])
    return {"reads": given, "alignments": aligned, "mapped": mapped, "mapq0": q0, "mapq30": q30,
            "identity": round(matches / block, 4) if block else 0.0, "identity_below_0.9": low}


def call_qc(body):
    """The calls' numbers: records, records with an indel of 10 letters or more, their least and median GQ."""
    gq = []
    big = 0
    for line in body:
        w = line.rstrip("\n").split("\t")
        if abs(max(len(a) for a in w[4].split(",")) - len(w[3])) < 10:
            continue
        big += 1
        f = dict(zip(w[8].split(":"), w[9].split(":")))
        gq.append(int(f["GQ"]) if f.get("GQ", ".").lstrip("-").isdigit() else 0)
    gq.sort()
    return {"records": len(body), "indel_records": big, "gq_min": gq[0] if gq else "",
            "gq_median": gq[len(gq) // 2] if gq else ""}


def sampled_unmatched(vg, gbz, g):
    """The sampled haplotypes no allele spells whole (recombinants, which no row of the MSA is), and whether CHM13's
    path still spells its window (test_arm.sh checks the same of chr20)."""
    comp = str.maketrans("ACGT", "TGCA")
    proc = subprocess.run([vg, "convert", "-f", gbz], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    seq, walks = {}, []
    for line in proc.stdout.splitlines():
        w = line.split("\t")
        if w[0] == "S":
            seq[w[1]] = w[2]
        elif w[0] == "W":
            walks.append((w[1], w[6]))
    spelled = {"".join(g["seq"][n] for n in walk).upper() for walk in g["walks"].values()}
    window = "".join(g["seq"][n] for n in g["walks"][g["ref"]]).upper()
    unmatched = 0
    chm13 = 0
    for sample, steps in walks:
        text, k = [], 0
        while k < len(steps):
            j = k + 1
            while j < len(steps) and steps[j] not in "<>":
                j += 1
            s = seq[steps[k + 1 : j]]
            text.append(s if steps[k] == ">" else s.translate(comp)[::-1])
            k = j
        hap = "".join(text).upper()
        if steps.startswith("<"):
            hap = hap.translate(comp)[::-1]
        if sample == "CHM13":
            chm13 = int(hap == window)
        else:
            unmatched += hap not in spelled
    return unmatched, chm13


def region_kff(prep_dir, rid):
    """The region's reads' 29-mer counts as KFF (jellyfish, then jf2kff, as kmers.sh makes chr20's), made once: vg
    haplotypes counts only the K-mers unique in its graph, and in a region's own graph a K-mer found elsewhere on chr20
    looks unique, so chr20's counts would inflate it."""
    out = os.path.join(prep_dir, rid, "reads.kff")
    if os.path.exists(out):
        return out
    jf = os.environ.get("JELLYFISH", os.path.expanduser("~/micromamba/envs/jellyfish/bin/jellyfish"))
    j2k = os.environ.get("JF2KFF", os.path.expanduser("~/.local/bin/jf2kff"))
    db = out + ".jf"
    reads = [os.path.join(prep_dir, rid, n + ".fq") for n in ("R1", "R2", "single")]
    run([jf, "count", "-m", "29", "-C", "-s", "10M", "-t", "1", "-o", db] + reads)
    dump = subprocess.Popen([jf, "dump", "-c", "-L", "2", db], stdout=subprocess.PIPE)
    subprocess.run([j2k, out + ".part", "29"], stdin=dump.stdout, check=True)
    dump.wait()
    os.replace(out + ".part", out)
    os.remove(db)
    return out


def sampled_alleles(vg, gbz, g):
    """The alleles whose walks vg haplotypes kept in `gbz` (its haplotypes' sequences matched to the alleles'; a walk
    stored in reverse read back forward), the reference's first."""
    comp = str.maketrans("ACGT", "TGCA")
    proc = subprocess.run([vg, "convert", "-f", gbz], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    seq, walks = {}, []
    for line in proc.stdout.splitlines():
        w = line.split("\t")
        if w[0] == "S":
            seq[w[1]] = w[2]
        elif w[0] == "W" and w[1] != "CHM13":
            walks.append(w[6])
    allele_of = {}
    for a, walk in g["walks"].items():
        allele_of.setdefault("".join(g["seq"][n] for n in walk).upper(), a)
    kept = [g["ref"]]
    for steps in walks:
        text, k = [], 0
        while k < len(steps):
            j = k + 1
            while j < len(steps) and steps[j] not in "<>":
                j += 1
            s = seq[steps[k + 1 : j]]
            text.append(s if steps[k] == ">" else s.translate(comp)[::-1])
            k = j
        hap = "".join(text).upper()
        if steps.startswith("<"):
            hap = hap.translate(comp)[::-1]
        a = allele_of.get(hap)
        if a is not None and a not in kept:
            kept.append(a)
    return kept


def score(order, calls, header, prep_dir, d, jobs, truvari):
    """The calls of the regions `order` ({region: VCF record lines on its local contig}, under one VCF `header`)
    scored together, as vg-call-eval scores a contig: one multi-contig reference, truth and calls under <d>/batch.
    Returns score_contig's result, the structural and small-variant counts, per region its structural tp/fn/fp and
    its truvari fn and fp lines."""
    import bench_wgs  # noqa: E402
    import bench_metrics as bm  # noqa: E402

    c = os.path.join(d, "batch")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(c)
    lengths = {}
    with open(os.path.join(c, "batch.fa"), "w") as fa, open(os.path.join(c, "batch.fa.fai"), "w") as fai:
        at = 0
        for rid in order:
            seq = "".join(line.strip() for line in open(os.path.join(prep_dir, rid, "chr20.fa"))
                          if not line.startswith(">"))
            lengths[rid] = len(seq)
            fa.write(">%s\n%s\n" % (rid, seq))
            fai.write("%s\t%d\t%d\t%d\t%d\n" % (rid, len(seq), at + len(rid) + 2, len(seq), len(seq) + 1))
            at += len(rid) + 2 + len(seq) + 1
    contigs = "".join("##contig=<ID=%s,length=%d>\n" % (rid, lengths[rid]) for rid in order)
    with open(os.path.join(c, "batch.vcf"), "w") as f:
        f.writelines(header[:-1])
        f.write(contigs)
        f.write(header[-1])
        for rid in order:
            f.writelines(calls[rid])
    run(["bgzip", "-f", os.path.join(c, "batch.vcf")])
    run(["tabix", "-f", "-p", "vcf", os.path.join(c, "batch.vcf.gz")])
    for kind in ("smvar", "stvar"):
        records, head = [], None
        for rid in order:
            with gzip.open(os.path.join(prep_dir, rid, "truth.chr20.%s.vcf.gz" % kind), "rt") as f:
                lines = f.readlines()
            if head is None:
                head = [x for x in lines if x.startswith("#") and not x.startswith("##contig")]
            # Only records wholly inside the window: one reaching past its end has a reference allele the local
            # reference no longer spells, which bcftools norm refuses.
            for x in lines:
                if x.startswith("#"):
                    continue
                w = x.split("\t", 4)
                if int(w[1]) - 1 + len(w[3]) <= lengths[rid]:
                    records.append("%s\t%s" % (rid, x.split("\t", 1)[1]))
        path = os.path.join(c, "truth.batch.%s.vcf" % kind)
        with open(path, "w") as f:
            f.writelines(head[:-1])
            f.write(contigs)
            f.write(head[-1])
            f.writelines(records)
        run(["bgzip", "-f", path])
        run(["tabix", "-f", "-p", "vcf", path + ".gz"])
        with open(os.path.join(c, "truth.batch.%s.bed" % kind), "w") as out:
            for rid in order:
                # The benchmark only where the MSA aligns: the flanks are CHM13's in every graph, so HG002's
                # variants there are no graph's to call.
                meta = json.load(open(os.path.join(prep_dir, rid, "region.json")))
                lo, hi = len(meta["left"]), lengths[rid] - len(meta["right"])
                for line in open(os.path.join(prep_dir, rid, "truth.chr20.%s.bed" % kind)):
                    w = line.rstrip("\n").split("\t")
                    a, b = max(int(w[1]), lo), min(int(w[2]), hi)
                    if a < b:
                        out.write("%s\t%d\t%d\n" % (rid, a, b))
    r = bench_wgs.score_contig(Path(d), "batch", "HG002", jobs, truvari)
    small = bm.aardvark_counts(bench_wgs.pick(r.get("aardvark"), "GT", "ALL"))
    sv = bm.truvari_counts(r.get("truvari"))
    per = collections.defaultdict(collections.Counter)
    split = collections.defaultdict(lambda: collections.defaultdict(list))
    tdir = os.path.join(d, "score", "batch.truvari")
    for part, key in (("tp-base", "tp"), ("fn", "fn"), ("fp", "fp")):
        path = os.path.join(tdir, part + ".vcf.gz")
        if os.path.exists(path):
            with gzip.open(path, "rt") as f:
                for line in f:
                    if not line.startswith("#"):
                        per[line.split("\t", 1)[0]][key] += 1
                        if key != "tp":
                            split[line.split("\t", 1)[0]][key].append(line)
    # Small variants per region: aardvark marks each record (BD: TP, FN, FP) in its truth and query.
    adir = os.path.join(d, "score", "batch.aardvark")
    for part, wanted in (("truth.vcf.gz", ("TP", "FN")), ("query.vcf.gz", ("FP",))):
        path = os.path.join(adir, part)
        if not os.path.exists(path):
            continue
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                w = line.rstrip("\n").split("\t")
                fields = dict(zip(w[8].split(":"), w[9].split(":")))
                if fields.get("BD") in wanted:
                    per[w[0]]["sm_" + fields["BD"].lower()] += 1
    return r, sv, small, per, split


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    jobs = 8
    if "--jobs" in args:
        jobs = int(args[args.index("--jobs") + 1])
        del args[args.index("--jobs") : args.index("--jobs") + 2]
    sample = 32
    if "--sample" in args:
        sample = int(args[args.index("--sample") + 1])
        del args[args.index("--sample") : args.index("--sample") + 2]
    kff = None
    if "--vg-haplotypes" in args:
        kff = args[args.index("--vg-haplotypes") + 1]
        kff = kff if kff == "region" else os.path.abspath(kff)
        del args[args.index("--vg-haplotypes") : args.index("--vg-haplotypes") + 2]
        sample = -1
    seed_rows = "--seed-rows" in args
    args = [a for a in args if a != "--seed-rows"]
    chosen = None
    if "--regions" in args:
        chosen = local_call.region_ids(args[args.index("--regions") + 1])
        del args[args.index("--regions") : args.index("--regions") + 2]
    prep_dir, packages, out_dir = args[:3]
    labels = [tuple(x.split("=", 1)) for x in args[3:]]
    ids = chosen or sorted(x for x in os.listdir(prep_dir) if os.path.exists(os.path.join(prep_dir, x, "region.json")))
    pkgs = local_call.packages_of(packages, set(ids))
    ids = [i for i in ids if i in pkgs]
    cache = os.path.join(out_dir, "cache")
    os.makedirs(cache, exist_ok=True)
    os.environ["PATH"] = os.pathsep.join(x for x in (os.path.dirname(os.environ["GAFBASE"]),
                                                     os.environ.get("BCFTOOLS_DIR"), os.environ.get("AARDVARK_DIR"),
                                                     os.environ["PATH"]) if x)
    vce = os.environ["VCE"]
    for sub in ("scripts/wgs", "scripts/tier2", "scripts", "src"):
        sys.path.insert(0, os.path.join(vce, sub))
    truvari = os.environ.get("TRUVARI", os.path.join(vce, "work", "truvari-venv", "bin", "truvari"))
    for label, source in labels:
        started = time.time()
        msas = {}
        if source.endswith(".jsonl.gz"):
            with gzip.open(source, "rt") as f:
                for line in f:
                    r = json.loads(line)
                    if r["id"] in pkgs:
                        msas[r["id"]] = r
        tasks = [(label, source, pkgs[i], {i: msas[i]} if i in msas else {}, prep_dir, cache,
                  os.path.join(out_dir, "%s.%s" % (label, i), "rows.json") if seed_rows else None, sample) for i in ids]
        graphs, failed = {}, {}
        with Pool(jobs) as pool:
            for rid, g, error in pool.imap_unordered(graph_of, tasks):
                if g is None:
                    failed[rid] = error
                else:
                    graphs[rid] = g
        t_msa = time.time() - started
        todo = {rid: g for rid, g in graphs.items() if not os.path.exists(os.path.join(cache, g["key"] + ".vcf"))}
        d = os.path.join(out_dir, label)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        times = {}
        header_path = os.path.join(cache, "header.vcf")
        if todo:
            t = time.time()
            tasks = [(rid, g, prep_dir, os.path.join(d, "work", rid), kff) for rid, g in todo.items()]
            with Pool(jobs) as pool:
                for rid, header, body, kept, qc in pool.imap_unordered(call_region, tasks):
                    if header is None:
                        failed[rid] = body
                        continue
                    with open(os.path.join(cache, todo[rid]["key"] + ".qc.json"), "w") as f:
                        json.dump(qc, f)
                    if not os.path.exists(header_path):
                        with open(header_path, "w") as f:
                            f.writelines(header)
                    with open(os.path.join(cache, todo[rid]["key"] + ".vcf"), "w") as f:
                        f.writelines(body)
                    if kept is not None:
                        with open(os.path.join(cache, todo[rid]["key"] + ".kept.json"), "w") as f:
                            json.dump(kept, f)
            # Where vg haplotypes refuses a region's graph (its partitioner can call a small graph's one top-level
            # chain a loop), that region sampled by `sampled` instead, and counted.
            refused = [rid for rid, e in failed.items() if kff and "haplotype sampling cannot be used" in str(e)]
            if refused:
                again = {}
                tasks = [(label, source, pkgs[i], {i: msas[i]} if i in msas else {}, prep_dir, cache, None, 32)
                         for i in refused]
                with Pool(jobs) as pool:
                    for rid, g, error in pool.imap_unordered(graph_of, tasks):
                        if g is not None:
                            again[rid] = g
                            graphs[rid] = g
                tasks = [(rid, g, prep_dir, os.path.join(d, "work", rid), None) for rid, g in again.items()
                         if not os.path.exists(os.path.join(cache, g["key"] + ".vcf"))]
                with Pool(jobs) as pool:
                    for rid, header, body, kept, qc in pool.imap_unordered(call_region, tasks):
                        if header is not None:
                            with open(os.path.join(cache, again[rid]["key"] + ".vcf"), "w") as f:
                                f.writelines(body)
                            qc["sampling"] = "greedy (vg refused)"
                            with open(os.path.join(cache, again[rid]["key"] + ".qc.json"), "w") as f:
                                json.dump(qc, f)
                for rid in again:
                    if os.path.exists(os.path.join(cache, again[rid]["key"] + ".vcf")):
                        failed.pop(rid, None)
                print("   vg haplotypes refused %d regions, sampled by the K-mer greedy instead" % len(refused))
            times["map and call"] = time.time() - t
        order_ok = [rid for rid in ids if rid in graphs and rid not in failed]
        # Score every region together: one multi-contig reference, truth and calls.
        t = time.time()
        order = order_ok
        r, sv, small, per, split = score(order, {rid: open(os.path.join(cache, graphs[rid]["key"] + ".vcf")).readlines()
                                                 for rid in order}, open(header_path).readlines(), prep_dir,
                                         os.path.join(d, "batch_score"), jobs, truvari)
        t_score = time.time() - t
        # Per region, the layout evaluate.py report reads: <out dir>/<label>.<region>/rows.json and its share of
        # truvari's missed and spurious calls.
        for rid in order:
            rd = os.path.join(out_dir, "%s.%s" % (label, rid))
            sd = os.path.join(rd, "score", "chr20.truvari")
            os.makedirs(sd, exist_ok=True)
            rows = graphs[rid]["rows"]
            kept_path = os.path.join(cache, graphs[rid]["key"] + ".kept.json")
            if os.path.exists(kept_path):
                # The rows vg haplotypes kept: the graph vg call had.
                kept = set(json.load(open(kept_path)))
                rows = {a: r for a, r in rows.items() if a in kept}
                width = len(next(iter(rows.values())))
                cols = [c for c in range(width) if any(r[c] != "-" for r in rows.values())]
                rows = {a: "".join(r[c] for c in cols) for a, r in rows.items()}
            with open(os.path.join(rd, "rows.json"), "w") as f:
                json.dump(rows, f)
            with open(os.path.join(rd, "key.txt"), "w") as f:
                f.write(graphs[rid]["key"] + "\n")
            for part in ("fn", "fp"):
                with gzip.open(os.path.join(sd, part + ".vcf.gz"), "wt") as f:
                    f.writelines(split[rid][part])
        # Every step's numbers per region (graph_of's invariants, call_region's sampling, mapping and calls, the
        # scores), and the regions whose invariants fail.
        cols = ["alleles", "rows", "rows_spell", "walks_spell", "chm13_window", "columns", "nodes", "edges", "sampling",
                "kept_alleles", "sampled_walks_unmatched", "chm13_path_kept", "reads", "alignments", "mapped", "mapq0", "mapq30",
                "identity", "identity_below_0.9", "records", "indel_records", "gq_min", "gq_median", "seconds"]
        broken = []
        with open(os.path.join(d, "qc.tsv"), "w") as f:
            f.write("region\t" + "\t".join(cols) + "\tsv_tp\tsv_fn\tsv_fp\tsm_tp\tsm_fn\tsm_fp\terror\n")
            for rid in ids:
                q = dict(graphs[rid]["qc"]) if rid in graphs else {}
                qp = os.path.join(cache, graphs[rid]["key"] + ".qc.json") if rid in graphs else ""
                if qp and os.path.exists(qp):
                    q.update(json.load(open(qp)))
                p = per.get(rid, collections.Counter())
                f.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" % (
                    rid, "\t".join(str(q.get(c, "")) for c in cols), p["tp"], p["fn"], p["fp"], p["sm_tp"], p["sm_fn"],
                    p["sm_fp"], failed.get(rid, "")[:200]))
                if q and (q.get("rows_spell") != q.get("rows") or q.get("walks_spell") != q.get("rows")
                          or not q.get("chm13_window") or (q.get("reads") and not q.get("mapped"))
                          or q.get("chm13_path_kept") == 0):
                    broken.append(rid)
        if broken:
            print("   invariants failing in %d regions (qc.tsv): %s" % (len(broken), " ".join(broken[:10])))
        with open(os.path.join(d, "batch_calls.tsv"), "w") as f:
            f.write("region\tsv_tp\tsv_fn\tsv_fp\terror\n")
            for rid in ids:
                p = per.get(rid, collections.Counter())
                f.write("%s\t%d\t%d\t%d\t%s\n" % (rid, p["tp"], p["fn"], p["fp"], failed.get(rid, "")))
        print("%-12s structural TP %d FN %d FP %d (F1 %.4f) | small TP %d FN %d FP %d (F1 %.4f) | %d regions, %d new, "
              "%d failed | seconds: MSAs and graphs %.0f, index, map and call %.0f, score %.0f, all %.0f" % (
                  label, sv.truth_tp, sv.truth_fn, sv.query_fp, sv.f1, small.truth_tp, small.truth_fn, small.query_fp,
                  small.f1, len(order), len(todo), len(failed), t_msa, times.get("map and call", 0), t_score,
                  time.time() - started), flush=True)
        for key in ("aardvark_error", "truvari_error", "error"):
            if key in r:
                print("   %s: %s" % (key, str(r[key])[:300]))


if __name__ == "__main__":
    main()
