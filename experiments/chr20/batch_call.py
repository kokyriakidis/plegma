#!/usr/bin/env python3
"""local_call.py's per-region genotyping, made fast: every region's graph indexed, mapped and called on its own in
parallel (seconds each), its calls cached by its MSA, so a rerun maps and calls only the regions whose MSA changed;
then all of them scored at once.

usage: python3 -I experiments/chr20/batch_call.py <prep dir> <packages.jsonl.gz> <out dir>
           <label>=<msas.jsonl.gz | plegma flags> [...] [--regions a,b|file] [--jobs n] [--seed-rows] [--sample n]

Prep is local_call.py prep's (each region's flanks, reads, CHM13 slice and truth). Per label:
1. each region's MSA (in parallel; as local_call.py makes it; plegma's cached by the binary, the flags and the alleles,
   and with `--seed-rows` taken from a previous run's <out dir>/<label>.<region>/rows.json when not cached yet),
   restricted (`--sample`, 32; 0 keeps all) to CHM13's row and the alleles the region's reads support best, as the
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
and spurious calls, as event_presence.py reads them. Needs what local_call.py run needs.
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


def read_kmers(prep_dir, rid):
    """The 29-mers of a region's reads, both strands."""
    comp = str.maketrans("ACGT", "TGCA")
    out = set()
    for name in ("R1", "R2", "single"):
        path = os.path.join(prep_dir, rid, name + ".fq")
        with open(path) as f:
            for k, line in enumerate(f):
                if k % 4 != 1:
                    continue
                seq = line.strip().upper()
                for s in (seq, seq.translate(comp)[::-1]):
                    for i in range(len(s) - K + 1):
                        out.add(s[i : i + K])
    return out


def sampled(rows, pkg, prep_dir, count):
    """The MSA restricted to CHM13's row and the `count` alleles the region's reads support best (an allele's distinct
    29-mers the reads hold less those they do not, then its haplotypes: an allele the reads fit whole outranks a
    shorter one whose fewer k-mers they also hold), as vg haplotypes samples the haplotypes the production graph keeps
    by the sample's k-mers; the columns no kept row uses dropped."""
    kmers = read_kmers(prep_dir, pkg["id"])
    weight = {a["id"]: a["weight"] for a in pkg["alleles"]}
    score = {}
    for a, r in rows.items():
        seq = r.replace("-", "")
        mine = {seq[i : i + K] for i in range(len(seq) - K + 1)}
        held = len(mine & kmers)
        score[a] = (held - (len(mine) - held), weight.get(a, 0))
    by_seq = {}
    for a in sorted(rows, key=lambda a: score[a], reverse=True):
        by_seq.setdefault(rows[a].replace("-", ""), a)
    keep = set(list(by_seq.values())[:count]) | {pkg.get("ref")}
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
    key = hashlib.sha1(json.dumps([rid, sorted(rows.items())]).encode()).hexdigest()[:20]
    weight = {a["id"]: a["weight"] for a in pkg["alleles"]}
    return rid, {"key": key, "rows": rows, "seq": g.seq, "edges": sorted(g.edges), "walks": g.walks,
                 "ref": pkg.get("ref"), "weight": weight}, None


def run(cmd, out=None, log=None):
    with open(out, "w") if out else open(os.devnull, "w") as o:
        proc = subprocess.run([str(c) for c in cmd], stdout=o, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(" ".join(str(c) for c in cmd[:3]) + ": " + proc.stderr[-400:])


def call_region(task):
    """One region's graph indexed, its reads mapped and its CHM13 path called, on a graph of its own (a call over a
    graph of many regions loads and decomposes all of them, however few it calls): (region, VCF header, records)."""
    rid, g, prep_dir, d = task
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
            for h in range(min(g["weight"].get(a, 1), MAX_COPIES)):
                f.write("W\t%s\t%d\t%s\t0\t%d\t%s\n" % (a, h + 1, rid, length, steps))
    src = os.path.join(prep_dir, rid)
    vg = os.environ["VG"]
    prefix = os.path.join(d, "map")
    try:
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
                proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
                out.writelines(x + "\n" for x in proc.stdout.splitlines() if x and not x.startswith("@"))
        gafbase = os.environ["GAFBASE"]
        run([gafbase, "sort", gaf, "-o", gaf + ".sorted"])
        db = os.path.join(d, "reads.gaf.db")
        run([gafbase, "construct", gaf + ".sorted", "-r", gbz, "-o", db, "--overwrite"])
        vcf = os.path.join(d, "calls.vcf")
        run([vg, "call", gbz, "-p", "CHM13#0#%s" % rid, "-d", "2", "-t", "1", "-s", "HG002", "--mismap-max", "0.95",
             "--read-likelihood", "--phased", "--gaf-base", db, "--gbz-base", gbz], vcf)
    except RuntimeError as error:
        return rid, None, str(error)
    header, body = [], []
    for line in open(vcf):
        if line.startswith("#"):
            if not line.startswith("##contig"):
                header.append(line)
        else:
            body.append(line)
    shutil.rmtree(d, ignore_errors=True)
    return rid, header, body


def main():
    args = sys.argv[1:]
    jobs = 8
    if "--jobs" in args:
        jobs = int(args[args.index("--jobs") + 1])
        del args[args.index("--jobs") : args.index("--jobs") + 2]
    sample = 32
    if "--sample" in args:
        sample = int(args[args.index("--sample") + 1])
        del args[args.index("--sample") : args.index("--sample") + 2]
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
    import bench_wgs  # noqa: E402
    import bench_metrics as bm  # noqa: E402

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
            tasks = [(rid, g, prep_dir, os.path.join(d, "work", rid)) for rid, g in todo.items()]
            with Pool(jobs) as pool:
                for rid, header, body in pool.imap_unordered(call_region, tasks):
                    if header is None:
                        failed[rid] = body
                        continue
                    if not os.path.exists(header_path):
                        with open(header_path, "w") as f:
                            f.writelines(header)
                    with open(os.path.join(cache, todo[rid]["key"] + ".vcf"), "w") as f:
                        f.writelines(body)
            times["map and call"] = time.time() - t
        order_ok = [rid for rid in ids if rid in graphs and rid not in failed]
        # Score every region together: one multi-contig reference, truth and calls.
        t = time.time()
        c = os.path.join(d, "batch_score", "batch")
        os.makedirs(c)
        order = order_ok
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
            header = open(header_path).readlines()
            f.writelines(header[:-1])
            f.write(contigs)
            f.write(header[-1])
            for rid in order:
                f.writelines(open(os.path.join(cache, graphs[rid]["key"] + ".vcf")).readlines())
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
                    for line in open(os.path.join(prep_dir, rid, "truth.chr20.%s.bed" % kind)):
                        out.write(rid + "\t" + line.split("\t", 1)[1])
        r = bench_wgs.score_contig(Path(os.path.join(d, "batch_score")), "batch", "HG002", jobs, truvari)
        t_score = time.time() - t
        a = bench_wgs.pick(r.get("aardvark"), "GT", "ALL")
        small = bm.aardvark_counts(a)
        sv = bm.truvari_counts(r.get("truvari"))
        per = collections.defaultdict(collections.Counter)
        tdir = os.path.join(d, "batch_score", "score", "batch.truvari")
        for part, key in (("tp-base", "tp"), ("fn", "fn"), ("fp", "fp")):
            path = os.path.join(tdir, part + ".vcf.gz")
            if os.path.exists(path):
                with gzip.open(path, "rt") as f:
                    for line in f:
                        if not line.startswith("#"):
                            per[line.split("\t", 1)[0]][key] += 1
        # Per region, the layout event_presence.py reads: <out dir>/<label>.<region>/rows.json and its share of
        # truvari's missed and spurious calls.
        split = collections.defaultdict(lambda: collections.defaultdict(list))
        for part in ("fn", "fp"):
            path = os.path.join(tdir, part + ".vcf.gz")
            if os.path.exists(path):
                with gzip.open(path, "rt") as f:
                    for line in f:
                        if not line.startswith("#"):
                            split[line.split("\t", 1)[0]][part].append(line)
        for rid in order:
            rd = os.path.join(out_dir, "%s.%s" % (label, rid))
            sd = os.path.join(rd, "score", "chr20.truvari")
            os.makedirs(sd, exist_ok=True)
            with open(os.path.join(rd, "rows.json"), "w") as f:
                json.dump(graphs[rid]["rows"], f)
            for part in ("fn", "fp"):
                with gzip.open(os.path.join(sd, part + ".vcf.gz"), "wt") as f:
                    f.writelines(split[rid][part])
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
