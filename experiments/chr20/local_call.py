#!/usr/bin/env python3
"""The local harness's prep and MSA helpers (evaluate.py runs them): each region on its own, its window of CHM13
and its truth, and its reads, so an MSA's graph can be genotyped in seconds per region instead of the chr20 arms' hours.

  prep:  python3 -I experiments/chr20/evaluate.py prep ... (or local_call.py prep <chr20 work dir> <packages.jsonl.gz>
             <regions a,b|file> <prep dir> <arm> [<arm> ...])

`prep`, once for a set of regions: the reads that the chr20 arms (`arm_<arm>`) align into each region or its flanks
(its nodes and CHM13's within READ_FLANK letters of them, as region_reads.py finds them), both mates, written per
region from the arms' FASTQs; CHM13's sequence over the region and FLANK letters each side; and the truth over that
interval (the T2T-Q100 small and structural VCFs and benchmark beds that score_calls.py sliced into
<work dir>/calls/truth), shifted to the interval's coordinates. `msa_rows` makes a region's MSA (plegma's, or a
pgrealign msas.jsonl.gz's) for batch_call.py.

Needs VG and TRUTH (prep); VNTR_PGREALIGN and, for flags, PLEGMA (msa_rows).
"""
import collections
import gzip
import json
import os
import subprocess
import sys
import tempfile

FLANK = 1500
READ_FLANK = 1000
"""The window's flanks are FLANK letters, the reads taken those mapped within READ_FLANK of the region: a read (and
most mates) then lies inside the window, where one reaching past its end would be forced through it."""
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
        # The reads the arms map within READ_FLANK letters of a region, as well as into it.
        flanks = region_reads.flank_nodes(packages, walks[0], walks, READ_FLANK, wanted_ids)
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


def msa_rows(source, pkg, msas, flanks=None):
    """{allele id: row} of the region's alleles, from a msas.jsonl.gz record or plegma's run on the distinct alleles.
    A bare `-F` among plegma's flags takes `flanks` (CHM13's left and right of the region) as its file."""
    alleles = pkg["alleles"]
    if source.endswith(".jsonl.gz"):
        rec = msas.get(pkg["id"])
        if rec is None or rec.get("status") != "ok":
            return None
        rows = {a["id"]: rec["rows"][a["id"]] for a in alleles if a["id"] in rec["rows"]}
    else:
        # CHM13's allele first (plegma ref anchors the MSA on the first sequence).
        ref_seq = next((a["seq"].upper() for a in alleles if a["id"] == pkg.get("ref")), None)
        distinct = list(dict.fromkeys(([ref_seq] if ref_seq else []) + [a["seq"].upper() for a in alleles]))
        if len(distinct) < 2:
            row_of = {distinct[0]: distinct[0]} if distinct else {}
        else:
            with tempfile.TemporaryDirectory() as tmp:
                path = os.path.join(tmp, "in.fa")
                flags = source.split()
                if "-F" in flags:
                    flanks_path = os.path.join(tmp, "flanks.fa")
                    with open(flanks_path, "w") as f:
                        f.write(">left\n%s\n>right\n%s\n" % flanks)
                    flags.insert(flags.index("-F") + 1, flanks_path)
                with open(path, "w") as f:
                    for k, s in enumerate(distinct):
                        f.write(">s%d\n%s\n" % (k, s))
                text = subprocess.run([os.environ["PLEGMA"]] + flags + ["-r", "1", path], capture_output=True,
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


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if args[0] == "prep":
        args = args[1:]
    prep(args)


if __name__ == "__main__":
    main()
