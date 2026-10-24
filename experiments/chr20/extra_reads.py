#!/usr/bin/env python3
"""HG002's reads that the chr20 FASTQs leave out (reads.sh takes GRCh38's chr20 by range: a read from sequence
GRCh38 lacks, or places on an alt, decoy or unplaced contig, or leaves unmapped, is not there), given to the regions
whose alleles they come from.

usage: python3 -I experiments/chr20/extra_reads.py <prep dir> <packages.jsonl.gz> <out prep dir> <genome.fa>
           <reads.fq.gz>...

<out prep dir> is a copy of <prep dir> whose regions' single.fq also hold every extra read sharing at least SHARE of its
K-mers with one region's alleles or flanks, as singles (their mates may be anywhere). Only K-mers that place a read in
the region count: those of fewer than three regions, and that the genome (CHM13, counted by jellyfish for these K-mers
alone) holds nowhere but in the region's own window: a satellite's copies on other chromosomes would otherwise pull
their reads into a chr20 region, where a whole-genome mapping would not put them. Prints per region the reads added.
"""
import collections
import gzip
import json
import os
import shutil
import subprocess
import sys

K = 31
SHARE = 0.3
STRIDE = 8
COMP = str.maketrans("ACGT", "TGCA")


def canonical(seq):
    rc = seq.translate(COMP)[::-1]
    n = len(seq)
    for i in range(n - K + 1):
        a = seq[i : i + K]
        b = rc[n - K - i : n - i]
        yield a if a <= b else b


def genome_counts(kmers, genome, work):
    """{K-mer: count in the genome} for `kmers` (canonical), by jellyfish counting only those."""
    jf = os.environ.get("JELLYFISH", os.path.expanduser("~/micromamba/envs/jellyfish/bin/jellyfish"))
    listed = work + ".kmers.fa"
    with open(listed, "w") as f:
        for i, km in enumerate(kmers):
            f.write(">%d\n%s\n" % (i, km))
    db = work + ".db"
    subprocess.run([jf, "count", "-m", str(K), "-C", "-s", str(2 * len(kmers)), "-t", "8", "--if", listed, "-o", db,
                    genome], check=True)
    out = {}
    dump = subprocess.run([jf, "dump", "-c", db], capture_output=True, text=True, check=True).stdout
    for line in dump.splitlines():
        km, n = line.split()
        out[km.upper()] = int(n)
    os.remove(listed)
    os.remove(db)
    return out


def main():
    prep_dir, packages, out_dir, genome = sys.argv[1:5]
    fastqs = sys.argv[5:]
    ids = sorted(x for x in os.listdir(prep_dir) if os.path.exists(os.path.join(prep_dir, x, "region.json")))
    wanted = set(ids)
    owners = collections.defaultdict(set)
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] not in wanted:
                continue
            meta = json.load(open(os.path.join(prep_dir, p["id"], "region.json")))
            for a in {x["seq"].upper() for x in p["alleles"]}:
                for km in canonical(meta["left"] + a + meta["right"]):
                    owners[km].add(p["id"])
    owners = {km: r for km, r in owners.items() if len(r) < 3}
    # The genome's counts of these K-mers, against the counts in the regions' own CHM13 windows.
    in_windows = collections.Counter()
    for rid in ids:
        window = "".join(x.strip() for x in open(os.path.join(prep_dir, rid, "chr20.fa")) if not x.startswith(">"))
        for km in canonical(window.upper()):
            if km in owners:
                in_windows[km] += 1
    in_genome = genome_counts(owners, genome, out_dir + ".jf")
    before = len(owners)
    owners = {km: next(iter(r)) if len(r) == 1 else r for km, r in owners.items()
              if in_genome.get(km, 0) <= in_windows.get(km, 0)}
    print("K-mers of %d regions: %d, %d of them the genome holds only in the regions" % (len(ids), before, len(owners)),
          flush=True)
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    shutil.copytree(prep_dir, out_dir, ignore=shutil.ignore_patterns("reads.kff*"))
    added = collections.defaultdict(list)
    seen = 0
    for path in fastqs:
        with gzip.open(path, "rt") as f:
            while True:
                head = f.readline()
                if not head:
                    break
                seq, plus, qual = f.readline(), f.readline(), f.readline()
                seen += 1
                s = seq.strip().upper()
                if len(s) < K:
                    continue
                # A screen on every STRIDE-th K-mer first: most reads share none with any region.
                probes = [s[i : i + K] for i in range(0, len(s) - K + 1, STRIDE)]
                comp = s.translate(COMP)[::-1]
                n = len(s)
                if not any(p in owners or comp[n - K - i : n - i] in owners
                           for i, p in zip(range(0, n - K + 1, STRIDE), probes)):
                    continue
                hits = collections.Counter()
                for km in canonical(s):
                    o = owners.get(km)
                    if o is None:
                        continue
                    if isinstance(o, str):
                        hits[o] += 1
                    else:
                        for r in o:
                            hits[r] += 1
                total = len(s) - K + 1
                for rid, n in hits.items():
                    if n >= SHARE * total:
                        added[rid].append(head + seq + plus + qual)
    for rid, recs in added.items():
        with open(os.path.join(out_dir, rid, "single.fq"), "a") as f:
            f.writelines(recs)
    print("extra reads read %d; added %d to %d regions" % (seen, sum(len(v) for v in added.values()), len(added)))
    for rid in sorted(added, key=lambda r: -len(added[r]))[:20]:
        print("  %-20s %d" % (rid, len(added[rid])))


if __name__ == "__main__":
    main()
