#!/usr/bin/env python3
"""HG002's two haplotypes of a region, from the truth and from each label's calls, set against the panel's alleles: did
the genotyper pick the alleles nearest the truth (then the MSA serves), or other ones (the graph misleads it)?

usage: python3 -I experiments/chr20/haps.py <batch_call out dir> <prep dir> <packages.jsonl.gz> <region> <label>...
       python3 -I experiments/chr20/haps.py <batch_call out dir> <prep dir> <packages.jsonl.gz> --regions <file>
           [--jobs n] <out.tsv> <label>...

For every label, the calls batch_call.py cached for the region's rows (<out>/cache/<key>.vcf) applied to the region's
CHM13 contig, haplotype by haplotype, and the truth's records the same; then the edit distances between the called and
the true haplotypes (paired the cheaper way), and, for each, the nearest allele of the label's (sampled) rows and of the
whole panel, with its distance. With --regions, a table of every region and label instead: the called haplotypes'
distance to the truth, the least the graph's alleles allow (the truth's nearest rows), and each called haplotype's
distance to its nearest row (above 0: a mosaic of rows, the genotyper's choices in neighbouring snarls disagreeing).
"""
import gzip
import hashlib
import json
import multiprocessing
import os
import sys


def distance(a, b):
    """Global edit distance by Myers's bit vectors (Hyyrö's form), `a` the pattern."""
    if not a:
        return len(b)
    m = len(a)
    peq = {}
    for i, c in enumerate(a):
        peq[c] = peq.get(c, 0) | (1 << i)
    mask = (1 << m) - 1
    top = 1 << (m - 1)
    pv, mv, score = mask, 0, m
    for c in b:
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & mask)
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        ph = ((ph << 1) | 1) & mask
        mh = (mh << 1) & mask
        pv = mh | (~(xv | ph) & mask)
        mv = ph & xv
    return score


def records(lines):
    out = []
    for line in lines:
        if line.startswith("#"):
            continue
        w = line.rstrip("\n").split("\t")
        gt = w[9].split(":")[0].replace("/", "|").split("|")
        out.append((int(w[1]), w[3], w[4].split(","), gt))
    return out


def apply(contig, recs, h):
    """Haplotype `h` of the contig under the phased records (a record overlapping one applied before skipped)."""
    out, at = [], 0
    for pos, ref, alts, gt in sorted(recs):
        if h >= len(gt) or gt[h] in (".", "0"):
            continue
        start = pos - 1
        if start < at:
            continue
        out.append(contig[at:start])
        out.append(alts[int(gt[h]) - 1])
        at = start + len(ref)
    out.append(contig[at:])
    return "".join(out)


def region(out_dir, prep_dir, pkg, rid):
    meta = json.load(open(os.path.join(prep_dir, rid, "region.json")))
    contig = "".join(x.strip() for x in open(os.path.join(prep_dir, rid, "chr20.fa")) if not x.startswith(">")).upper()
    truth = []
    for kind in ("smvar", "stvar"):
        with gzip.open(os.path.join(prep_dir, rid, "truth.chr20.%s.vcf.gz" % kind), "rt") as f:
            truth += records(f)
    return meta, contig, [apply(contig, truth, h) for h in (0, 1)]


def called(out_dir, rid, label, meta, contig, true_haps):
    """A label's called haplotypes paired with the truth's: (rows, haplotypes, the truth's index for each)."""
    rows = json.load(open(os.path.join(out_dir, "%s.%s" % (label, rid), "rows.json")))
    key = hashlib.sha1(json.dumps([rid, sorted(rows.items())]).encode()).hexdigest()[:20]
    calls = records(open(os.path.join(out_dir, "cache", key + ".vcf")))
    haps = [apply(contig, calls, h) for h in (0, 1)]
    straight = distance(haps[0], true_haps[0]) + distance(haps[1], true_haps[1])
    crossed = distance(haps[0], true_haps[1]) + distance(haps[1], true_haps[0])
    sampled = {}
    for a, r in rows.items():
        sampled.setdefault(meta["left"] + r.replace("-", "") + meta["right"], []).append(a)
    return calls, sampled, haps, [0, 1] if straight <= crossed else [1, 0]


def nearest(seq, pool, weight):
    return min((distance(seq, s), ids[0], sum(weight.get(i, 0) for i in ids)) for s, ids in pool.items())


def table_row(task):
    out_dir, prep_dir, pkg, labels = task
    rid = pkg["id"]
    weight = {a["id"]: a["weight"] for a in pkg["alleles"]}
    out = []
    try:
        meta, contig, true_haps = region(out_dir, prep_dir, pkg, rid)
        for label in labels:
            _, sampled, haps, pair = called(out_dir, rid, label, meta, contig, true_haps)
            got = [distance(haps[h], true_haps[pair[h]]) for h in (0, 1)]
            best = [nearest(t, sampled, weight)[0] for t in true_haps]
            mosaic = [nearest(x, sampled, weight)[0] for x in haps]
            out.append("\t".join(map(str, [rid, label, sum(got), sum(best), max(mosaic), min(mosaic),
                                            len(true_haps[0]) + len(true_haps[1]) - 2 * len(contig)])))
    except (OSError, ValueError) as error:
        out.append("%s\terror\t%s" % (rid, error))
    return out


def table(out_dir, prep_dir, packages, args):
    regions = set(open(args[args.index("--regions") + 1]).read().split())
    jobs = int(args[args.index("--jobs") + 1]) if "--jobs" in args else 8
    rest = [a for i, a in enumerate(args) if a not in ("--regions", "--jobs") and (i == 0 or args[i - 1] not in (
        "--regions", "--jobs"))]
    out_path, labels = rest[0], rest[1:]
    pkgs = []
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] in regions:
                pkgs.append(p)
    with multiprocessing.Pool(jobs) as pool, open(out_path, "w") as out:
        out.write("region\tlabel\tcalled_vs_truth\tleast\tmosaic_max\tmosaic_min\ttruth_length_change\n")
        tasks = [(out_dir, prep_dir, p, labels) for p in sorted(pkgs, key=lambda p: -len(p["alleles"]))]
        for lines in pool.imap_unordered(table_row, tasks):
            out.writelines(x + "\n" for x in lines)
            out.flush()
    summarize(out_path, labels)


def summarize(path, labels):
    rows = [x.rstrip("\n").split("\t") for x in open(path)][1:]
    rows = [r for r in rows if r[1] != "error"]
    by = {}
    for r in rows:
        by.setdefault(r[0], {})[r[1]] = r
    shared = [rid for rid in by if all(label in by[rid] for label in labels)]
    print("%d regions" % len(shared))
    for label in labels:
        got = sum(int(by[r][label][2]) for r in shared)
        least = sum(int(by[r][label][3]) for r in shared)
        exact = sum(1 for r in shared if by[r][label][2] == by[r][label][3])
        mosaics = sum(1 for r in shared if int(by[r][label][4]) > 10)
        print("  %-13s called vs truth %7d  least %6d  regions at the least %3d  with a mosaic (>10) %3d" % (
            label, got, least, exact, mosaics))


def main():
    out_dir, prep_dir, packages = sys.argv[1:4]
    if "--regions" in sys.argv:
        table(out_dir, prep_dir, packages, sys.argv[4:])
        return
    rid = sys.argv[4]
    labels = sys.argv[5:]
    pkg = None
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] == rid:
                pkg = p
                break
    meta, contig, true_haps = region(out_dir, prep_dir, pkg, rid)
    panel = {}
    for a in pkg["alleles"]:
        panel.setdefault(meta["left"] + a["seq"].upper() + meta["right"], []).append(a["id"])
    weight = {a["id"]: a["weight"] for a in pkg["alleles"]}
    print("%s: contig %d; truth haplotypes %d and %d letters" % (rid, len(contig), *map(len, true_haps)))
    for h, t in enumerate(true_haps):
        d, a, w = nearest(t, panel, weight)
        print("  truth %d: nearest panel allele %s (%d haplotypes) at %d" % (h, a, w, d))
    for label in labels:
        calls, sampled, haps, pair = called(out_dir, rid, label, meta, contig, true_haps)
        print("%s: %d calls; called vs truth %d" % (
            label, len(calls), sum(distance(haps[h], true_haps[pair[h]]) for h in (0, 1))))
        for h in (0, 1):
            t = true_haps[pair[h]]
            d, a, w = nearest(haps[h], sampled, weight)
            dt, at, _ = nearest(t, sampled, weight)
            print("  hap %d (%d letters) vs truth %d: %d; called is nearest %s (%d haplotypes) at %d; the truth's "
                  "nearest in the graph %s at %d" % (h, len(haps[h]), pair[h], distance(haps[h], t), a, w, d, at, dt))

if __name__ == "__main__":
    main()
