#!/usr/bin/env python3
"""plegma poa against shasta2's own local assembly (LocalAssembly7, with or without msaRepair and its homopolymer
model) on the local assemblies of a real shasta2 run, scored against a truth assembly.

usage: python3 -I experiments/msa/shasta2_la7.py <dump dir> <truth fasta> <minimap2 index> <plegma> <karan>
           [--minimap2 path] [--jobs n] [--min-mapq 10] [--limit n: an even sample of n]

The dump comes from a shasta2 run with SHASTA2_LA7_DUMP set (branch la7-dump of the shasta2 fork): per local
assembly between anchors A and B, `A_B.reads.fa` (every oriented read's sequence used, tagged anchor=both, left or
right), `A_B.consensus.fa` (shasta2's consensus, and before msaRepair when it ran) and `A_B.windows.fa` (each read on
both anchors 5000 bases past either anchor, the anchors' positions in it).

Truth: every window is mapped to the truth (minimap2 map-ont, primary alignments at --min-mapq or more) and its two
anchor positions lifted through the CIGAR; the truth interval (and strand) most of a local assembly's reads agree
on, held by at least two and half of them, gives its true sequence (reverse-complemented for reads on the reverse
strand). A local assembly whose reads do not agree is left out.

Contestants, each scored by its edit distance to the truth: shasta2's consensus as the run made it; shasta2's
consensus before msaRepair (where msaRepair ran); `plegma poa` on every read, anchored as tagged; `plegma poa` on the
reads on both anchors alone. Prints totals per contestant, per shasta2 method (FastPath, Abpoa, TheseusOnly,
TheseusAll) and per truth length, and writes every assembly's row to <dump dir>/../la7_eval.tsv.
"""
import collections
import os
import re
import subprocess
import sys
import tempfile
from multiprocessing import Pool

PLEGMA = KARAN = None


def records(path):
    out = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            out.append([line[1:], []])
        elif out:
            out[-1][1].append(line.strip())
    return [(h, "".join(s)) for h, s in out]


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))[::-1]


def lift(cigar, start, q):
    """The reference position of query position `q` (0-based, the query as aligned), None inside an insertion."""
    r, at = start, 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        consumes_q = op in "MIS=X"
        consumes_r = op in "MDN=X"
        if consumes_q and at <= q < at + n:
            return r + (q - at) if consumes_r else None
        if consumes_q:
            at += n
        if consumes_r:
            r += n
    return None


def fetch(fasta, fai, contig, begin, end):
    """`contig`'s letters begin ..< end, through the FASTA's index (no samtools: its library can be out of step)."""
    length, offset, line_bases, line_width = fai[contig]
    end = min(end, length)
    first = offset + (begin // line_bases) * line_width + begin % line_bases
    last = offset + ((end - 1) // line_bases) * line_width + (end - 1) % line_bases
    with open(fasta, "rb") as f:
        f.seek(first)
        raw = f.read(last - first + 1)
    return raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def option(args, name, default):
    if name in args:
        k = args.index(name)
        value = args[k + 1]
        del args[k:k + 2]
        return value
    return default


def score(task):
    """Every contestant's edit distance to the truth for one local assembly."""
    prefix, truth, method = task
    out = {"assembly": os.path.basename(prefix), "method": method, "truth_length": len(truth)}
    with tempfile.TemporaryDirectory() as tmp:
        truth_path = os.path.join(tmp, "truth.fa")
        with open(truth_path, "w") as f:
            f.write(">truth\n%s\n" % truth)
        cons = records(prefix + ".consensus.fa")
        candidates = {"shasta2": cons[0][1]}
        if len(cons) > 1:
            candidates["shasta2_before_repair"] = cons[1][1]
        reads = records(prefix + ".reads.fa")
        both = os.path.join(tmp, "both.fa")
        with open(both, "w") as f:
            for h, s in reads:
                if "anchor=both" in h:
                    f.write(">%s\n%s\n" % (h, s))
        for name, path in (("plegma", prefix + ".reads.fa"), ("plegma_both_only", both)):
            proc = subprocess.run([PLEGMA, "poa", "-r", "0", path], capture_output=True, text=True)
            got = records_text(proc.stdout)
            candidates[name] = got[0] if got else None
        for name, seq in candidates.items():
            if seq is None:
                out[name] = None
                continue
            path = os.path.join(tmp, name + ".fa")
            with open(path, "w") as f:
                f.write(">c\n%s\n" % seq)
            lines = subprocess.run([KARAN, "--distance", truth_path, path], capture_output=True,
                                   text=True).stdout.strip().splitlines()
            out[name] = int(lines[-1].split("\t")[2]) if len(lines) > 1 else None
    out["reads"] = sum(1 for _ in reads)
    out["both"] = sum(1 for h, _ in reads if "anchor=both" in h)
    return out


def records_text(text):
    seqs, cur = [], None
    for line in text.splitlines():
        if line.startswith(">"):
            cur = []
            seqs.append(cur)
        elif cur is not None:
            cur.append(line.strip())
    return ["".join(c) for c in seqs]


def main():
    global PLEGMA, KARAN
    args = sys.argv[1:]
    minimap2 = option(args, "--minimap2", "minimap2")
    jobs = int(option(args, "--jobs", "8"))
    min_mapq = int(option(args, "--min-mapq", "10"))
    limit = int(option(args, "--limit", "0"))
    dump, truth_fasta, index, PLEGMA, KARAN = args[:5]
    prefixes = sorted(os.path.join(dump, f[:-len(".reads.fa")]) for f in os.listdir(dump) if f.endswith(".reads.fa"))
    if limit and len(prefixes) > limit:
        # An even sample: every so many in anchor-id order.
        step = len(prefixes) / limit
        prefixes = [prefixes[int(k * step)] for k in range(limit)]
    work = os.path.join(os.path.dirname(os.path.abspath(dump)), "la7_eval_work")
    os.makedirs(work, exist_ok=True)
    # Every window, named assembly|read, mapped at once.
    windows = os.path.join(work, "windows.fa")
    anchors = {}
    with open(windows, "w") as f:
        for prefix in prefixes:
            name = os.path.basename(prefix)
            for header, seq in records(prefix + ".windows.fa"):
                fields = header.split()
                a = int(fields[1].split("=")[1])
                b = int(fields[2].split("=")[1])
                key = "%s|%s" % (name, fields[0])
                anchors[key] = (a, b, len(seq))
                f.write(">%s\n%s\n" % (key, seq))
    sam = os.path.join(work, "windows.sam")
    if not os.path.exists(sam) or os.path.getmtime(sam) < os.path.getmtime(windows):
        with open(sam, "w") as out:
            subprocess.run([minimap2, "-a", "--eqx", "-x", "map-ont", "-t", str(jobs), index, windows], stdout=out,
                           stderr=subprocess.DEVNULL, check=True)
    intervals = collections.defaultdict(list)
    for line in open(sam):
        if line.startswith("@"):
            continue
        w = line.split("\t")
        flag = int(w[1])
        if flag & 0x904 or int(w[4]) < min_mapq:
            continue
        key = w[0]
        a, b, length = anchors[key]
        reverse = bool(flag & 16)
        qa, qb = (length - 1 - a, length - 1 - b) if reverse else (a, b)
        start = int(w[3]) - 1
        ra, rb = lift(w[5], start, qa), lift(w[5], start, qb)
        if ra is None or rb is None:
            continue
        begin, end = (min(ra, rb) + 1, max(ra, rb) + 1) if reverse else (min(ra, rb), max(ra, rb))
        if begin < end:
            intervals[key.split("|")[0]].append((w[2], begin, end, reverse))
    fai = {}
    for line in open(truth_fasta + ".fai"):
        w = line.split("\t")
        fai[w[0]] = (int(w[1]), int(w[2]), int(w[3]), int(w[4]))
    tasks = []
    unresolved = 0
    for prefix in prefixes:
        name = os.path.basename(prefix)
        votes = collections.Counter(intervals.get(name, []))
        if not votes:
            unresolved += 1
            continue
        (contig, begin, end, reverse), n = votes.most_common(1)[0]
        if n < 2 or n * 2 < sum(votes.values()):
            unresolved += 1
            continue
        seq = fetch(truth_fasta, fai, contig, begin, end)
        if reverse:
            seq = revcomp(seq)
        method = records(prefix + ".consensus.fa")[0][0].split()[1].split("=")[1]
        tasks.append((prefix, seq, method))
    print("%d local assemblies, %d with an agreed truth, %d without" % (len(prefixes), len(tasks), unresolved))
    rows = []
    with Pool(jobs, initializer=init, initargs=(PLEGMA, KARAN)) as pool:
        for row in pool.imap_unordered(score, tasks, chunksize=4):
            rows.append(row)
    contestants = ["shasta2", "shasta2_before_repair", "plegma", "plegma_both_only"]
    with open(os.path.join(os.path.dirname(os.path.abspath(dump)), "la7_eval.tsv"), "w") as f:
        cols = ["assembly", "method", "truth_length", "reads", "both"] + contestants
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    def report(title, subset):
        print("\n%s (%d local assemblies, %d truth letters)" % (title, len(subset), sum(r["truth_length"] for r in subset)))
        for c in contestants:
            have = [r for r in subset if r.get(c) is not None]
            if not have:
                continue
            errors = sum(r[c] for r in have)
            exact = sum(1 for r in have if r[c] == 0)
            letters = sum(r["truth_length"] for r in have)
            print("  %-22s errors %7d (%.3f per kb)  exact %5d/%d" % (c, errors, 1000.0 * errors / max(1, letters),
                                                                     exact, len(have)))
        # Head to head, where both were scored.
        pairs = [("plegma", "shasta2"), ("plegma", "shasta2_before_repair")]
        for a, b in pairs:
            both = [r for r in subset if r.get(a) is not None and r.get(b) is not None]
            if both:
                better = sum(1 for r in both if r[a] < r[b])
                worse = sum(1 for r in both if r[a] > r[b])
                print("  %s vs %s: better in %d, worse in %d, equal in %d" % (a, b, better, worse,
                                                                             len(both) - better - worse))

    report("All", rows)
    for method in sorted({r["method"] for r in rows}):
        report("shasta2 method " + method, [r for r in rows if r["method"] == method])
    for lo, hi in ((0, 500), (500, 2000), (2000, 5000), (5000, 10 ** 9)):
        report("truth length %d-%s" % (lo, hi if hi < 10 ** 9 else "more"),
               [r for r in rows if lo <= r["truth_length"] < hi])


def init(plegma, karan):
    # Workers may not fork from this process (forkserver is Python 3.14's default), so set the globals there too.
    global PLEGMA, KARAN
    PLEGMA, KARAN = plegma, karan


if __name__ == "__main__":
    main()
