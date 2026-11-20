#!/usr/bin/env python3
"""A grid of simulated local assemblies (reads anchored at both ends, one or the other; overlaps between the two
sides from none to whole, tandem arrays across the join, one side alone, reads past the far anchor, another copy's
one-sided reads) through `plegma poa`, every output checked against the truth: each MSA row spells its read; sides
overlapping 200 letters or more are joined, without a duplication or a deletion at the join; one side alone no
longer than the min_cover-th longest read; the consensus's edit distance to the truth.

usage: python3 -I experiments/validate/anchored_grid.py <work dir> <plegma> <karan CLI> [--jobs n] [--only regex]

Prints one line per case (OK, or FAIL and why). Each case's reads are seeded by its name, so two binaries run into
different work directories see the same inputs: diff their outputs to compare them. It found the joins refused when
the sides overlapped wholly, and a 20 kb assembly's accuracy lost to unpinned one-sided reads."""
import os, random, re, subprocess, sys, zlib
from multiprocessing import Pool
PLEGMA = KARAN = None

def mutated(rng, s, rate):
    out = []
    for c in s:
        r = rng.random()
        if r < rate / 3: out.append(rng.choice("ACGT"))
        elif r < 2 * rate / 3: continue
        elif r < rate: out.append(c); out.append(rng.choice("ACGT"))
        else: out.append(c)
    return "".join(out)

def make(path, sc):
    rng = random.Random(zlib.crc32(os.path.basename(path).encode()))
    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    L, rate = sc["L"], sc["rate"]
    if sc.get("tr"):
        unit, copies = sc["tr"]
        u = rnd(unit); arr = "".join(mutated(rng, u, 0.01) for _ in range(copies))
        mid = sc.get("tr_at", L // 2) - len(arr) // 2
        truth = rnd(mid) + arr + rnd(L - len(arr) - mid)
    else:
        truth = rnd(L)
    L = len(truth)
    over = sc.get("over", 0)
    before, after = rnd(over), rnd(over)
    reads = []
    # Each one-sided read's extent along the truth from its anchor (the TRIM check's measure).
    extents = []
    for _ in range(sc.get("nb", 0)):
        reads.append(("both", mutated(rng, truth, rate)))
    # one-sided: left ends; right starts (letters), first read at the extreme
    ov = sc.get("ov")
    nl, nr = sc.get("nl", 0), sc.get("nr", 0)
    m = sc.get("mid", L // 2)
    for k in range(nl):
        if ov is not None:
            e = min(L, m + ov // 2) if k < sc.get("deep", 1) else rng.randint(max(1, m - L // 4), min(L, m + ov // 2))
        else:
            e = rng.randint(int(sc["lend"][0] * L), int(sc["lend"][1] * L))
        reads.append(("left", mutated(rng, truth[:e] + (after if e >= L else ""), rate)))
        extents.append(e)
    for k in range(nr):
        if ov is not None:
            s = max(0, m - (ov - ov // 2)) if k < sc.get("deep", 1) else rng.randint(max(0, m - (ov - ov // 2)), min(L - 1, m + L // 4))
        else:
            s = rng.randint(int(sc["rstart"][0] * L), int(sc["rstart"][1] * L))
        reads.append(("right", mutated(rng, (before if s <= 0 else "") + truth[s:], rate)))
        extents.append(L - s)
    if sc.get("nf"):
        cut = sc.get("fcut", L // 3)
        other = truth[:cut] + mutated(rng, truth[cut:], 0.15)
        for _ in range(sc["nf"]):
            e = rng.randint(int(0.6 * L), L)
            reads.append(("left", mutated(rng, other[:e], rate)))
            extents.append(e)
    open(path + ".truth.fa", "w").write(">truth\n%s\n" % truth)
    with open(path + ".fa", "w") as f:
        for k, (a, s) in enumerate(reads):
            f.write(">r%d anchor=%s\n%s\n" % (k, a, s))
    return truth, reads, extents

def fasta(p):
    out = []; n = None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"): out.append([l[1:], ""])
        elif out: out[-1][1] += l
    return out

def dist(tpath, seq, mode, tmp):
    open(tmp, "w").write(">c\n%s\n" % seq)
    o = subprocess.run([KARAN, "--distance", "--mode", mode, tpath, tmp], capture_output=True, text=True).stdout
    last = o.strip().splitlines()[-1].split("\t")
    return int(last[2])

def run(task):
    path, sc = task
    truth, reads, extents = make(path, sc)
    v = []
    env = dict(os.environ, PLEGMA_TRACE="poa")
    extra = sc.get("args", [])
    p = subprocess.run([PLEGMA, "poa", "-r", "2"] + extra + [path + ".fa"], capture_output=True, text=True, env=env, timeout=900)
    open(path + ".out", "w").write(p.stdout); open(path + ".trace", "w").write(p.stderr)
    if p.returncode != 0:
        return path, sc, ["CRASH rc=%d %s" % (p.returncode, p.stderr.strip().splitlines()[-1:])], {}
    recs = []; n = None
    for l in p.stdout.splitlines():
        if l.startswith(">"): recs.append([l[1:], ""])
        else: recs[-1][1] += l
    cons = [r for r in recs if r[0].startswith("Consensus")]
    rows = [r for r in recs if not r[0].startswith("Consensus")]
    info = {"ncons": len(cons), "lens": [len(c[1]) for c in cons], "L": len(truth)}
    # MSA rows spell reads
    if len(rows) != len(reads):
        v.append("rows %d != reads %d" % (len(rows), len(reads)))
    else:
        bad = [k for k in range(len(reads)) if rows[k][1].replace("-", "").upper() != reads[k][1].upper()]
        if bad: v.append("MSA rows not spelling reads %s" % bad[:5])
    nb = sum(1 for a, _ in reads if a == "both")
    nl = sum(1 for a, _ in reads if a == "left"); nr = sum(1 for a, _ in reads if a == "right")
    tr = p.stderr
    m = re.search(r"joined, (\d+) letters", tr)
    info["joined"] = int(m.group(1)) if m else None
    jm = re.search(r"(join at [^\n]*)", tr); info["how"] = jm.group(1) if jm else ""
    if nb == 0 and nl and nr:
        ov = sc.get("true_ov")
        if ov is not None and ov >= 200 and len(cons) != 1:
            v.append("NOT JOINED (true overlap %d)" % ov)
    if len(cons) == 1:
        c = cons[0][1]
        mode = "global"
        if nb == 0 and nl and not nr: mode = "prefix"
        if nb == 0 and nr and not nl: mode = "suffix"
        d = dist(path + ".truth.fa", c, mode, path + ".c.fa")
        info["dist"] = d
        if mode == "global" and abs(len(c) - len(truth)) > max(20, 0.01 * len(truth)):
            v.append("LENGTH %d vs truth %d" % (len(c), len(truth)))
        if d > max(20, 0.01 * len(truth)):
            v.append("DIST %d" % d)
        if mode != "global":
            mc = 5
            if "-c" in extra: mc = int(extra[extra.index("-c") + 1])
            # The cut: where fewer than -c reads go on, along the truth (a read's own length stands for its reach
            # only up to its insertions and deletions). A few letters of slack for the consensus's own errors.
            lens = sorted(extents, reverse=True)
            cap = lens[min(mc, len(lens)) - 1]
            info["cap"] = cap
            slack = max(3, len(c) // 500)
            if len(c) > cap + slack:
                v.append("TRIM len %d > %d-th furthest read's reach %d" % (len(c), mc, cap))
            elif len(c) < cap - slack:
                v.append("TRIM len %d < %d-th furthest read's reach %d" % (len(c), mc, cap))
    if sc.get("nf"):
        info["foreign"] = re.search(r"another copy's[^\n]*", tr).group(0) if "another copy" in tr else "none flagged"
    return path, sc, v, info

def init(a, b):
    global PLEGMA, KARAN
    PLEGMA, KARAN = a, b

def grid():
    cases = []
    for L in (1000, 5000, 20000):
        for rate in (0.01, 0.05, 0.10):
            for seed in (1, 2):
                # no spanning, controlled overlap
                for ov in (0, 100, 199, 200, 500, 1900, 2100, 5000):
                    if ov >= L * 0.9: continue
                    cases.append(("ov%d" % ov, dict(L=L, rate=rate, nl=8, nr=8, ov=ov, deep=2, true_ov=ov)))
                # full overlap (every one-sided read reaches far anchor), w/ and w/o overhang
                for over in (0, 300):
                    cases.append(("full_o%d" % over, dict(L=L, rate=rate, nl=8, nr=8, lend=(1.0, 1.0), rstart=(0.0, 0.0), over=over, true_ov=L)))
                # tandem repeat at the join
                if L >= 5000:
                    for ov in (500, 2100):
                        cases.append(("tr_ov%d" % ov, dict(L=L, rate=rate, nl=8, nr=8, ov=ov, deep=2, tr=(40, 25), true_ov=ov)))
                # telomere-like
                for n in (3, 8, 15):
                    cases.append(("lonly%d" % n, dict(L=L, rate=rate, nl=n, lend=(0.3, 0.95))))
                    cases.append(("ronly%d" % n, dict(L=L, rate=rate, nr=n, rstart=(0.05, 0.7))))
                cases.append(("lonly_c1", dict(L=L, rate=rate, nl=6, lend=(0.3, 0.95), args=["-c", "1"])))
                # spanning + one-sided + overhang + foreign
                cases.append(("span_oh", dict(L=L, rate=rate, nb=3, nl=6, nr=6, lend=(1.0, 1.0), rstart=(0.0, 0.0), over=int(0.1 * L) + 30)))
                cases.append(("span_foreign", dict(L=L, rate=rate, nb=4, nl=4, nr=4, lend=(0.3, 0.95), rstart=(0.05, 0.7), nf=4)))
                cases.append(("few_both", dict(L=L, rate=rate, nb=1, nl=8, nr=8, lend=(0.3, 0.95), rstart=(0.05, 0.7))))
    out = []
    for name, sc in cases:
        sc = dict(sc)
        out.append(("%s_L%d_r%02d_s%d" % (name, sc["L"], int(sc["rate"] * 100), 0), sc))
    return out

def main():
    a = sys.argv[1:]
    jobs = int(a[a.index("--jobs") + 1]) if "--jobs" in a else 8
    only = a[a.index("--only") + 1] if "--only" in a else ""
    work, pl, ka = a[:3]
    os.makedirs(work, exist_ok=True)
    tasks = []
    seen = {}
    for name, sc in grid():
        k = seen.get(name, 0); seen[name] = k + 1
        name2 = name[:-1] + str(k + 1)
        if only and not re.search(only, name2): continue
        tasks.append((os.path.join(work, name2), sc))
    with Pool(jobs, initializer=init, initargs=(pl, ka)) as pool:
        for path, sc, v, info in pool.imap_unordered(run, tasks):
            print("%s\t%s\t%s" % (os.path.basename(path), "OK" if not v else "FAIL: " + "; ".join(v), info), flush=True)

if __name__ == "__main__":
    main()
