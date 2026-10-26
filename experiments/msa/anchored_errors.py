#!/usr/bin/env python3
"""Every consensus error of an anchored benchmark run (experiments/msa/anchored.py), profiled: its kind, region (the
tandem array or a flank), coverage, and whether fixing it alone lowers the reads' summed edit distance, each read
aligned as polishing aligns it (`missed`: polishing should have taken the fix; `tie`; `limited`: the reads prefer the
error, as the information they hold allows).

usage: python3 -I experiments/msa/anchored_errors.py <work dir> [tool (plegma)]
"""
import os, re, subprocess, sys
from multiprocessing import Pool
K = os.path.expanduser("~/karan/build/cli/karan")
work = os.path.expanduser(sys.argv[1]); tool = sys.argv[2] if len(sys.argv) > 2 else "plegma"
def records(p):
    out = []; cur = None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"): cur = [l[1:], []]; out.append(cur)
        elif cur: cur[1].append(l.strip())
    return [(h, "".join(s)) for h, s in out]
def dist(ref, q, mode, tag):
    a = "/tmp/ep_%s_r.fa" % tag; b = "/tmp/ep_%s_q.fa" % tag
    open(a, "w").write(">r\n%s\n" % ref); open(b, "w").write(">q\n%s\n" % q)
    o = subprocess.run([K, "--distance", "--mode", mode, a, b], capture_output=True, text=True).stdout.strip().splitlines()
    return int(o[-1].split("\t")[2])
def total(cons, reads, tag):
    t = 0
    for anchor, r in reads:
        if anchor == "both": t += dist(cons, r, "global", tag)
        elif anchor == "left": t += min(dist(cons, r, "prefix", tag), dist(cons, r, "ends:0,0,0,%d" % len(r), tag))
        else: t += min(dist(cons, r, "suffix", tag), dist(cons, r, "ends:0,0,%d,0" % len(r), tag))
    return t
def case(name):
    p = os.path.join(work, name)
    if not os.path.exists(p + "." + tool + ".out"): return name, []
    truth = records(p + ".truth.fa")[0][1]
    cons = records(p + "." + tool + ".out")
    if not cons: return name, []
    cons = cons[0][1]
    reads = []
    for h, s in records(p + ".fa"):
        m = re.search(r"anchor=(\w+)", h); reads.append((m.group(1) if m else "both", s))
    scen = name.split("_L")[0]
    mode = {"left_only": "prefix", "right_only": "suffix"}.get(scen, "global")
    a = "/tmp/ep_%s_t.fa" % name; b = "/tmp/ep_%s_c.fa" % name
    open(a, "w").write(">t\n%s\n" % truth); open(b, "w").write(">c\n%s\n" % cons)
    f = subprocess.run([K, "--mode", mode, a, b], capture_output=True, text=True).stdout.strip().splitlines()[-1].split("\t")
    i = int(f[5]); j = 0; errs = []
    for c, op in re.findall(r"(\d+)([=XID])", f[9]):
        c = int(c)
        if op == "=": i += c; j += c; continue
        errs.append((i, j, op, c))
        if op in "XD": i += c
        if op in "XI": j += c
    L = len(truth); arr0 = (L - 1200) // 2; arr1 = arr0 + 1200
    base = total(cons, reads, name)
    out = []
    for (ti, cj, op, c) in errs:
        if op == "D": fixed = cons[:cj] + truth[ti:ti + c] + cons[cj:]
        elif op == "I": fixed = cons[:cj] + cons[cj + c:]
        else: fixed = cons[:cj] + truth[ti:ti + c] + cons[cj + c:]
        delta = total(fixed, reads, name) - base
        # coverage: reads whose truth span holds ti (by length from their anchor)
        cov = 0
        for anchor, r in reads:
            if anchor == "both" or (anchor == "left" and len(r) > ti) or (anchor == "right" and len(r) > L - ti): cov += 1
        region = "array" if arr0 <= ti < arr1 else "flank"
        near_end = min(ti, L - ti) < 200
        out.append((scen, region, op, c, cov, delta, near_end, ti))
    return name, out
if __name__ == "__main__":
    names = sorted(n[:-3] for n in os.listdir(work) if n.endswith(".fa") and not n.endswith(".truth.fa") and ".pericles" not in n and "." not in n[:-3])
    allrows = []
    with Pool(8) as pool:
        for name, rows in pool.imap(case, names):
            allrows += rows
            for r in rows:
                print(name, "pos", r[7], r[1], r[2] + str(r[3]), "cov", r[4], "delta", r[5], "missed" if r[5] < 0 else ("tie" if r[5] == 0 else "limited"), "end" if r[6] else "")
    print("\nsummary:")
    from collections import Counter
    c = Counter((r[0], "missed" if r[5] < 0 else ("tie" if r[5] == 0 else "limited")) for r in allrows)
    for k, v in sorted(c.items()): print(k, v)
    c2 = Counter((r[1], "missed" if r[5] < 0 else ("tie" if r[5] == 0 else "limited")) for r in allrows)
    for k, v in sorted(c2.items()): print(k, v)
