#!/usr/bin/env python3
"""The brute-force seed census (build/anchors --brute) over vntr-eval's regions: per pair, whether plegma's anchors
(split_anchors) and the flanks alone force any excess over the optimum, the narrowest corridor around them (and
around the straight line, no seeds) that holds an optimal alignment, in periods, and how many of the MUMs and the
longest SMEM-or-MUM lie on an optimal alignment, flank ones apart from those in the array.

usage: python3 -I experiments/anchors/brute_panel.py <vntr-eval> <anchors binary> <placement binary> <out dir>
           [--pairs n] [--jobs n] [--max-cells millions] [--prices x,o,e]

Corridors: per anchors (and the straight line), the narrowest half-width B holding an optimal alignment, the B the
bound proves, and stopping rules that need no proof: the first B whose best equals 2B's (`2x`), and 2B's and 4B's
(`4x`), each with its excess and cells against the full DP's.

Writes each pair's census to <out dir>/<region>.<a>_<b>.txt (and reuses one already there) and <out dir>/summary.tsv.
"""
import os
import random
import re
import subprocess
import sys
from multiprocessing import Pool


def distinct(path):
    seqs = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    return list(dict.fromkeys("".join(s) for s in seqs))


def one(task):
    ve, anchors, placement, out, rid, pairs, max_cells, prices = task
    hap = os.path.join(ve, "regions", rid, "hap32.fa")
    tandem = subprocess.run([placement, "--tandem", hap], capture_output=True, text=True).stdout
    m = re.search(r"^layout period (\d+)", tandem, re.M)
    period = int(m.group(1)) if m else 0
    seqs = distinct(hap)
    rng = random.Random(rid)
    all_pairs = [(i, j) for i in range(len(seqs)) for j in range(i + 1, len(seqs))
                 if len(seqs[i]) * len(seqs[j]) <= max_cells]
    rng.shuffle(all_pairs)
    rows = []
    for i, j in all_pairs[:pairs]:
        a, b = (i, j) if len(seqs[i]) >= len(seqs[j]) else (j, i)
        path = os.path.join(out, "%s.%d_%d.txt" % (rid, a, b))
        if os.path.exists(path):
            txt = open(path).read()
        else:
            txt = subprocess.run([anchors, "--brute", hap, str(a), str(b), str(max(period, 1)), "--prices", prices],
                                 capture_output=True, text=True).stdout
            with open(path, "w") as f:
                f.write(txt)
        rows.append((rid, period, a, b, txt))
    return rows


def corridor(text):
    """A `_corridor` report: excess per half-width, the proving width and its cells (percent of the full DP)."""
    bands = [(int(x), int(e)) for x, e in re.findall(r"(\d+):(-?\d+)", text.split("|")[0])]
    out = {"bands": bands, "band": next((x for x, e in bands if e == 0), 9999)}
    m = re.search(r"proved at (-?\d+), cells (\d+)%", text)
    if m:
        out["proved"], out["proved_cells"] = int(m.group(1)), int(m.group(2))
    return out


def stop(bands, k):
    """The first width whose best equals the next k - 1 doublings' (the last width when none does)."""
    for i in range(len(bands) - k + 1):
        if all(bands[i][1] == bands[i + d][1] for d in range(1, k)):
            return i
    return len(bands) - 1


def parse(rid, period, a, b, txt):
    r = {"rid": rid, "period": period}
    m = re.search(r"^pair \S+ \S+ lengths (\d+) (\d+) optimum (\d+)", txt, re.M)
    if not m:
        return None
    r["n"], r["m"], r["opt"] = map(int, m.groups())
    for key, label in (("split", "plegma split_anchors"), ("flanks", "plegma flanks only"), ("mem", "all MEMs 0.0"),
                       ("mumsmem", "MUMs \\+ SMEMs 0.0"), ("mum", "MUMs 0.0")):
        m = re.search(r"^\s*%s\s*\| (\d+) \((\d+)\) \| (\d+) \| (\d+) of (\d+) \| (-?\d+) \| corridor (.*)$"
                      % label, txt, re.M)
        if m:
            r[key] = {"anchors": int(m.group(1)), "right": int(m.group(2)), "excess": int(m.group(6)),
                      "forced_right": int(m.group(4)), "forced": int(m.group(5))}
            r[key].update(corridor(m.group(7)))
    m = re.search(r"straight corridor \(no seeds\) \| (.*)$", txt, re.M)
    if m:
        r["straight"] = corridor(m.group(1))
    # MUMs: flank ones (in the first or last 5% of a) apart from the array's.
    mums = [(int(x[0]), int(x[1]), int(x[2])) for x in
            re.findall(r"^\s*MUM a (\d+) b \d+ len (\d+) .* excess (\d+)", txt, re.M)]
    runs = [(int(x[0]), int(x[1]), int(x[2]), int(x[3])) for x in
            re.findall(r"^\s*MUM a (\d+) b \d+ len (\d+) .* excess (\d+) .*optimal exact run (\d+)", txt, re.M)]
    edge = 0.05 * r["n"]
    r["mum_flank"] = [e == 0 for p, L, e in mums if p + L <= edge or p >= r["n"] - edge]
    r["mum_inner"] = [e == 0 for p, L, e in mums if not (p + L <= edge or p >= r["n"] - edge)]
    inner = [(L, e == 0) for p, L, e in mums if not (p + L <= edge or p >= r["n"] - edge)]
    r["longest_inner_mum"] = max(inner)[1] if inner else None
    edge_ = edge
    r["inner_runs"] = [(L, e, run) for p, L, e, run in runs if not (p + L <= edge_ or p >= r["n"] - edge_)]
    return r


def main():
    args = sys.argv[1:]
    opts = {"--pairs": 3, "--jobs": 12, "--max-cells": 60}
    for k in list(opts):
        if k in args:
            opts[k] = int(args[args.index(k) + 1])
            del args[args.index(k):args.index(k) + 2]
    prices = "1,0,1"
    if "--prices" in args:
        prices = args[args.index("--prices") + 1]
        del args[args.index("--prices"):args.index("--prices") + 2]
    ve, anchors, placement, out = args
    os.makedirs(out, exist_ok=True)
    rids = sorted(r for r in os.listdir(os.path.join(ve, "regions"))
                  if os.path.exists(os.path.join(ve, "regions", r, "hap32.fa")))
    tasks = [(ve, anchors, placement, out, rid, opts["--pairs"], opts["--max-cells"] * 1_000_000, prices)
             for rid in rids]
    rows = []
    with Pool(opts["--jobs"]) as pool:
        for res in pool.imap_unordered(one, tasks):
            for row in res:
                p = parse(*row)
                if p:
                    rows.append(p)
    print("pairs", len(rows), "regions", len({r["rid"] for r in rows}))
    opt = sum(r["opt"] for r in rows)
    for key in ("split", "flanks", "mem", "mumsmem", "mum"):
        rs = [r for r in rows if key in r]
        ex = sum(r[key]["excess"] for r in rs)
        zero = sum(1 for r in rs if r[key]["excess"] == 0)
        anchors_ = sum(r[key]["anchors"] for r in rs)
        right = sum(r[key]["right"] for r in rs)
        fr = sum(r[key]["forced_right"] for r in rs)
        fo = sum(r[key]["forced"] for r in rs)
        print("%-8s excess %d of optima %d (%.3f%%), pairs with none %d/%d, anchors on an optimum %d/%d, forced letters"
              " right %.1f%%" % (key, ex, opt, 100.0 * ex / max(opt, 1), zero, len(rs), right, anchors_,
                                 100.0 * fr / max(fo, 1)))
    # Corridors: the narrowest holding an optimum; the proved width; and the stopping rules.
    for key in ("split", "flanks", "mem", "straight"):
        rs = [r for r in rows if r.get(key)]
        within = {cap: sum(1 for r in rs if r[key]["band"] <= cap) for cap in (27, 54, 108, 216, 432, 864, 1728)}
        print("corridor %-8s narrowest with an optimum: %s of %d" % (key, " ".join("<=%d:%d" % kv for kv in within.items()),
                                                               len(rs)))
        proved = [r for r in rs if "proved" in r[key]]
        if proved:
            cells = sorted(r[key]["proved_cells"] for r in proved)
            print("   proved: median width %d, cells median %d%% (90%% %d%%) of the full DP"
                  % (sorted(r[key]["proved"] for r in proved)[len(proved) // 2], cells[len(cells) // 2],
                     cells[int(0.9 * len(cells))]))
        for k in (2, 3):
            ex = 0
            short = 0
            cells = 0
            full = 0
            for r in rs:
                bands = r[key]["bands"]
                at = stop(bands, k)
                ex += bands[at][1]
                short += bands[at][1] > 0
                # Cells: every width tried up to the stop and the k - 1 doublings that confirm it.
                for x, _ in bands[:min(len(bands), at + k)]:
                    cells += min(r["m"] + 1, 2 * x + 1 + abs(r["n"] - r["m"]) * 0) * (r["n"] + 1)
                full += (r["n"] + 1) * (r["m"] + 1)
            print("   stop when %d widths agree: excess %d of optima %d, pairs short %d, cells <= %.1f%% of the full DP"
                  % (k, ex, opt, short, 100.0 * cells / full))
    fl = [x for r in rows for x in r["mum_flank"]]
    inn = [x for r in rows for x in r["mum_inner"]]
    li = [r["longest_inner_mum"] for r in rows if r["longest_inner_mum"] is not None]
    print("MUMs on an optimum: flanks %d/%d, in the array %d/%d; the array's longest MUM %d/%d"
          % (sum(fl), len(fl), sum(inn), len(inn), sum(li), len(li)))
    wrong = [(L, run) for r in rows for L, e, run in r["inner_runs"] if e > 0]
    if wrong:
        longer = sum(1 for L, run in wrong if L > run)
        print("the array's wrong MUMs: %d; longer than the longest exact run an optimal alignment has in their span:"
              " %d (%.0f%%); median length %d against that run %d" % (len(wrong), longer, 100.0 * longer / len(wrong),
              sorted(L for L, _ in wrong)[len(wrong) // 2], sorted(run for _, run in wrong)[len(wrong) // 2]))
    with open(os.path.join(out, "summary.tsv"), "w") as f:
        f.write("region\tperiod\tn\tm\topt\tsplit_excess\tsplit_band\tflank_band\tstraight_band\n")
        for r in rows:
            f.write("%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n" % (r["rid"], r["period"], r["n"], r["m"], r["opt"],
                    r.get("split", {}).get("excess"), r.get("split", {}).get("band"),
                    r.get("flanks", {}).get("band"), r.get("straight", {}).get("band")))


if __name__ == "__main__":
    main()
