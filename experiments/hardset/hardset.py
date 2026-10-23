#!/usr/bin/env python3
"""vntr-eval's hard set (its 45 chr20 regions whose realignment agrees worst with the pairwise optima) as a fast
loop: plegma's MSAs scored by pgrealign qc's pair F1, the same pairs and optima vntr-eval's strategies are scored
on (tools/hardset_strategies.py's `optima` and `score_region`), without mapping or genotyping.

usage: python3 -I experiments/hardset/hardset.py <vntr-eval> <packages.jsonl.gz> <work dir> <label> [--flags "..."]
           [--regions a,b,... | --pick n] [--jobs n] [--timeout s] [--star] [--gapnorm] [--diagnose] [--sites n]
           [--against <label>]

- Regions: `--regions`, or the `--pick` (6) cheapest of the 45 (distinct alleles x the longest's length squared).
- Alleles: each region's, from pgrealign extract's packages (the same graph's); a region's optima are cached in
  <work dir>/optima/<id>.jsonl.gz (qc's 20 random allele pairs and CHM13's allele against the 20 heaviest, aligned
  by pgrealign's star.align_pair).
- Runs `plegma <flags> -r 1` (`-A` by default) on the region's distinct alleles, as experiments/chr20/realign_plegma.py
  does, into <work dir>/msa/<label>/<id>.fa, and with `--star` pgrealign's medoid star too (label `star`), which
  scored on the same alleles and pairs (the hard set's `pair_f1_realigned` comes from vntr-eval's own extract, whose
  alleles differ); with `--gapnorm` each MSA is also scored with its gaps left-normalised (vntr-eval's gap_norm.py,
  its `star_gapnorm`), as F1 counts a co-optimal gap placement as wrong.
- With `--diagnose` (implies `--star`): why each MSA scores what it does (experiments/hardset/diagnose.py): the F1 a
  pair would have counting any optimal alignment right (what tie-breaking could recover), the letters no optimal
  alignment aligns so by kind, allele and position, and the MSA's sum of pairs against the star's (objective or
  search), and the sites it gets wrong: every wrong letter on the region's CHM13 coordinate, nearby ones merged, each
  site with its interval, letters, pairs, alleles, kinds, shift in periods, the star's wrong letters there and its
  worst pair spelled out (the `--sites` (5) worst per region printed, all in <work dir>/sites.<label>.tsv, the wrong
  letters' coordinates in <work dir>/wrong/<label>/). Needs ANCHORS (build/anchors).
- With `--against <label>` (implies `--diagnose`): every site of that earlier run re-scored under this one (`fixed`
  under a tenth of its wrong letters left, `better`, `same`, `worse`) and this run's new sites: whether a change to
  the MSA fixed the sites it was meant to. With `--regions` only those regions run, in seconds.
- Prints per region the pair F1 of the original graph and of the star (from hardset/chr20.tsv), and this run's pair
  F1, ref F1, cost ratio and seconds; then the medians.

Needs VNTR_ABPOA (abPOA, for the optima and the star) and VNTR_PGREALIGN (pgrealign's checkout) in the environment.
"""
import collections
import csv
import gzip
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # diagnose.py (python -I leaves it off the path)
HS = None


def load(vntr_eval):
    global HS
    if HS is not None:
        return
    sys.path.insert(0, os.path.join(vntr_eval, "tools"))
    import hardset_strategies  # noqa: E402

    HS = hardset_strategies


def regions(vntr_eval):
    with open(os.path.join(vntr_eval, "hardset", "chr20.tsv")) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def packages(path, ids):
    out = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] in ids:
                out[p["id"]] = p
                if len(out) == len(ids):
                    break
    return out


def alleles_of(package):
    """[(id, seq, weight)], distinct sequences only (the first allele spelling each), as the hard set's FASTA."""
    seen = set()
    out = []
    for a in package["alleles"]:
        s = a["seq"].upper()
        if s in seen:
            continue
        seen.add(s)
        out.append((a["id"], s, a["weight"]))
    return out


def cache_optima(task):
    vntr_eval, region, alleles, path = task
    load(vntr_eval)
    if not os.path.exists(path):
        with tempfile.TemporaryDirectory() as wd:
            HS.optima(region, alleles, path, wd, 9)
    return region["id"]


def plegma_msa(plegma, flags, alleles, timeout):
    """plegma's rows ({allele id: row}) for the distinct `alleles`, and its seconds."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "in.fa")
        with open(path, "w") as f:
            for k, (_, s, _) in enumerate(alleles):
                f.write(">s%d\n%s\n" % (k, s))
        start = time.time()
        out = subprocess.run([plegma] + flags + ["-r", "1", path], capture_output=True, text=True, timeout=timeout)
        seconds = time.time() - start
        if out.returncode != 0:
            raise RuntimeError(out.stderr.strip()[-200:])
    rows = {}
    name = None
    for line in out.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            rows[name] = []
        elif name is not None:
            rows[name].append(line.strip().upper())
    return {a[0]: "".join(rows["s%d" % k]) for k, a in enumerate(alleles)}, seconds


def star_msa(region, alleles):
    with tempfile.TemporaryDirectory() as wd:
        start = time.time()
        rows = HS.s_star(alleles, wd, region["ref_allele"])
        return rows, time.time() - start


def one(task):
    vntr_eval, label, plegma, flags, region, alleles, work, timeout, gapnorm = task
    load(vntr_eval)
    rid = region["id"]
    out = os.path.join(work, "msa", label, rid + ".fa")
    try:
        if label == "star" and os.path.exists(out) and os.path.exists(out + ".seconds"):
            # The star is deterministic: its MSA, once made, is reused.
            rows = read_rows(out)
            seconds = float(open(out + ".seconds").read())
        elif label == "star":
            rows, seconds = star_msa(region, alleles)
            if isinstance(rows, list):
                rows = {a[0]: r for a, r in zip(alleles, rows)}
        else:
            rows, seconds = plegma_msa(plegma, flags, alleles, timeout)
    except subprocess.TimeoutExpired:
        return rid, {"status": "timeout"}
    except Exception as error:  # noqa: BLE001
        return rid, {"status": "failed: %s" % str(error)[-120:]}
    for a, s, _ in alleles:
        if rows[a].replace("-", "") != s:
            return rid, {"status": "row %s does not spell its allele" % a}
    with open(out, "w") as f:
        for a, _, _ in alleles:
            f.write(">%s\n%s\n" % (a, rows[a]))
    with open(out + ".seconds", "w") as f:
        f.write("%.3f" % seconds)
    with gzip.open(os.path.join(work, "optima", rid + ".jsonl.gz"), "rt") as f:
        opt = [json.loads(line) for line in f]
    out = {"status": "ok", "seconds": seconds}
    out.update(scores(rows, opt))
    if gapnorm:
        # The same MSA, every row's gap runs shifted left against the column consensus (vntr-eval's gap_norm.py).
        import gap_norm  # noqa: E402

        names = [a for a, _, _ in alleles]
        moved, _ = gap_norm.normalise([rows[a] for a in names], "left")
        out["gapnorm"] = scores(dict(zip(names, moved)), opt)
    return rid, out


def read_rows(path):
    rows, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            rows[name] = ""
        elif name is not None:
            rows[name] += line
    return rows


def diagnosed(task):
    """One (label, region)'s summary, sites and wrong letters' coordinates (diagnose.py), the star's beside it."""
    vntr_eval, anchors, work, label, region, alleles, ref, ref_start = task
    load(vntr_eval)
    import diagnose  # noqa: E402
    from pgrealign import qc  # noqa: E402

    rid = region["id"]
    mine = os.path.join(work, "msa", label, rid + ".fa")
    star = os.path.join(work, "msa", "star", rid + ".fa")
    if not os.path.exists(mine) or not os.path.exists(star):
        return label, rid, None
    with gzip.open(os.path.join(work, "optima", rid + ".jsonl.gz"), "rt") as f:
        opt = [json.loads(line) for line in f]
    star_rows = read_rows(star)
    rows = star_rows if label == "star" else read_rows(mine)
    return label, rid, diagnose.diagnose(HS, qc, anchors, rows, opt, int(region["period"] or 0), alleles, ref,
                                         ref_start, star_rows)


def scores(rows, opt):
    sc = HS.score_region(rows, opt)
    return {
        "pair_f1": HS.f1(sc["pair_tp"], sc["pair_opt"], sc["pair_msa"]),
        "ref_f1": HS.f1(sc["ref_tp"], sc["ref_opt"], sc["ref_msa"]) if sc["ref_opt"] else None,
        "cost_ratio": sc["pair_cost"] / sc["pair_cost_opt"] if sc["pair_cost_opt"] else None,
    }


def main():
    args = sys.argv[1:]

    def option(name, default):
        if name in args:
            k = args.index(name)
            value = args[k + 1]
            del args[k : k + 2]
            return value
        return default

    flags = option("--flags", "-A").split()
    pick = int(option("--pick", "6"))
    chosen = option("--regions", "")
    jobs = int(option("--jobs", "6"))
    timeout = int(option("--timeout", "1800"))
    base = option("--against", "")
    shown_sites = int(option("--sites", "5"))
    diagnosing = "--diagnose" in args or bool(base)
    with_star = "--star" in args or diagnosing
    gapnorm = "--gapnorm" in args
    args = [a for a in args if a not in ("--star", "--gapnorm", "--diagnose")]
    vntr_eval, package_path, work, label = args[:4]
    plegma = os.environ.get("PLEGMA", os.path.join(os.path.dirname(__file__), "..", "..", "build", "plegma"))
    load(vntr_eval)
    every = regions(vntr_eval)
    if chosen:
        wanted = chosen.split(",")
        picked = [r for r in every if r["id"] in wanted]
    else:
        picked = sorted(every, key=lambda r: int(r["alleles"]) * int(r["longest"]) ** 2)[:pick]
    found = packages(package_path, {r["id"] for r in picked})
    picked = [r for r in picked if r["id"] in found]
    alleles = {r["id"]: alleles_of(found[r["id"]]) for r in picked}
    os.makedirs(os.path.join(work, "optima"), exist_ok=True)
    labels = [label] + (["star"] if with_star else [])
    for lab in labels:
        os.makedirs(os.path.join(work, "msa", lab), exist_ok=True)
    with Pool(jobs) as pool:
        tasks = [(vntr_eval, r, alleles[r["id"]], os.path.join(work, "optima", r["id"] + ".jsonl.gz")) for r in picked]
        list(pool.imap_unordered(cache_optima, tasks))
        results = {}
        tasks = [(vntr_eval, lab, plegma, flags, r, alleles[r["id"]], work, timeout, gapnorm) for lab in labels
                 for r in picked]
        for (rid, res), task in zip(pool.imap(one, tasks), tasks):
            results[(task[1], rid)] = res
    shown = [(lab, False) for lab in labels] + ([(lab, True) for lab in labels] if gapnorm else [])
    names = [lab + ("+gapnorm" if g else "") for lab, g in shown]
    print("%-20s %-9s %7s %6s | %6s %6s | %s" % ("region", "class", "alleles", "graph", "star*", "", " ".join(
        "%-26s" % name for name in names)))
    medians = {name: [] for name in names}
    for r in picked:
        cells = []
        for (lab, g), name in zip(shown, names):
            res = results[(lab, r["id"])]
            if res["status"] == "ok":
                sc = res["gapnorm"] if g else res
                medians[name].append(sc["pair_f1"])
                cells.append("F1 %.4f cost %.3f %6.1fs" % (sc["pair_f1"], sc["cost_ratio"] or 0, res["seconds"]))
            else:
                cells.append("%-26s" % res["status"][:26])
        print("%-20s %-9s %7s %6s | %6s %6s | %s" % (r["id"], r["classes"][:9], len(alleles[r["id"]]),
                                                   r["pair_f1_graph"], r["pair_f1_realigned"], "", " ".join(cells)))
    print("median pair F1: graph %.4f, star* %.4f | %s" % (
        statistics.median(float(r["pair_f1_graph"]) for r in picked),
        statistics.median(float(r["pair_f1_realigned"]) for r in picked),
        ", ".join("%s %.4f (%d)" % (lab, statistics.median(v), len(v)) for lab, v in medians.items() if v)))
    print("(star*: the hard set's pair_f1_realigned, pgrealign's medoid star on vntr-eval's own extract)")
    if diagnosing:
        refs = {}
        for r in picked:
            pkg = found[r["id"]]
            ref_seq = next((a["seq"].upper() for a in pkg["alleles"] if a["id"] == pkg.get("ref")), None)
            ref = next((a for a, seq, _ in alleles[r["id"]] if seq == ref_seq), None)
            refs[r["id"]] = (ref, int(pkg.get("ref_start") or 0))
        report(vntr_eval, work, labels, picked, alleles, refs, jobs, base, shown_sites)


SITE_KEYS = ["region", "start", "end", "wrong", "ties", "pairs", "alleles", "kinds", "shift", "star_wrong", "worst",
             "pair_list"]


def report(vntr_eval, work, labels, picked, alleles, refs, jobs, base, shown_sites):
    """Every label's diagnosis on every region: the summary table, the worst sites, <work>/sites.<label>.tsv and
    <work>/wrong/<label>/<region>.json (the wrong letters' CHM13 coordinates); with `base`, each of its sites
    re-scored under the first label (fixed, better, same, worse) and the first label's sites new since."""
    import diagnose  # noqa: E402

    anchors = os.environ.get("ANCHORS", os.path.join(os.path.dirname(__file__), "..", "..", "build", "anchors"))
    tasks = [(vntr_eval, anchors, work, lab, r, alleles[r["id"]]) + refs[r["id"]] for lab in labels for r in picked]
    with Pool(jobs) as pool:
        found = {(lab, rid): d for lab, rid, d in pool.imap_unordered(diagnosed, tasks)}
    keys = ["f1", "co_f1", "wrong", "ties", "near", "copy", "far", "gapped", "top3", "excess", "sp"]
    print()
    print("diagnosis: f1 (qc); co_f1 (a letter pair right on any optimal alignment: co_f1 - f1, what tie-breaking "
          "could recover); wrong: letters aligned as no optimal alignment aligns them, by kind (shares); ties: "
          "aligned as an optimum does but not as the reference; top3: the wrong letters' share in the 3 alleles most "
          "involved; excess: induced cost over the optima; sp: the MSA's sum of pairs over the star's (< 1: the "
          "objective prefers it)")
    print("%-20s %-12s %6s %6s %6s %6s %5s %5s %5s %5s %5s %7s %5s" % (
        "region", "label", "f1", "co_f1", "wrong", "ties", "near", "copy", "far", "gap", "top3", "excess", "sp"))
    for r in picked:
        for lab in labels:
            d = found.get((lab, r["id"]))
            if d is None:
                print("%-20s %-12s %s" % (r["id"], lab, "no MSA"))
                continue
            sm = d[0]
            w = max(sm["wrong"], 1)
            print("%-20s %-12s %6.3f %6.3f %6d %6d %5.2f %5.2f %5.2f %5.2f %5.2f %7d %5.3f" % (
                r["id"], lab, sm["f1"], sm["co_f1"], sm["wrong"], sm["ties"], sm["near"] / w, sm["copy"] / w,
                sm["far"] / w, sm["gapped"] / w, sm["top3"], sm["excess"], sm.get("sp") or 0))
    print()
    for lab in labels:
        ds = [found[(lab, r["id"])][0] for r in picked if found.get((lab, r["id"]))]
        if not ds:
            continue
        med = lambda k: statistics.median(d[k] for d in ds)  # noqa: E731
        wrong = sum(d["wrong"] for d in ds) or 1
        print("%-12s median f1 %.3f, co_f1 %.3f (ties could recover %.3f), sp vs star %.3f; wrong letters %d (ties "
              "%d): near %.0f%%, copy %.0f%%, far %.0f%%, gapped %.0f%%; top3 alleles hold %.0f%% (median)" % (
                  lab, med("f1"), med("co_f1"), med("co_f1") - med("f1"), med("sp"), wrong,
                  sum(d["ties"] for d in ds), 100 * sum(d["near"] for d in ds) / wrong,
                  100 * sum(d["copy"] for d in ds) / wrong, 100 * sum(d["far"] for d in ds) / wrong,
                  100 * sum(d["gapped"] for d in ds) / wrong, 100 * med("top3")))
        os.makedirs(os.path.join(work, "wrong", lab), exist_ok=True)
        with open(os.path.join(work, "sites.%s.tsv" % lab), "w") as f:
            f.write("\t".join(SITE_KEYS) + "\n")
            for r in picked:
                d = found.get((lab, r["id"]))
                if d is None:
                    continue
                with open(os.path.join(work, "wrong", lab, r["id"] + ".json"), "w") as g:
                    json.dump(d[2], g)
                for site in d[1]:
                    f.write("\t".join(str(site.get(k, r["id"])) if k != "region" else r["id"] for k in SITE_KEYS)
                            + "\n")
    lab = labels[0]
    print()
    print("the worst sites of %s (CHM13 interval; wrong letters; ties nearby; pairs; alleles most involved; kinds; "
          "copy shifts in periods; the star's wrong letters there; its worst pair):" % lab)
    for r in picked:
        d = found.get((lab, r["id"]))
        if d is None:
            continue
        for site in d[1][:shown_sites]:
            print("  %s:%d-%d wrong %d ties %d pairs %d | %s | %s | shift %s | star wrong %s | %s" % (
                r["id"], site["start"], site["end"], site["wrong"], site["ties"], site["pairs"], site["alleles"],
                site["kinds"], site["shift"] or "-", site["star_wrong"], site["worst"]))
    if base:
        path = os.path.join(work, "sites.%s.tsv" % base)
        with open(path) as f:
            before = [dict(x, start=int(x["start"]), end=int(x["end"]), wrong=int(x["wrong"]))
                      for x in csv.DictReader(f, delimiter="\t")]
        changes = []
        tally = collections.Counter()
        letters = [0, 0]
        new_sites = []
        for r in picked:
            d = found.get((lab, r["id"]))
            wrong_path = os.path.join(work, "wrong", base, r["id"] + ".json")
            if d is None or not os.path.exists(wrong_path):
                continue
            wrong_before = json.load(open(wrong_path))
            mine = [x for x in before if x["region"] == r["id"]]
            for site, was, now, status in diagnose.against(mine, wrong_before, d[2]):
                tally[status] += 1
                letters[0] += was
                letters[1] += now
                changes.append((now - was, status, r["id"], site, was, now))
            for site in d[1]:
                if not any(x["start"] <= site["end"] and site["start"] <= x["end"] and x["kinds"] == site["kinds"]
                           and x["shift"] == site["shift"] for x in mine):
                    new_sites.append((r["id"], site))
        print()
        print("%s's %d sites under %s (each site's own pairs, its interval): %s; their wrong letters %d -> %d; new "
              "sites %d (%d wrong letters)" % (base, len(before), lab, ", ".join("%s %d" % kv for kv in
                                                                                 tally.most_common()),
                                              letters[0], letters[1], len(new_sites),
                                              sum(x[1]["wrong"] for x in new_sites)))
        changes.sort(key=lambda x: x[0])
        for delta, status, rid, site, was, now in changes[:shown_sites] + changes[-shown_sites:]:
            print("  %-7s %s:%d-%d %s%s wrong %d -> %d" % (status, rid, site["start"], site["end"], site["kinds"],
                                                          " " + site["shift"] if site["shift"] else "", was, now))
        for rid, site in sorted(new_sites, key=lambda x: -x[1]["wrong"])[:shown_sites]:
            print("  new     %s:%d-%d %s%s wrong %d | %s" % (rid, site["start"], site["end"], site["kinds"],
                                                           " " + site["shift"] if site["shift"] else "", site["wrong"],
                                                           site["worst"]))


if __name__ == "__main__":
    main()
