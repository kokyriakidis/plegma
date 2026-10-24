#!/usr/bin/env python3
"""pgrealign's realign stage done by plegma: an msas.jsonl.gz for pgrealign replace, from extract's packages.

usage: python3 -I experiments/chr20/realign_plegma.py <plegma> <karan> <packages.jsonl.gz> <out msas.jsonl.gz>
           [--jobs n] [--timeout s] [--flags "..."] [--max-refined n] [--resume]

Every package whose status is `ok` and whose alleles are not all one (pgrealign's `--select any`) is aligned:
- its distinct alleles, the reference allele (CHM13's, else the heaviest) first, by `plegma <flags> -r 1` (`-L -I
  20,64,1` by default; `-A` anchors the MSA on that first allele); past `--max-refined` distinct alleles (100,000 by
  default, so none), or past the time limit, without -I's refinement (a timed-out region then gets three times the
  limit);
- each fragment, a partial walk at the region's edge or inside it, into the row of an allele holding it exactly where
  its side says (a prefix for `left`, a suffix for `right`, anywhere for `internal`), the allele's letters outside it
  gapped; else aligned (karan, unit costs) to the heaviest allele as a prefix, a suffix or an infix of it and threaded
  into the MSA, its extra letters in columns of their own.
Each record holds `id`, `method` ("plegma"), `status`, `seconds` and, when `ok`, `rows` ({member id: row}) and
`columns`, every row checked to spell its member; replace reads `id`, `status` and `rows`. A region past `--timeout`
seconds (900) is `timeout`, one plegma fails on `aligner_failed`; replace leaves both as they are. With `--resume`,
the `ok` records an interrupted run left in `<out>.part` (read up to where it was cut) are kept and their regions
skipped; the others are tried again.
"""
import gzip
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from multiprocessing import Pool

PLEGMA = KARAN = None
FLAGS = ["-L", "-I", "20,64,1"]
TIMEOUT = 900
MAX_REFINED = 100000


def unrefined(flags):
    """`flags` without -I's refinement: -I and its value, the other flags' values kept."""
    out, skip = [], False
    for x in flags:
        if skip:
            skip = False
            if re.match(r"^\d+(,\d+)*$", x):
                continue
        if x == "-I":
            skip = True
            continue
        out.append(x)
    return out


def plegma_rows(seqs, timeout, copies=None, flags=None):
    """plegma's MSA of the distinct sequences `seqs` (a list), as rows in the same order; with `copies`, each written
    as many times as it says (so plegma weighs a sequence by the haplotypes spelling it); `flags` for FLAGS."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "in.fa")
        with open(path, "w") as f:
            for k, s in enumerate(seqs):
                f.write(">s%d\n%s\n" % (k, s))
                for c in range((copies or {}).get(s, 1) - 1):
                    f.write(">s%d_%d\n%s\n" % (k, c, s))
        flags = flags if flags is not None else FLAGS if len(seqs) <= MAX_REFINED else unrefined(FLAGS)
        out = subprocess.run([PLEGMA] + flags + ["-r", "1", path], capture_output=True, text=True, timeout=timeout)
        if out.returncode != 0:
            raise RuntimeError("plegma exited %d: %s" % (out.returncode, out.stderr.strip()[-200:]))
    rows = {}
    name = None
    for line in out.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            rows[name] = []
        elif name is not None:
            rows[name].append(line.strip().upper())
    return ["".join(rows["s%d" % k]) for k in range(len(seqs))]


def karan_pairing(reference, query, mode):
    """For each query letter, the reference letter it is aligned to (= or X), or -1 (inserted)."""
    with tempfile.TemporaryDirectory() as tmp:
        tsv = os.path.join(tmp, "p.tsv")
        with open(tsv, "w") as f:
            f.write("q\t%s\t%s\n" % (reference, query))
        out = subprocess.run([KARAN, "--pairs", tsv, "--mode", mode], capture_output=True, text=True).stdout
    line = out.strip().splitlines()[-1].split("\t")
    i = int(line[5])
    j = int(line[7])
    pair = [-1] * len(query)
    for n, op in re.findall(r"(\d+)([=XIDM])", line[9]):
        for _ in range(int(n)):
            if op == "D":
                i += 1
            elif op == "I":
                j += 1
            else:
                pair[j] = i
                i += 1
                j += 1
    return pair


def place_fragment(rows, order, alleles, fragment):
    """A row for `fragment` in the MSA `rows` (dict, rewritten in place when columns are added)."""
    seq, side = fragment["seq"].upper(), fragment["side"]
    heaviest = sorted(alleles, key=lambda a: -a["weight"])
    for a in heaviest:
        s = a["seq"].upper()
        at = 0 if side == "left" and s.startswith(seq) else (len(s) - len(seq) if side == "right" and s.endswith(seq)
                                                             else (s.find(seq) if side == "internal" else -1))
        if at < 0 or (side == "left" and at != 0):
            continue
        row = rows[a["id"]]
        out, k = [], 0
        for x in row:
            if x == "-":
                out.append("-")
            else:
                out.append(x if at <= k < at + len(seq) else "-")
                k += 1
        return "".join(out)
    # Aligned to the heaviest allele, then threaded: its letters take the allele's columns, the rest new ones.
    a = heaviest[0]
    s = a["seq"].upper()
    mode = {"left": "prefix", "right": "suffix"}.get(side, "infix")
    pair = karan_pairing(s, seq, mode)
    row = rows[a["id"]]
    column_of = [c for c, x in enumerate(row) if x != "-"]
    # The fragment's letters, each with its column, or None (a new column after the last placed column).
    width = len(row)
    placed = []
    last = -1
    for j, i in enumerate(pair):
        if i >= 0 and column_of[i] > last:
            last = column_of[i]
            placed.append(last)
        else:
            placed.append(None)
    # New columns: after column `last` before each unplaced letter.
    insert_after = {}
    last = -1
    for j, c in enumerate(placed):
        if c is None:
            insert_after.setdefault(last, []).append(j)
        else:
            last = c
    new_rows = {m: [] for m in order}
    frag = []
    # The fragment's letter in each column, once (scanning the fragment for every column costs their product).
    letter_at = {}
    for j, p in enumerate(placed):
        if p is not None and p not in letter_at:
            letter_at[p] = seq[j]
    for c in range(-1, width):
        if c >= 0:
            for m in order:
                new_rows[m].append(rows[m][c])
            frag.append(letter_at.get(c, "-"))
        for j in insert_after.get(c, []):
            for m in order:
                new_rows[m].append("-")
            frag.append(seq[j])
    for m in order:
        rows[m] = "".join(new_rows[m])
    return "".join(frag)


def realign(package):
    start = time.time()
    record = {"id": package["id"], "method": "plegma"}
    try:
        alleles = package["alleles"]
        # The reference allele (CHM13's) first, else the heaviest: plegma -A -Q anchors the MSA on the first sequence.
        first = next((a for a in alleles if a["id"] == package.get("ref")), None)
        if first is None and alleles:
            first = max(alleles, key=lambda a: a["weight"])
        distinct = list(dict.fromkeys(([first["seq"].upper()] if first else []) + [a["seq"].upper() for a in alleles]))
        if len(distinct) == 0:
            record.update(status="no_alleles", message="no allele")
            return record
        if len(distinct) == 1:
            msa = [distinct[0]]
        else:
            # WEIGHTED=1: each distinct sequence once per haplotype spelling it (-A -Q's centre tree needs weights).
            copies = None
            if os.environ.get("WEIGHTED") == "1":
                copies = {}
                for a in alleles:
                    copies[a["seq"].upper()] = copies.get(a["seq"].upper(), 0) + a["weight"]
            try:
                msa = plegma_rows(distinct, TIMEOUT, copies)
            except subprocess.TimeoutExpired:
                # Past the time limit refined, the unrefined MSA (three times the limit): the region keeps a plegma
                # MSA rather than its original alignment.
                if unrefined(FLAGS) == FLAGS:
                    raise
                msa = plegma_rows(distinct, 3 * TIMEOUT, copies, unrefined(FLAGS))
                record["refined"] = False
        row_of = dict(zip(distinct, msa))
        order = [a["id"] for a in alleles]
        rows = {a["id"]: row_of[a["seq"].upper()] for a in alleles}
        for frag in package.get("fragments", []):
            order_now = list(rows)
            rows[frag["id"]] = place_fragment(rows, order_now, alleles, frag)
            order.append(frag["id"])
        # Columns gapped in every row go; every row must spell its member.
        width = len(next(iter(rows.values())))
        keep = [c for c in range(width) if any(r[c] != "-" for r in rows.values())]
        rows = {m: "".join(r[c] for c in keep) for m, r in rows.items()}
        members = {a["id"]: a["seq"].upper() for a in alleles}
        members.update({f["id"]: f["seq"].upper() for f in package.get("fragments", [])})
        for m, r in rows.items():
            if r.replace("-", "") != members[m]:
                raise RuntimeError("row %s does not spell its member" % m)
        record.update(status="ok", rows=rows, columns=len(keep))
    except subprocess.TimeoutExpired:
        record.update(status="timeout", message="over %d s" % TIMEOUT)
    except RuntimeError as error:
        record.update(status="aligner_failed", message=str(error))
    except Exception as error:  # a bug: kept, as pgrealign keeps its own
        record.update(status="error", message="%s: %s" % (type(error).__name__, error))
    record["seconds"] = round(time.time() - start, 3)
    return record


def chosen(path, skip=frozenset()):
    with gzip.open(path, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p.get("status") != "ok" or p["id"] in skip:
                continue
            if len({a["seq"].upper() for a in p["alleles"]}) < 2 and not p.get("fragments"):
                continue
            yield p


def kept_records(path):
    """The whole records an interrupted run wrote to `path` (a gzip stream cut anywhere)."""
    out = []
    if not os.path.exists(path):
        return out
    try:
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.endswith("\n"):
                    out.append(json.loads(line))
    except (EOFError, OSError, json.JSONDecodeError):
        pass
    return out


def init(plegma, karan, flags, timeout, max_refined):
    global PLEGMA, KARAN, FLAGS, TIMEOUT, MAX_REFINED
    PLEGMA, KARAN, FLAGS, TIMEOUT, MAX_REFINED = plegma, karan, flags, timeout, max_refined


def main():
    args = sys.argv[1:]

    def option(name, default):
        if name in args:
            k = args.index(name)
            value = args[k + 1]
            del args[k:k + 2]
            return value
        return default

    jobs = int(option("--jobs", "4"))
    timeout = int(option("--timeout", "900"))
    flags = option("--flags", "-L -I 20,64,1").split()
    max_refined = int(option("--max-refined", "100000"))
    resume = "--resume" in args
    args = [a for a in args if a != "--resume"]
    plegma, karan, packages, out = args
    counts = {}
    # The regions done right are kept; a timeout or a failure is tried again.
    kept = [r for r in kept_records(out + ".part") if r.get("status") == "ok"] if resume else []
    if kept:
        os.replace(out + ".part", out + ".kept")
    done = frozenset(r["id"] for r in kept)
    print("resuming: %d records kept" % len(kept), file=sys.stderr, flush=True)
    with Pool(jobs, initializer=init, initargs=(plegma, karan, flags, timeout, max_refined)) as pool, \
            gzip.open(out + ".part", "wt") as f:
        for record in kept:
            f.write(json.dumps(record) + "\n")
            counts[record["status"]] = counts.get(record["status"], 0) + 1
        for n, record in enumerate(pool.imap(realign, chosen(packages, done), chunksize=4), len(kept) + 1):
            f.write(json.dumps(record) + "\n")
            counts[record["status"]] = counts.get(record["status"], 0) + 1
            if n % 500 == 0:
                print("%d regions: %s" % (n, counts), file=sys.stderr, flush=True)
    os.replace(out + ".part", out)
    print("done: %s" % counts, file=sys.stderr)


if __name__ == "__main__":
    main()
