#!/usr/bin/env python3
"""pgrealign's realign stage done by plegma: an msas.jsonl.gz for pgrealign replace, from extract's packages.

usage: python3 -I experiments/chr20/realign_plegma.py <plegma> <karan> <packages.jsonl.gz> <out msas.jsonl.gz>
           [--jobs n] [--timeout s] [--flags "..."] [--max-refined n] [--resume] [--reference fasta]

Every package whose status is `ok` and whose alleles are not all one (pgrealign's `--select any`) is aligned:
- its distinct alleles, the reference allele (CHM13's, else the heaviest) first, by `plegma <flags> -r 1` (`-L -I
  20,64,1` by default; `-A` anchors the MSA on that first allele); past `--max-refined` distinct alleles (100,000 by
  default, so none), or past the time limit, without -I's refinement (a timed-out region then gets three times the
  limit);
- each fragment, a partial walk at the region's edge or inside it, into the row of an allele holding it exactly where
  its side says (a prefix for `left`, a suffix for `right`, anywhere for `internal`), the allele's letters outside it
  gapped; else aligned (karan, unit costs) to the heaviest allele as a prefix, a suffix or an infix of it and threaded
  into the MSA, its extra letters in columns of their own.
A bare `-F` among the flags (-A's pinned alignment, which seeds with the reference's flanks) takes, per package,
FLANK letters of `--reference` each side of its span (`contig`, `ref_start`, `ref_end`); a package without them is
aligned without -F.
Each record holds `id`, `method` ("plegma"), `status`, `seconds` and, when `ok`, `rows` ({member id: row}) and
`columns`, every row checked to spell its member; replace reads `id`, `status` and `rows`. A region past `--timeout`
seconds (900) is `timeout`, one plegma fails on `aligner_failed`; replace leaves both as they are. With `--resume`,
the `ok` records an interrupted run left in `<out>.part` (read up to where it was cut) are kept and their regions
skipped; the others are tried again.
"""
import collections
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
REFERENCE = {}
FLANK = 1500
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


def plegma_rows(seqs, timeout, copies=None, flags=None, flanks=None):
    """plegma's MSA of the distinct sequences `seqs` (a list), as rows in the same order; with `copies`, each written
    as many times as it says (so plegma weighs a sequence by the haplotypes spelling it); `flags` for FLAGS; `flanks`
    (left, right) the file of a bare -F, which is dropped without them."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "in.fa")
        with open(path, "w") as f:
            for k, s in enumerate(seqs):
                f.write(">s%d\n%s\n" % (k, s))
                for c in range((copies or {}).get(s, 1) - 1):
                    f.write(">s%d_%d\n%s\n" % (k, c, s))
        flags = list(flags if flags is not None else FLAGS if len(seqs) <= MAX_REFINED else unrefined(FLAGS))
        if "-F" in flags:
            at = flags.index("-F")
            if flanks:
                flanks_path = os.path.join(tmp, "flanks.fa")
                with open(flanks_path, "w") as f:
                    f.write(">left\n%s\n>right\n%s\n" % flanks)
                flags.insert(at + 1, flanks_path)
            else:
                del flags[at]
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


def karan_pairings(reference, queries, mode):
    """For each query, for each of its letters, the reference letter it is aligned to (= or X), or -1 (inserted): one
    karan run for all of them."""
    if not queries:
        return []
    with tempfile.TemporaryDirectory() as tmp:
        tsv = os.path.join(tmp, "p.tsv")
        with open(tsv, "w") as f:
            for k, q in enumerate(queries):
                f.write("q%d\t%s\t%s\n" % (k, reference, q))
        out = subprocess.run([KARAN, "--pairs", tsv, "--mode", mode], capture_output=True, text=True, check=True).stdout
    pairs = [None] * len(queries)
    for line in out.strip().splitlines():
        w = line.split("\t")
        if not w[0].startswith("q") or not w[0][1:].isdigit() or len(w) < 10:
            continue
        k = int(w[0][1:])
        i, j = int(w[5]), int(w[7])
        pair = [-1] * len(queries[k])
        for n, op in re.findall(r"(\d+)([=XIDM])", w[9]):
            for _ in range(int(n)):
                if op == "D":
                    i += 1
                elif op == "I":
                    j += 1
                else:
                    pair[j] = i
                    i += 1
                    j += 1
        pairs[k] = pair
    if any(p is None for p in pairs):
        raise RuntimeError("karan returned %d of %d pairings" % (sum(p is not None for p in pairs), len(pairs)))
    return pairs


def place_fragments(rows, order, alleles, fragments):
    """Rows for all of a region's `fragments` in the MSA `rows` (dict, rewritten once when columns are added): a
    fragment an allele holds whole where its side says takes that allele's columns; any other is aligned (karan, one
    run per mode) to the heaviest allele, its letters taking the allele's columns and the rest columns of their own
    after the last column it placed. Every row is built once: placing fragments one by one rebuilt every row, the
    fragments placed before included, per fragment (hundreds of fragments times hundreds of rows times kilobases)."""
    heaviest = sorted(alleles, key=lambda a: -a["weight"])
    width = len(next(iter(rows.values())))
    # Per fragment: the column of each placed letter, and the letters after each column that need columns of their own.
    placements = {}
    aligned = collections.defaultdict(list)
    for frag in fragments:
        seq, side = frag["seq"].upper(), frag["side"]
        for a in heaviest:
            s = a["seq"].upper()
            at = 0 if side == "left" and s.startswith(seq) else (len(s) - len(seq) if side == "right" and
                                                                 s.endswith(seq) else (s.find(seq) if side ==
                                                                                       "internal" else -1))
            if at < 0 or (side == "left" and at != 0):
                continue
            column_of = [c for c, x in enumerate(rows[a["id"]]) if x != "-"]
            placements[frag["id"]] = ({column_of[at + j]: seq[j] for j in range(len(seq))}, {})
            break
        else:
            aligned[{"left": "prefix", "right": "suffix"}.get(side, "infix")].append(frag)
    base = heaviest[0]
    column_of = [c for c, x in enumerate(rows[base["id"]]) if x != "-"]
    for mode, frags in aligned.items():
        seqs = [f["seq"].upper() for f in frags]
        for frag, seq, pair in zip(frags, seqs, karan_pairings(base["seq"].upper(), seqs, mode)):
            letter_at, insert_after = {}, collections.defaultdict(list)
            last = -1
            for j, i in enumerate(pair):
                if i >= 0 and column_of[i] > last:
                    last = column_of[i]
                    letter_at[last] = seq[j]
                else:
                    insert_after[last].append(seq[j])
            placements[frag["id"]] = (letter_at, insert_after)
    # The columns: each old column, then after it every fragment's own letters there, fragment by fragment.
    ids = [f["id"] for f in fragments]
    out = {m: [] for m in list(rows) + ids}
    for c in range(-1, width):
        if c >= 0:
            for m in rows:
                out[m].append(rows[m][c])
            for f in ids:
                out[f].append(placements[f][0].get(c, "-"))
        for f in ids:
            extra = placements[f][1].get(c)
            if not extra:
                continue
            for m in out:
                if m == f:
                    out[m].extend(extra)
                else:
                    out[m].extend("-" * len(extra))
    return {m: "".join(x) for m, x in out.items()}


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
            flanks = None
            contig = REFERENCE.get(package.get("contig"))
            if contig is not None and package.get("ref_end"):
                s0, e0 = int(package["ref_start"]), int(package["ref_end"])
                flanks = (contig[max(0, s0 - FLANK) : s0].upper(), contig[e0 : e0 + FLANK].upper())
            try:
                msa = plegma_rows(distinct, TIMEOUT, copies, flanks=flanks)
            except subprocess.TimeoutExpired:
                # Past the time limit refined, the unrefined MSA (three times the limit): the region keeps a plegma
                # MSA rather than its original alignment.
                if unrefined(FLAGS) == FLAGS:
                    raise
                msa = plegma_rows(distinct, 3 * TIMEOUT, copies, unrefined(FLAGS), flanks)
                record["refined"] = False
        row_of = dict(zip(distinct, msa))
        order = [a["id"] for a in alleles]
        rows = {a["id"]: row_of[a["seq"].upper()] for a in alleles}
        fragments = package.get("fragments", [])
        if fragments:
            rows = place_fragments(rows, order, alleles, fragments)
            order += [f["id"] for f in fragments]
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


def init(plegma, karan, flags, timeout, max_refined, reference):
    global PLEGMA, KARAN, FLAGS, TIMEOUT, MAX_REFINED, REFERENCE
    PLEGMA, KARAN, FLAGS, TIMEOUT, MAX_REFINED = plegma, karan, flags, timeout, max_refined
    REFERENCE = {}
    if reference and "-F" in flags:
        name, parts = None, []
        for line in open(reference):
            if line.startswith(">"):
                if name is not None:
                    REFERENCE[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
        if name is not None:
            REFERENCE[name] = "".join(parts)


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
    reference = option("--reference", "")
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
    with Pool(jobs, initializer=init, initargs=(plegma, karan, flags, timeout, max_refined, reference)) as pool, \
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
