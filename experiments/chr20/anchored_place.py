#!/usr/bin/env python3
"""A prototype placement for one region: each allele placed hard onto one reliable partner already in the MSA, the
alternative to -A's placement by many partners' votes and a band.

usage: python3 -I experiments/chr20/anchored_place.py <packages.jsonl.gz> <region> <out.msas.jsonl.gz> [--reference]

The heaviest allele (most haplotypes) starts the MSA; then each time the allele nearest (shared 15-mers) one already
placed, aligned to that partner by pgrealign's pairwise aligner (star.align_pair under its scores: seeds chained and
filled, the same alignment qc and replace's star use), its letters taking the partner's columns; the letters the
partner lacks go into the columns between, where the partner has gaps (left first), new columns past those. Every
induced pair with its partner is then that pair's own alignment. Writes a msas.jsonl.gz record ({id, status, rows})
that fragmentation.py and local_call.py read. Needs VNTR_PGREALIGN and VNTR_ABPOA.
"""
import gzip
import json
import os
import sys
import tempfile

K = 15


def kmers(s):
    return {s[i : i + K] for i in range(0, max(0, len(s) - K + 1), 3)}


def main():
    packages, rid, out = sys.argv[1:4]
    # --reference: every allele placed onto CHM13's allele (a star on the reference), not its nearest placed partner.
    on_reference = "--reference" in sys.argv[4:]
    sys.path.insert(0, os.environ["VNTR_PGREALIGN"])
    from pgrealign import star  # noqa: E402

    pkg = None
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] == rid:
                pkg = p
                break
    weight = {}
    for a in pkg["alleles"]:
        weight[a["seq"].upper()] = weight.get(a["seq"].upper(), 0) + a["weight"]
    distinct = sorted(weight, key=lambda s: (-weight[s], -len(s)))
    ref_seq = next(a["seq"].upper() for a in pkg["alleles"] if a["id"] == pkg["ref"])
    if on_reference:
        distinct.remove(ref_seq)
        distinct.insert(0, ref_seq)
    marks = {s: kmers(s) for s in distinct}
    rows = {distinct[0]: list(distinct[0])}
    placed = [distinct[0]]
    left = set(distinct[1:])
    # The nearest placed partner of every allele not yet placed, kept up to date (Prim's order).
    best = {s: (len(marks[s] & marks[distinct[0]]) / max(1, len(marks[s] | marks[distinct[0]])), distinct[0])
            for s in left}
    with tempfile.TemporaryDirectory() as wd, star.Abpoa(os.environ["VNTR_ABPOA"], star.SCORES, workdir=wd) as align:
        while left:
            if on_reference:
                q = max(left, key=lambda s: (weight[s], -len(s)))
                partner = ref_seq
            else:
                q = max(left, key=lambda s: (best[s][0], weight[s]))
                partner = best[q][1]
            left.discard(q)
            pa, pq = star.align_pair(partner, q, align, banded=len(partner) * len(q) > 4e8)
            prow = rows[partner]
            cols = [c for c, x in enumerate(prow) if x != "-"]
            width = len(prow)
            # q's letters: on a partner letter, its column; between partner letters, the partner's gap columns there
            # (left first), then new columns.
            new_row = ["-"] * width
            inserts = []  # (after partner letter index, letters)
            i = -1
            pending = []
            for x, y in zip(pa, pq):
                if x != "-":
                    if pending:
                        inserts.append((i, pending))
                        pending = []
                    i += 1
                    if y != "-":
                        new_row[cols[i]] = y
                elif y != "-":
                    pending.append(y)
            if pending:
                inserts.append((i, pending))
            # Room for each insert: the partner's gap columns between its letters i and i + 1.
            added = []  # (column before which new columns go, count)
            for after, letters in inserts:
                lo = cols[after] + 1 if after >= 0 else 0
                hi = cols[after + 1] if after + 1 < len(cols) else width
                room = list(range(lo, hi))
                for c, y in zip(room, letters):
                    new_row[c] = y
                if len(letters) > len(room):
                    added.append((hi, letters[len(room) :]))
            # New columns, rightmost first so the earlier positions stay valid.
            for at, letters in sorted(added, key=lambda x: -x[0]):
                for r in rows.values():
                    r[at:at] = ["-"] * len(letters)
                new_row[at:at] = letters
            rows[q] = new_row
            placed.append(q)
            for s in left:
                j = len(marks[s] & marks[q]) / max(1, len(marks[s] | marks[q]))
                if j > best[s][0]:
                    best[s] = (j, q)
    record = {"id": rid, "status": "ok", "rows": {a["id"]: "".join(rows[a["seq"].upper()]) for a in pkg["alleles"]}}
    for a in pkg["alleles"]:
        assert record["rows"][a["id"]].replace("-", "") == a["seq"].upper()
    with gzip.open(out, "wt") as f:
        f.write(json.dumps(record) + "\n")
    print("placed %d distinct alleles; %d columns" % (len(placed), len(next(iter(rows.values())))))


if __name__ == "__main__":
    main()
