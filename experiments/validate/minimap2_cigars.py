#!/usr/bin/env python3
"""`plegma ref -F` against minimap2 itself: each allele of a region aligned to the region's reference (its first
record) between the reference's flanks, as pgrealign runs it, and the pair's CIGAR over flank + sequence + flank
compared with minimap2's own alignment of the two (-x asm5 at dipcall's settings), which is what dipcall-made truth
sets write.

usage: python3 -I experiments/validate/minimap2_cigars.py <plegma> <region.fa> <flanks.fa> [max alleles]

minimap2 is $MINIMAP2 or the one on the PATH (2.31 is what plegma ports). Prints per allele `same`, `differs` (both
CIGARs) or `skip` (minimap2 gave no single end-to-end primary alignment), then the counts. Differences left are ties
at equal cost (minimap2's ksw2 against karan's ties). It found plegma leaving each gap where its window between two
pins put it, where minimap2 shifts it left last (mm_fix_cigar).
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile


def records(path):
    out = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            out.append([line[1:].split()[0], ""])
        elif out:
            out[-1][1] += line.strip().upper()
    return out


def runs(ops):
    out = []
    for o in ops:
        if out and out[-1][1] == o:
            out[-1][0] += 1
        else:
            out.append([1, o])
    return out


def merged(cigar):
    out = []
    for n, o in cigar:
        if n == 0:
            continue
        if out and out[-1][1] == o:
            out[-1][0] += n
        else:
            out.append([n, o])
    return out


def text(cigar):
    return "".join("%d%s" % (n, o) for n, o in cigar)


def main():
    plegma, region, flanks_path = sys.argv[1], sys.argv[2], sys.argv[3]
    most = int(sys.argv[4]) if len(sys.argv) > 4 else 1 << 30
    minimap2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")
    if not minimap2:
        sys.exit("minimap2 not found: set MINIMAP2")
    alleles = records(region)
    left, right = [s for _, s in records(flanks_path)]
    reference = alleles[0][1]
    work = tempfile.mkdtemp()
    tally = {"same": 0, "differs": 0, "skip": 0}
    for name, allele in alleles[1 : 1 + most]:
        pair = os.path.join(work, "pair.fa")
        with open(pair, "w") as f:
            f.write(">reference\n%s\n>%s\n%s\n" % (reference, name, allele))
        out = subprocess.run([plegma, "ref", "-r", "1", "-F", flanks_path, pair], capture_output=True, text=True)
        rows = [r[1] for r in records_text(out.stdout)]
        ops = []
        for a, b in zip(rows[0], rows[1]):
            if a == "-" and b == "-":
                continue
            ops.append("M" if a != "-" and b != "-" else ("D" if b == "-" else "I"))
        mine = merged([[len(left), "M"]] + runs(ops) + [[len(right), "M"]])
        t = os.path.join(work, "t.fa")
        q = os.path.join(work, "q.fa")
        with open(t, "w") as f:
            f.write(">t\n%s\n" % (left + reference + right))
        with open(q, "w") as f:
            f.write(">q\n%s\n" % (left + allele + right))
        got = subprocess.run(
            [minimap2, "-cx", "asm5", "-z200000,10000", "-r2k", "--secondary=no", "-t1", t, q],
            capture_output=True,
            text=True,
        ).stdout
        primary = [line.split("\t") for line in got.splitlines() if "tp:A:P" in line]
        whole_q = len(left) + len(allele) + len(right)
        whole_t = len(left) + len(reference) + len(right)
        if (
            len(primary) != 1
            or primary[0][2] != "0"
            or primary[0][3] != str(whole_q)
            or primary[0][7] != "0"
            or primary[0][8] != str(whole_t)
        ):
            tally["skip"] += 1
            print("skip", name)
            continue
        tag = [x for x in primary[0] if x.startswith("cg:Z:")][0][5:]
        theirs = merged([[int(n), o] for n, o in re.findall(r"(\d+)([MID])", tag)])
        if theirs == mine:
            tally["same"] += 1
            print("same", name)
        else:
            tally["differs"] += 1
            print("differs", name, "plegma", text(mine)[:300], "minimap2", text(theirs)[:300])
    print("%d same, %d differ, %d skipped" % (tally["same"], tally["differs"], tally["skip"]))


def records_text(fasta):
    out = []
    for line in fasta.splitlines():
        if line.startswith(">"):
            out.append([line[1:].split()[0], ""])
        elif out:
            out[-1][1] += line.strip()
    return out


if __name__ == "__main__":
    main()
