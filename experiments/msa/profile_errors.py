#!/usr/bin/env python3
"""Why each consensus that is still wrong is wrong: a search the polishing left short, or an objective that prefers it.

usage: python3 -I experiments/msa/profile_errors.py <karan> <set dir> <config> <jobs>

For every `<set>.fa` with a `<set>.truth.fa` and a `<set>.fa.<config>.cons.fa` whose consensus differs from the
truth, prints the set, its errors by kind (errors.py), and the reads' summed edit distance to the consensus and to
the truth. When the truth's is lower, the polishing stopped at a local optimum (`search`); when it is not, unit
edit distance itself prefers the consensus (`objective`), and only a different objective would mend it. Ends with
the totals of each verdict by kind.
"""
import collections
import glob
import os
import subprocess
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))


def summed(karan, reference, reads):
    out = subprocess.run([karan, reference, reads, "--distance"], capture_output=True, text=True).stdout
    return sum(int(line.split("\t")[2]) for line in out.strip().splitlines()[1:])


def profile(task):
    karan, fasta, config = task
    truth = fasta[:-3] + ".truth.fa"
    cons = fasta + "." + config + ".cons.fa"
    kinds = subprocess.run(["python3", "-I", os.path.join(HERE, "errors.py"), karan, truth, cons],
                           capture_output=True, text=True).stdout.split()
    if len(kinds) != 5 or kinds[0] == "0":
        return None
    return os.path.basename(fasta), [int(x) for x in kinds], summed(karan, cons, fasta), summed(karan, truth, fasta)


def main():
    karan, folder, config, jobs = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    sets = [f for f in sorted(glob.glob(os.path.join(folder, "*.fa")))
            if not f.endswith(".truth.fa") and ".cons." not in f
            and os.path.exists(f[:-3] + ".truth.fa") and os.path.exists(f + "." + config + ".cons.fa")]
    with Pool(jobs) as pool:
        rows = [r for r in pool.map(profile, [(karan, f, config) for f in sets]) if r]
    totals = collections.defaultdict(lambda: [0] * 5)
    print("| set | errors | substitution | homopolymer | tandem | other | reads to consensus | to truth | verdict |")
    print("| :-- |" + " --: |" * 7 + " :-- |")
    for name, kinds, to_cons, to_truth in sorted(rows, key=lambda r: -r[1][0]):
        verdict = "search" if to_truth < to_cons else "objective"
        for k in range(5):
            totals[verdict][k] += kinds[k]
        print("| %s | %s | %d | %d | %s |" % (name, " | ".join(str(x) for x in kinds), to_cons, to_truth, verdict))
    print("\n| verdict | sets | errors | substitution | homopolymer | tandem | other |")
    print("| :-- | --: |" + " --: |" * 5)
    for verdict in ("search", "objective"):
        sets_of = sum(1 for _, k, c, t in rows if ("search" if t < c else "objective") == verdict)
        print("| %s | %d | %s |" % (verdict, sets_of, " | ".join(str(x) for x in totals[verdict])))


if __name__ == "__main__":
    main()
