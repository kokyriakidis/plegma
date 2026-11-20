#!/usr/bin/env python3
"""Golden outputs: every command and library call over a fixed set of random inputs, written to a directory, so a
change meant to keep behaviour (a refactor, a move) can be checked byte for byte against the outputs before it.

usage: python3 -I experiments/validate/golden.py <plegma binary> <libplegma.so> <out dir> [cases] [seed] [model.csv]
       diff -r <before> <after>

The inputs are fuzz.py's cases (reads, repeats, anchored local assemblies, lower case, N, reverse strands) and
diploid read sets for the clustered and phased calls. Each output is one file named for its case and call, so a
diff names what changed. With a homopolymer model (shasta2's csv) the calls are run with it as well.
"""
import os
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))
import fuzz  # noqa: E402

POA = {
    "default": [],
    "N": ["-N"],
    "Nn": ["-N", "-n"],
    "L": ["-L"],
    "s": ["-s"],
    "m1": ["-m", "1"],
    "m2": ["-m", "2"],
    "a1": ["-a", "1"],
    "d2": ["-d", "2"],
    "P0": ["-P", "0"],
    "P3": ["-P", "3,2"],
    "B": ["-b", "-1"],
}
REF = {"default": [], "t3": ["-t", "3"], "F": ["-F"]}


def diploid(rng):
    """Two haplotypes a few edits apart and reads of each, labelled (some unlabelled)."""
    length = rng.choice([40, 150, 600])
    a = "".join(rng.choice("ACGT") for _ in range(length))
    if rng.random() < 0.5:
        at = rng.randint(0, length)
        a = a[:at] + rng.choice("ACGT") * rng.randint(4, 14) + a[at:]
    b = list(a)
    for _ in range(rng.randint(1, 4)):
        i = rng.randrange(len(b))
        b[i] = rng.choice("ACGT")
    b = "".join(b)
    reads, labels = [], []
    for i in range(rng.choice([6, 12, 24])):
        hap = i % 2
        reads.append(fuzz.mutate(rng, a if hap == 0 else b, rng.choice([0.0, 0.01, 0.05])))
        labels.append(hap if rng.random() < 0.8 else -1)
    return reads, labels


def run(cmd, out):
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    with open(out, "w") as f:
        f.write("exit %d\n" % proc.returncode)
        if proc.returncode != 0:
            f.write((proc.stderr.strip().splitlines() or ["?"])[-1] + "\n")
        else:
            f.write(proc.stdout)


def main():
    binary, library, out = sys.argv[1], sys.argv[2], sys.argv[3]
    cases = int(sys.argv[4]) if len(sys.argv) > 4 else 60
    rng = random.Random(int(sys.argv[5]) if len(sys.argv) > 5 else 11)
    model_path = sys.argv[6] if len(sys.argv) > 6 else None
    os.environ["PLEGMA_LIBRARY"] = os.path.abspath(library)
    import plegma

    inputs = os.path.join(out, "inputs")
    os.makedirs(inputs, exist_ok=True)
    jobs = []
    library_cases = []
    for c in range(cases):
        kind, seqs = fuzz.make_case(rng)
        fasta = os.path.join(inputs, "case%d.fa" % c)
        with open(fasta, "w") as f:
            for name, s in seqs:
                f.write(">%s strand=%d\n%s\n" % (name, rng.randint(0, 1), s))
        flanks = os.path.join(inputs, "flanks%d.fa" % c)
        with open(flanks, "w") as f:
            for side in ("left", "right"):
                f.write(">%s\n%s\n" % (side, "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 300)))))
        for name, opts in POA.items():
            for result in ("0", "1", "4"):
                jobs.append(([binary, "poa"] + opts + ["-r", result, fasta], "poa_%d_%s_r%s" % (c, name, result)))
        if model_path:
            jobs.append(([binary, "poa", "-N", "-H", model_path, fasta], "poa_%d_H" % c))
        for name, opts in REF.items():
            opts = [x for o in opts for x in ([o, flanks] if o == "-F" else [o])]
            jobs.append(([binary, "ref"] + opts + ["-r", "2", fasta], "ref_%d_%s" % (c, name)))
        anchors = [w.split("=")[1] for name, _ in seqs for w in name.split()[1:] if w.startswith("anchor=")]
        library_cases.append(([s for _, s in seqs], anchors if len(anchors) == len(seqs) else None))
    # A real local assembly (HG00512 HLA-DRB1, ONT reads anchored at both ends, one or the other).
    real = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tests", "data", "hg00512_drb1")
    real = os.path.join(real, "reads.fa")
    real_options = {"default": [], "N": ["-N"], "d2": ["-d", "2"]}
    if model_path:
        real_options["NH"] = ["-N", "-H", model_path]
    for name, opts in real_options.items():
        jobs.append(([binary, "poa"] + opts + ["-r", "2", real], "drb1_%s" % name))
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as pool:
        list(pool.map(lambda job: run(job[0], os.path.join(out, job[1])), jobs))

    model = plegma.Model.load(model_path) if model_path else None
    with open(os.path.join(out, "library"), "w") as f:
        for c, (reads, anchors) in enumerate(library_cases):
            for pinned in (False, True):
                try:
                    got = plegma.consensus(reads, anchors=anchors, pinned=pinned)
                except plegma.Error as e:
                    got = "error " + str(e)
                f.write("consensus %d pinned=%s %r\n" % (c, pinned, got))
                if model:
                    strands = [i % 2 for i in range(len(reads))]
                    got = plegma.consensus(reads, anchors=anchors, strands=strands, model=model, pinned=pinned)
                    f.write("consensus %d H pinned=%s %r\n" % (c, pinned, got))
        for c in range(cases):
            reads, labels = diploid(rng)
            for m in [None, model] if model else [None]:
                tag = "H" if m else "-"
                two = plegma.clustered_consensus(reads, model=m, pinned=False)
                f.write("clustered %d %s %r %r %r\n" % (c, tag, two.consensus, two.members, two.placements))
                two = plegma.phased_consensus(reads, labels, model=m, pinned=False)
                f.write("phased %d %s %r %r %r\n" % (c, tag, two.consensus, two.members, two.placements))
    print("%d command runs, %d library cases written to %s" % (len(jobs), cases, out))


if __name__ == "__main__":
    main()
