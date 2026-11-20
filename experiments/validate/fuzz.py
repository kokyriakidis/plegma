#!/usr/bin/env python3
"""Runs build/plegma's commands (poa, and ref with and without flanks) over random inputs with many option
combinations and checks every output's invariants.

usage: python3 -I experiments/validate/fuzz.py <plegma binary> <work dir> [cases] [seed]

Checks: the run succeeds; an MSA has one row per input record, all of one length, each spelling its record
(case-insensitive, reverse-complemented when its name says so); a GFA's segments, links and paths are
consistent and each path spells its record; a consensus is non-empty when the input is, and uses only letters
the input has. Prints each failure with the command and keeps its input.
"""
import os
import random
import subprocess
import sys

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def mutate(rng, s, rate):
    out = []
    for c in s:
        r = rng.random()
        if r < rate / 3:
            out.append(rng.choice("ACGT"))
        elif r < 2 * rate / 3:
            continue
        elif r < rate:
            out.append(c)
            out.append(rng.choice("ACGT"))
        else:
            out.append(c)
    return "".join(out) or rng.choice("ACGT")


def make_case(rng):
    kind = rng.choice(["reads", "reads", "repeat", "mixed", "tiny", "lower", "n", "dups", "strand", "one", "anchored",
                       "anchored"])
    length = rng.choice([1, 2, 5, 20, 60, 200, 700, 1500])
    base = "".join(rng.choice("ACGT") for _ in range(length))
    if kind == "repeat":
        unit = "".join(rng.choice("ACGT") for _ in range(rng.choice([2, 7, 18, 34])))
        base = "".join(rng.choice("ACGT") for _ in range(30)) + unit * rng.randint(3, 40) + "".join(
            rng.choice("ACGT") for _ in range(30))
    n = rng.choice([1, 2, 3, 5, 10, 30]) if kind != "one" else 1
    seqs = []
    for i in range(n):
        s = mutate(rng, base, rng.choice([0.0, 0.01, 0.05, 0.15, 0.3]))
        if kind == "repeat":
            unit_cut = rng.randint(0, len(s) // 3)
            s = s[: len(s) - unit_cut] if rng.random() < 0.5 else s
        if kind == "mixed" and rng.random() < 0.3:
            s = s[rng.randint(0, len(s) // 2):]
        if kind == "tiny":
            s = s[: rng.randint(1, 4)]
        if kind == "lower" and rng.random() < 0.5:
            s = s.lower()
        if kind == "n":
            s = "".join("N" if rng.random() < 0.05 else c for c in s)
        if kind == "dups" and i > 0 and rng.random() < 0.5:
            s = seqs[0][1]
        if kind == "strand" and rng.random() < 0.4:
            s = rc(s)
        seqs.append(("s%d" % i, s or "A"))
    if kind == "anchored":
        # A local assembly: reads spanning it, fixed at its left or its right end only, the spanning ones first or
        # none at all (each tagged anchor=... in its header, which poa reads; the other command ignores it).
        seqs = []
        for i in range(rng.choice([1, 3, 8, 20])):
            side = rng.choice(["both", "left", "right"] if rng.random() < 0.6 else ["left", "right"])
            cut = rng.randint(0, len(base))
            part = base if side == "both" else (base[:cut] if side == "left" else base[cut:])
            seqs.append(("s%d anchor=%s" % (i, side), mutate(rng, part, rng.choice([0.0, 0.02, 0.1])) or "A"))
    return kind, seqs


def read_fasta(path):
    out = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            out.append([line[1:].split()[0], ""])
        elif out:
            out[-1][1] += line.strip()
    return out


def check_msa(seqs, path):
    rows = read_fasta(path)
    if len(rows) != len(seqs):
        return "%d rows for %d records" % (len(rows), len(seqs))
    width = len(rows[0][1]) if rows else 0
    for (name, seq), (rname, row) in zip(seqs, rows):
        if len(row) != width:
            return "row %s has %d columns, not %d" % (rname, len(row), width)
        spelled = row.replace("-", "")
        want = rc(seq) if rname.endswith("_reverse_complement") else seq
        if not rname.startswith(name.split()[0]):
            return "row %s for record %s" % (rname, name)
        if spelled.upper() != want.upper():
            return "row %s spells %r, not %r" % (rname, spelled[:40], want[:40])
    return None


def check_gfa(seqs, path, with_consensus):
    segs, links, paths = {}, set(), []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "S":
            segs[f[1]] = f[2]
        elif f[0] == "L":
            links.add((f[1], f[3]))
        elif f[0] == "P":
            paths.append((f[1], f[2].split(",") if f[2] else []))
    expected = len(seqs) + (1 if with_consensus else 0)
    if len(paths) != expected:
        return "%d paths for %d records" % (len(paths), expected)
    for (name, seq), (pname, steps) in zip(seqs, paths):
        ids = [s[:-1] for s in steps]
        for a, b in zip(ids, ids[1:]):
            if (a, b) not in links:
                return "path %s steps %s -> %s without a link" % (pname, a, b)
        spelled = "".join(segs[i] for i in ids)
        want = rc(seq) if pname.endswith("_reverse_complement") else seq
        if spelled.upper() != want.upper():
            return "path %s spells %r, not %r" % (pname, spelled[:40], want[:40])
    return None


def check_cli(binary, work):
    """Fixed checks of the command line itself: gzipped input (headers' anchor= and strand= words included) and FASTQ
    read as the plain FASTA file, CRLF line ends as LF; bad options and values fail with one line on stderr and status 1."""
    import gzip

    failures = 0
    plain = os.path.join(work, "cli.fa")
    rng = random.Random(5)
    base = "".join(rng.choice("ACGT") for _ in range(300))
    with open(plain, "w") as f:
        for i in range(8):
            side = ["both", "left", "right"][i % 3]
            part = base if side == "both" else (base[:200] if side == "left" else base[100:])
            f.write(">r%d anchor=%s strand=%d\n%s\n" % (i, side, i % 2, mutate(rng, part, 0.03)))
    with open(plain, "rb") as f:
        text = f.read()
    with gzip.open(plain + ".gz", "wb") as f:
        f.write(text)
    with open(plain + ".crlf", "wb") as f:
        f.write(text.replace(b"\n", b"\r\n"))
    # The same records as FASTQ, header words and all; a quality line starting with '@' is no header.
    with open(plain + ".fq", "w") as f:
        lines = text.decode().split("\n")
        for k in range(0, len(lines) - 1, 2):
            f.write("@%s\n%s\n+\n%s\n" % (lines[k][1:], lines[k + 1], "@" + "I" * (len(lines[k + 1]) - 1)))
    want = subprocess.run([binary, "poa", "-N", "-r", "2", plain], capture_output=True).stdout
    for variant in (plain + ".gz", plain + ".crlf", plain + ".fq"):
        got = subprocess.run([binary, "poa", "-N", "-r", "2", variant], capture_output=True).stdout
        if got != want or not want:
            failures += 1
            print("FAIL cli: %s does not read as %s" % (variant, plain))
    bad = [["poa", "-r", "7"], ["poa", "-d", "0"], ["poa", "-P", "-1"], ["poa", "-b", "x"], ["poa", "-O", "4,x"],
           ["poa", "-q", "2"], ["poa", "-m", "3"], ["poa", "-Z"], ["ref", "-t", "0"], ["ref", "-m", "1"],
           ["poa", "-l"], ["frob"], ["poa", os.path.join(work, "missing.fa")]]
    for args in bad:
        # The plain file as input (with -l a FASTA, not a list of files), but where the input is the fault.
        cmd = [binary] + args + ([] if args == ["frob"] or args[-1].endswith("missing.fa") else [plain])
        proc = subprocess.run(cmd, capture_output=True, text=True)
        lines = proc.stderr.strip().splitlines()
        if proc.returncode != 1 or len(lines) != 1 or not lines[0].startswith("plegma: "):
            failures += 1
            print("FAIL cli: %s exits %d with %r" % (" ".join(cmd), proc.returncode, proc.stderr[:200]))
    return failures


def main():
    binary, work = sys.argv[1], sys.argv[2]
    cases = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    os.makedirs(work, exist_ok=True)
    poa_options = [[], ["-p"], ["-L"], ["-s"], ["-m", "1"], ["-m", "2"], ["-b", "-1"], ["-a", "1"], ["-d", "3"],
                   ["-O", "0", "-E", "2"], ["-O", "4,0", "-E", "2,0"], ["-b", "0", "-f", "0"], ["-P", "0"], ["-P", "1"],
                   ["-d", "3", "-P", "3"], ["-L", "-p"], ["-N"], ["-n"], ["-N", "-n"], ["-N", "-s"], ["-n", "-s"],
                   ["-N", "-L", "-d", "3"]]
    # A bare -F takes a flanks file made per case.
    ref_options = [[], ["-t", "2"], ["-a", "1"], ["-d", "3"], ["-P", "0"], ["-P", "1"], ["-F"], ["-F", "-t", "3"]]
    failures = check_cli(binary, work)
    runs = 0
    for c in range(cases):
        kind, seqs = make_case(rng)
        inp = os.path.join(work, "case%d.fa" % c)
        with open(inp, "w") as f:
            for name, s in seqs:
                f.write(">%s\n%s\n" % (name, s))
        letters = set("".join(s for _, s in seqs).upper())
        command = "ref" if rng.random() < 0.25 else "poa"
        pool = ref_options if command == "ref" else poa_options
        opts = list(rng.choice(pool))
        if rng.random() < 0.3:
            opts += list(rng.choice(pool))
        if "-F" in opts:
            flanks = os.path.join(work, "flanks%d.fa" % c)
            with open(flanks, "w") as f:
                for side in ("left", "right"):
                    f.write(">%s\n%s\n" % (side, "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 400)))))
            opts = [x for o in opts for x in ([o, flanks] if o == "-F" else [o])]
        if rng.random() < 0.15 and len(seqs) > 2 and command == "poa":
            # Start from the MSA of the first records, then add the rest.
            half = len(seqs) // 2
            start_in = os.path.join(work, "start%d.fa" % c)
            with open(start_in, "w") as f:
                for name, s in seqs[:half]:
                    f.write(">%s\n%s\n" % (name, s))
            start_msa = os.path.join(work, "start%d.msa.fa" % c)
            subprocess.run([binary, "poa", "-r", "1", "-o", start_msa, start_in], capture_output=True)
            with open(inp, "w") as f:
                for name, s in seqs[half:]:
                    f.write(">%s\n%s\n" % (name, s))
            opts += ["-i", start_msa]
        for result in ["0", "1", "3", "4"]:
            out = os.path.join(work, "out%d_%s" % (c, result))
            cmd = [binary, command] + opts + ["-r", result, "-o", out, inp]
            runs += 1
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            except subprocess.TimeoutExpired:
                print("TIMEOUT", kind, " ".join(cmd))
                failures += 1
                continue
            problem = None
            if proc.returncode != 0:
                problem = "exit %d: %s" % (proc.returncode, (proc.stderr.strip().splitlines() or ["?"])[-1][:200])
            elif result == "1":
                problem = check_msa(seqs, out)
            elif result in ("3", "4"):
                problem = check_gfa(seqs, out, result == "4")
            elif result == "0":
                cons = read_fasta(out)
                if not cons or any(not s for _, s in cons):
                    problem = "empty consensus"
                elif any(set(s.upper()) - letters for _, s in cons):
                    problem = "consensus letters outside the input"
            if problem:
                failures += 1
                print("FAIL", kind, problem, "|", " ".join(cmd))
    print("%d runs, %d failures" % (runs, failures))


if __name__ == "__main__":
    main()
