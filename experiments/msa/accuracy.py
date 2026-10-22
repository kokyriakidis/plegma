#!/usr/bin/env python3
"""Consensus accuracy against the truth over the benchmark's simulated sets, for several configurations.

usage: python3 -I experiments/msa/accuracy.py <data dir> <karan> <jobs> <max length> name=command ...
           [--only regex] [--cache dir]

Each `name=command` is a shell command with `{in}` for the input FASTA (its consensus on stdout). For every
simulated set up to `max length` letters that has a truth file, every command runs once; the consensus is scored
by its unit edit distance to the truth with karan. Prints, per configuration, the total errors, the sets it
reproduces exactly, its worst set and its total time, and per set the errors of each.

With `--cache dir`, a run whose command, input and programs are unchanged is not run again: its consensus and
time come from the cache, keyed by the command (with `{in}` unfilled), the input's bytes and the bytes of every
file the command names (the program, a wrapper script). A rival's runs are so made once, and plegma's again
whenever its binary changes. Cached times are those of the run that filled the cache.
"""
import glob
import hashlib
import json
import shutil
import os
import re
import shlex
import subprocess
import sys
import time
from multiprocessing import Pool


def length_of(path):
    m = re.search(r"_L(\d+)_", os.path.basename(path))
    return int(m.group(1)) if m else 0


def digest_of(path, memo={}):
    stat = os.stat(path)
    key = (path, stat.st_size, stat.st_mtime_ns)
    if key not in memo:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                h.update(block)
        memo[key] = h.hexdigest()
    return memo[key]


def cache_key(command, fasta):
    h = hashlib.sha256(command.encode())
    h.update(digest_of(fasta).encode())
    for token in shlex.split(command.replace("{in}", "")):
        if os.path.isfile(token):
            h.update(token.encode())
            h.update(digest_of(token).encode())
    return h.hexdigest()


def run(task):
    name, command, fasta, truth, karan, cache = task
    out = fasta + "." + name + ".cons.fa"
    key = cache_key(command, fasta) if cache else None
    stored = os.path.join(cache, key) if cache else None
    if stored and os.path.exists(stored + ".json"):
        meta = json.load(open(stored + ".json"))
        seconds, returncode = meta["seconds"], meta["returncode"]
        if returncode == 0:
            shutil.copyfile(stored + ".cons.fa", out)
    else:
        start = time.time()
        with open(out, "w") as f:
            proc = subprocess.run(command.replace("{in}", shlex.quote(fasta)), shell=True, stdout=f,
                                  stderr=subprocess.DEVNULL, timeout=3600)
        seconds, returncode = time.time() - start, proc.returncode
        if stored:
            if returncode == 0:
                shutil.copyfile(out, stored + ".cons.fa")
            with open(stored + ".json.part", "w") as f:
                json.dump({"seconds": seconds, "returncode": returncode, "command": command, "input": fasta}, f)
            os.replace(stored + ".json.part", stored + ".json")
    if returncode != 0:
        return name, fasta, None, seconds
    lines = subprocess.run([karan, truth, out], capture_output=True, text=True).stdout.strip().splitlines()
    try:
        errors = int(lines[-1].split("\t")[3])
    except (IndexError, ValueError):
        errors = None
    return name, fasta, errors, seconds


def main():
    data, karan, jobs, max_len = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    rest = sys.argv[5:]
    only = rest[rest.index("--only") + 1] if "--only" in rest else None
    cache = os.path.expanduser(rest[rest.index("--cache") + 1]) if "--cache" in rest else None
    if cache:
        os.makedirs(cache, exist_ok=True)
    configs = [a.split("=", 1) for a in rest if "=" in a and not a.startswith("--")]
    sets = []
    for fasta in sorted(glob.glob(os.path.join(data, "sim_*", "*.fa"))):
        if fasta.endswith(".truth.fa") or ".cons." in fasta:
            continue
        truth = fasta[:-3] + ".truth.fa"
        if os.path.exists(truth) and length_of(fasta) <= max_len and (not only or re.search(only, fasta)):
            sets.append((fasta, truth))
    tasks = [(name, cmd, fasta, truth, karan, cache) for name, cmd in configs for fasta, truth in sets]
    with Pool(jobs) as pool:
        results = pool.map(run, tasks)
    table = {}
    for name, fasta, errors, seconds in results:
        table.setdefault(name, {})[fasta] = (errors, seconds)
    print("| config | total errors | exact sets | failed | worst set | time s |")
    print("| :-- | --: | --: | --: | :-- | --: |")
    for name, _ in configs:
        rows = table[name]
        ok = [(e, f) for f, (e, _) in rows.items() if e is not None]
        worst = max(ok, default=(0, ""))
        print("| %s | %d | %d/%d | %d | %s (%d) | %.1f |" % (
            name, sum(e for e, _ in ok), sum(e == 0 for e, _ in ok), len(rows), len(rows) - len(ok),
            os.path.basename(worst[1]), worst[0], sum(s for _, s in rows.values())))
    print()
    print("| set | " + " | ".join(n for n, _ in configs) + " |")
    print("| :-- |" + " --: |" * len(configs))
    for fasta, _ in sets:
        cells = []
        for name, _ in configs:
            e = table[name][fasta][0]
            cells.append("-" if e is None else str(e))
        print("| %s | %s |" % (os.path.basename(fasta)[:-3], " | ".join(cells)))


if __name__ == "__main__":
    main()
