#!/usr/bin/env python3
"""Runs build/pairwise on every vntr-eval region package, CHM13 as the partner.

usage: python3 -I experiments/pairwise/run.py <vntr-eval>/regions <out dir> [jobs] [region ...]

Writes <out dir>/<region>.census.tsv and .chains.tsv per region, and <out dir>/runs.tsv: each region's
stratum, period, flanks, status and wall time. Regions already done are skipped.
"""
import csv
import json
import os
import subprocess
import sys
import time
from multiprocessing import Pool

BINARY = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'build', 'pairwise')


def run(task):
    regions, out, region = task
    path = os.path.join(regions, region)
    with open(os.path.join(path, 'region.json')) as f:
        meta = json.load(f)
    period = max(1, round(meta.get('period') or 1))
    prefix = os.path.join(out, region)
    row = [region, meta['stratum'], period, meta['flank_left_bp'], meta['flank_right_bp']]
    if os.path.exists(prefix + '.chains.tsv') and os.path.exists(prefix + '.done'):
        return row + ['done', 0.0]
    start = time.time()
    result = subprocess.run(
        [BINARY, os.path.join(path, 'hap32.fa'), region, 'CHM13', str(meta['flank_left_bp']),
         str(meta['flank_right_bp']), str(period), prefix],
        capture_output=True, text=True)
    status = 'ok' if result.returncode == 0 else 'failed: ' + result.stderr.strip().splitlines()[-1][:200]
    if result.returncode == 0:
        open(prefix + '.done', 'w').close()
    return row + [status, round(time.time() - start, 1)]


def main():
    regions, out = sys.argv[1], sys.argv[2]
    jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    os.makedirs(out, exist_ok=True)
    names = sys.argv[4:] or sorted(d for d in os.listdir(regions) if os.path.exists(os.path.join(regions, d, 'region.json')))
    # Largest first, so the long regions do not finish last on one core.
    names.sort(key=lambda r: -os.path.getsize(os.path.join(regions, r, 'hap32.fa')))
    with Pool(jobs) as pool, open(os.path.join(out, 'runs.tsv'), 'a') as log:
        writer = csv.writer(log, delimiter='\t')
        for row in pool.imap_unordered(run, [(regions, out, r) for r in names]):
            writer.writerow(row)
            log.flush()
            print('\t'.join(map(str, row)), flush=True)


if __name__ == '__main__':
    main()
