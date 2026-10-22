#!/usr/bin/env python3
"""Runs build/pairwise on every vntr-eval region package, CHM13 as the partner.

usage: python3 -I experiments/pairwise/run.py <vntr-eval>/regions <out dir> [jobs] [--chains] [--stratum S] [region ...]

`--chains` skips E1; `--stratum` keeps one stratum's regions.

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
    regions, out, region, chains_only = task
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
         str(meta['flank_right_bp']), str(period), prefix] + (['chains'] if chains_only else []),
        capture_output=True, text=True)
    status = 'ok' if result.returncode == 0 else 'failed: ' + result.stderr.strip().splitlines()[-1][:200]
    if result.returncode == 0:
        open(prefix + '.done', 'w').close()
    return row + [status, round(time.time() - start, 1)]


def main():
    regions, out = sys.argv[1], sys.argv[2]
    jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    rest = sys.argv[4:]
    chains_only = '--chains' in rest
    stratum = rest[rest.index('--stratum') + 1] if '--stratum' in rest else None
    rest = [a for i, a in enumerate(rest) if a not in ('--chains', '--stratum') and (i == 0 or rest[i - 1] != '--stratum')]
    os.makedirs(out, exist_ok=True)
    names = rest or sorted(d for d in os.listdir(regions) if os.path.exists(os.path.join(regions, d, 'region.json')))
    if stratum:
        names = [r for r in names if json.load(open(os.path.join(regions, r, 'region.json')))['stratum'] == stratum]
    # Largest first, so the long regions do not finish last on one core.
    names.sort(key=lambda r: -os.path.getsize(os.path.join(regions, r, 'hap32.fa')))
    with Pool(jobs) as pool, open(os.path.join(out, 'runs.tsv'), 'a') as log:
        writer = csv.writer(log, delimiter='\t')
        for row in pool.imap_unordered(run, [(regions, out, r, chains_only) for r in names]):
            writer.writerow(row)
            log.flush()
            print('\t'.join(map(str, row)), flush=True)


if __name__ == '__main__':
    main()
