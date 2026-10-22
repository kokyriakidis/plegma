#!/usr/bin/env python3
"""Compares vntr-eval scores of plegma variants with the stored baselines, per stratum.

usage: python3 -I experiments/progressive/summarise.py <vntr-eval dir> [--md file] [--only-common] method ...

For each stratum and method: regions scored, median cost/opt (unit) with how many at or below 1.1, median
gap-affine cost/opt, median k-mer redundancy, nodes per kb and SV pieces per path. With --only-common, every
method is summarised over the regions all of them scored. Paired against the first method, the regions a method
beats it in and loses in.
"""
import json
import os
import statistics
import sys

STRATA = ['hotspot_vntr', 'hotspot_other', 'control_vntr_matched', 'control_vntr_correct', 'control_nontr_sv']


def scores(vntr, method):
    out = {}
    folder = os.path.join(vntr, 'results', method)
    if not os.path.isdir(folder):
        return out
    for f in os.listdir(folder):
        if not f.endswith('.json'):
            continue
        try:
            d = json.load(open(os.path.join(folder, f)))
        except ValueError:
            continue
        a = d.get('alignment') or {}
        if a.get('all_cost_over_opt') is None:
            continue
        out[f[:-5]] = {
            'cost': a['all_cost_over_opt'],
            'affine': a.get('affine_all_cost_over_opt'),
            'kmer': (d.get('redundancy') or {}).get('kmer_frac_extra'),
            'nodes': (d.get('size') or {}).get('nodes_per_kb'),
            'sv': (d.get('inflation') or {}).get('sv_pieces_median'),
        }
    return out


def median(values):
    values = [v for v in values if v is not None]
    return statistics.median(values) if values else float('nan')


def main():
    vntr = sys.argv[1]
    rest = sys.argv[2:]
    md = rest[rest.index('--md') + 1] if '--md' in rest else None
    common = '--only-common' in rest
    methods = [a for i, a in enumerate(rest) if not a.startswith('--') and (i == 0 or rest[i - 1] != '--md')]
    strata = {}
    for r in os.listdir(os.path.join(vntr, 'regions')):
        path = os.path.join(vntr, 'regions', r, 'region.json')
        if os.path.exists(path):
            strata[r] = json.load(open(path))['stratum']
    table = {m: scores(vntr, m) for m in methods}
    lines = []
    for stratum in STRATA:
        regions = sorted(r for r, s in strata.items() if s == stratum)
        if common:
            regions = [r for r in regions if all(r in table[m] for m in methods)]
        lines.append(f'### {stratum} ({len(regions)} regions{" scored by every method" if common else ""})\n')
        lines.append('| method | regions | cost/opt median (<= 1.1) | affine median | k-mer extra | nodes/kb | SV pieces '
                     f'| vs {methods[0]}: better / worse |')
        lines.append('| :-- | --: | --: | --: | --: | --: | --: | --: |')
        for m in methods:
            got = [table[m][r] for r in regions if r in table[m]]
            if not got:
                continue
            good = sum(g['cost'] <= 1.1 for g in got)
            paired = [(table[m][r]['cost'], table[methods[0]][r]['cost']) for r in regions
                      if r in table[m] and r in table[methods[0]]]
            better = sum(a < b - 1e-9 for a, b in paired)
            worse = sum(a > b + 1e-9 for a, b in paired)
            lines.append(f'| {m} | {len(got)} | {median([g["cost"] for g in got]):.3f} ({good}) '
                         f'| {median([g["affine"] for g in got]):.3f} | {median([g["kmer"] for g in got]):.3f} '
                         f'| {median([g["nodes"] for g in got]):.0f} | {median([g["sv"] for g in got]):.1f} '
                         f'| {better} / {worse} |')
        lines.append('')
    text = '\n'.join(lines)
    print(text)
    if md:
        open(md, 'w').write(text + '\n')


if __name__ == '__main__':
    main()
