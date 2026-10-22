#!/usr/bin/env python3
"""Summarises experiments/pairwise/run.py output per stratum.

usage: python3 -I experiments/pairwise/summarise.py <out dir> [--md file]

E1: per anchor kind and minimum length, the anchors per pair, the share in the repeat array, and the share of
tested anchors on some optimal alignment, in the array and in the flanks.
E2-E5: per backbone configuration, pooled filled cost over pooled optimum (1.000 = no loss), the median and 90th
percentile of the per-pair ratio, pairs filled optimally, the largest window left, and time against the optimum's.
"""
import collections
import csv
import glob
import os
import statistics
import sys

STRATA = ['hotspot_vntr', 'hotspot_other', 'control_vntr_matched', 'control_vntr_correct', 'control_nonrepeat_sv']


def quantile(values, q):
    if not values:
        return float('nan')
    values = sorted(values)
    return values[min(len(values) - 1, int(q * len(values)))]


def load(out):
    strata = {}
    with open(os.path.join(out, 'runs.tsv')) as f:
        for row in csv.reader(f, delimiter='\t'):
            strata[row[0]] = row[1]
    census, chains = [], []
    for path in glob.glob(os.path.join(out, '*.census.tsv')):
        with open(path) as f:
            census.extend(csv.DictReader(f, delimiter='\t'))
    for path in glob.glob(os.path.join(out, '*.chains.tsv')):
        with open(path) as f:
            chains.extend(csv.DictReader(f, delimiter='\t'))
    return strata, census, chains


def main():
    out = sys.argv[1]
    md = sys.argv[sys.argv.index('--md') + 1] if '--md' in sys.argv else None
    strata, census, chains = load(out)
    lines = []
    present = [s for s in STRATA if any(strata.get(r['region']) == s for r in chains)]
    present += sorted({strata.get(r['region'], '?') for r in chains} - set(present))

    lines.append('## E1: anchors on some optimal alignment (CHM13 partner, unit costs)\n')
    for stratum in present:
        rows = [r for r in census if strata.get(r['region']) == stratum]
        pairs = len({(r['region'], r['query']) for r in rows})
        lines.append(f'### {stratum} ({pairs} pairs)\n')
        lines.append('| kind | L | anchors/pair | in array | optimal | optimal in array | optimal in flanks |')
        lines.append('| :-- | --: | --: | --: | --: | --: | --: |')
        groups = collections.defaultdict(list)
        for r in rows:
            groups[(r['kind'], int(r['min_len']))].append(r)
        order = ['mum', 'mumreference', 'rare2', 'rare8', 'rare50', 'smem', 'all']
        for (kind, min_len) in sorted(groups, key=lambda k: (order.index(k[0]), k[1])):
            g = groups[(kind, min_len)]
            total = lambda key: sum(int(r[key]) for r in g)
            anchors, in_array = total('anchors'), total('in_array')
            tested, tested_in = total('tested'), total('tested_in_array')
            optimal, optimal_in = total('optimal'), total('optimal_in_array')
            share = lambda a, b: f'{a / b:.2f}' if b else '-'
            lines.append(
                f'| {kind} | {min_len} | {anchors / len(g):.1f} | {share(in_array, anchors)} | {share(optimal, tested)} '
                f'| {share(optimal_in, tested_in)} | {share(optimal - optimal_in, tested - tested_in)} |')
        lines.append('')

    lines.append('## E2-E5: filled backbones against the optimum\n')
    for stratum in present:
        rows = [r for r in chains if strata.get(r['region']) == stratum]
        pairs = len({(r['region'], r['query']) for r in rows})
        lines.append(f'### {stratum} ({pairs} pairs)\n')
        lines.append('| config | pooled cost/opt | median | p90 | optimal pairs | anchors/pair | anchored | '
                     'largest window median / p90 / max | windows > 1 kb | backbone+fill / optimum time |')
        lines.append('| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: |')
        groups = collections.defaultdict(list)
        for r in rows:
            groups[r['config']].append(r)
        for config in sorted(groups, key=lambda c: sum(int(r['filled']) for r in groups[c]) / max(1, sum(int(r['optimum']) for r in groups[c]))):
            g = groups[config]
            filled = sum(int(r['filled']) for r in g)
            optimum = sum(int(r['optimum']) for r in g)
            ratios = [int(r['filled']) / int(r['optimum']) for r in g if int(r['optimum']) > 0]
            exact = sum(int(r['filled']) == int(r['optimum']) for r in g)
            anchored = sum(int(r['anchored']) for r in g) / max(1, sum(int(r['query_len']) for r in g))
            largest = [int(r['largest_window']) for r in g]
            time_ratio = sum(float(r['backbone_ms']) + float(r['fill_ms']) for r in g) / max(1e-9, sum(float(r['optimum_ms']) for r in g))
            lines.append(
                f'| {config} | {filled / max(1, optimum):.3f} | {statistics.median(ratios):.3f} | {quantile(ratios, 0.9):.3f} '
                f'| {exact / len(g):.2f} | {sum(int(r["anchors"]) for r in g) / len(g):.1f} | {anchored:.2f} '
                f'| {statistics.median(largest)} / {quantile(largest, 0.9)} / {max(largest)} '
                f'| {sum(int(r["windows_over_1kb"]) for r in g) / len(g):.2f} | {time_ratio:.2f} |')
        lines.append('')
    text = '\n'.join(lines)
    print(text)
    if md:
        with open(md, 'w') as f:
            f.write(text + '\n')


if __name__ == '__main__':
    main()
