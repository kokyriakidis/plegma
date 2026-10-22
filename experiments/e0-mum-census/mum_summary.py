#!/usr/bin/env python3
"""Summarise mum_driver.py output. usage: python3 -I mum_summary.py mum.tsv region_stats.tsv"""
import csv, statistics, sys, collections


def med(xs):
    return statistics.median(xs) if xs else float('nan')


def q(xs, p):
    xs = sorted(xs)
    if not xs:
        return float('nan')
    return xs[min(len(xs) - 1, int(p * (len(xs) - 1) + 0.5))]


rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
rs = {r['region_id']: r for r in csv.DictReader(open(sys.argv[2]), delimiter='\t')}
for r in rows:
    for k in list(r):
        if k not in ('region', 'stratum', 'kind', 'id'):
            r[k] = float(r[k])
strata = ['hotspot_vntr', 'hotspot_other', 'control_vntr_matched', 'control_vntr_correct', 'control_nontr_sv']
print('| stratum | kind | pairs (DP) | pairs w/ core MUM | core MUMs / pair (median) | core MUM bp / array bp (median) |'
      ' MUMs on an optimal path: core / flank | mean excess per core MUM | top chain cost/opt pooled (median) |'
      ' recursive cost/opt pooled (median) | rec max window median / p90 / max | anchored bp / len (median) |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|')
for s in strata:
    for kind in ['ref', 'med', 'rnd']:
        rr = [r for r in rows if r['stratum'] == s and r['kind'] == kind]
        if not rr:
            continue
        dp = [r for r in rr if r['opt'] >= 0 and r['opt'] > 0]
        arr = []
        for r in rr:
            st = rs[r['region']]
            core_len = max(1.0, r['lenA'] - float(st['flank_left']) - float(st['flank_right']))
            arr.append(r['core_bp'] / core_len)
        withcore = sum(1 for r in rr if r['n_core'] > 0)
        zc = sum(r['zero_core'] for r in dp); nc = sum(r['n_core'] for r in dp)
        zf = sum(r['zero_flank'] for r in dp); nf = sum(r['n_flank'] for r in dp)
        exc = sum(r['excess_core'] for r in dp)
        chain_p = sum(r['chain_cost'] for r in dp) / max(1, sum(r['opt'] for r in dp))
        chain_m = med([r['chain_cost'] / r['opt'] for r in dp])
        rec_p = sum(r['rec_cost'] for r in dp) / max(1, sum(r['opt'] for r in dp))
        rec_m = med([r['rec_cost'] / r['opt'] for r in dp])
        win = [r['rec_maxwin'] for r in dp]
        anc = med([r['rec_anchor_bp'] / max(r['lenA'], r['lenB']) for r in dp])
        print('| %s | %s | %d (%d) | %.0f%% | %g | %.2f | %s / %s | %s | %.3f (%.3f) | %.3f (%.3f) | %g / %g / %g | %.2f |' % (
            s, kind, len(rr), len(dp), 100.0 * withcore / len(rr), med([r['n_core'] for r in rr]), med(arr),
            '%.0f%%' % (100.0 * zc / nc) if nc else '-', '%.0f%%' % (100.0 * zf / nf) if nf else '-',
            '%.1f' % (exc / nc) if nc else '-', chain_p, chain_m, rec_p, rec_m,
            med(win), q(win, .9), max(win) if win else float('nan'), anc))
print()
# per-region view for hotspot VNTR ref pairs
print('per hotspot_vntr region (ref pairs): region period copies hap32max | pairs | med core MUMs | %core MUMs on opt | rec cost/opt pooled | rec max window max')
by = collections.defaultdict(list)
for r in rows:
    if r['stratum'] == 'hotspot_vntr' and r['kind'] == 'ref':
        by[r['region']].append(r)
for rid, rr in sorted(by.items(), key=lambda kv: -max(x['rec_maxwin'] for x in kv[1])):
    st = rs[rid]
    dp = [r for r in rr if r['opt'] > 0]
    nc = sum(r['n_core'] for r in dp); zc = sum(r['zero_core'] for r in dp)
    rp = sum(r['rec_cost'] for r in dp) / max(1, sum(r['opt'] for r in dp)) if dp else float('nan')
    print(rid, st['period'], st['copies'], st['hap32_max'], '|', len(rr), '|', med([r['n_core'] for r in rr]), '|',
          '%.0f%%' % (100 * zc / nc) if nc else '-', '| %.3f' % rp, '|', max(r['rec_maxwin'] for r in rr))
