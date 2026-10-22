import csv, statistics, sys
rs = {r['region_id']: r for r in csv.DictReader(open(sys.argv[2]), delimiter='\t')}
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
for s in ['hotspot_vntr', 'hotspot_other', 'control_vntr_matched']:
    for kind in ['ref', 'med']:
        rr = [r for r in rows if r['stratum'] == s and r['kind'] == kind and float(r['opt']) > 0]
        cov, none, ex = [], 0, []
        for r in rr:
            st = rs[r['region']]
            core = max(1.0, float(r['lenA']) - float(st['flank_left']) - float(st['flank_right']))
            cov.append(float(r['chain_core_bp']) / core)
            none += float(r['chain_core_n']) == 0
            ex.append(float(r['rec_cost']) - float(r['opt']))
        # by period bin (hotspot only)
        print(s, kind, 'pairs', len(rr), 'chain-anchored fraction of A array: median %.2f' % statistics.median(cov),
              'pairs with no chained array anchor %.0f%%' % (100 * none / len(rr)),
              'pairs with rec excess 0: %.0f%%' % (100 * sum(1 for e in ex if e == 0) / len(rr)))
rr = [r for r in rows if r['stratum'] == 'hotspot_vntr' and r['kind'] == 'ref' and float(r['opt']) > 0]
bins = {}
for r in rr:
    p = float(rs[r['region']]['period'])
    b = '7-15' if p <= 15 else '16-30' if p <= 30 else '31-60' if p <= 60 else '>60'
    bins.setdefault(b, []).append(r)
for b, xs in sorted(bins.items()):
    opt = sum(float(x['opt']) for x in xs)
    print('hotspot ref period', b, 'pairs', len(xs), 'rec cost/opt pooled %.3f' % (sum(float(x['rec_cost']) for x in xs) / opt),
          'top chain %.3f' % (sum(float(x['chain_cost']) for x in xs) / opt),
          'median rec max window', statistics.median(float(x['rec_maxwin']) for x in xs))
