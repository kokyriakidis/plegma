#!/usr/bin/env python3
"""Summary statistics of the vntr-eval region packages.

usage: python3 -I region_stats.py <vntr-eval/regions> <out.tsv>
Reads only FASTA / JSON / TSV data (no code from the repo). Writes one row per region and prints
per-stratum summaries in Markdown.
"""
import csv, gzip, json, os, statistics, sys


def read_fasta(path):
    op = gzip.open if path.endswith('.gz') else open
    name, seq, out = None, [], []
    with op(path, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ''.join(seq)))
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        out.append((name, ''.join(seq)))
    return out


def q(xs, p):
    xs = sorted(xs)
    if not xs:
        return float('nan')
    k = (len(xs) - 1) * p
    lo = int(k)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def main():
    rdir, out = sys.argv[1], sys.argv[2]
    rows = []
    with open(os.path.join(rdir, 'regions.tsv')) as f:
        idx = list(csv.DictReader(f, delimiter='\t'))
    for r in idx:
        d = os.path.join(rdir, r['region_id'])
        js = json.load(open(os.path.join(d, 'region.json')))
        h = read_fasta(os.path.join(d, 'hap32.fa'))
        hl = [len(s) for _, s in h]
        chm = [s for n, s in h if n.startswith('CHM13')][0]
        hp = os.path.join(d, 'hprc.fa.gz')
        fl = read_fasta(hp) if os.path.exists(hp) else []
        fll = [len(s) for _, s in fl]
        frag = os.path.join(d, 'hap32.fragments.fa')
        nfrag = len(read_fasta(frag)) if os.path.exists(frag) and os.path.getsize(frag) else 0
        period = js.get('period')
        row = dict(
            region_id=r['region_id'], stratum=r['stratum'], contig=r['contig'],
            period=period if period not in (None, '.') else '',
            copies=js.get('copies') if js.get('copies') not in (None, '.') else '',
            span_bp=js['span_bp'], core_bp=js['core_bp'],
            flank_left=js['flank_left_bp'], flank_right=js['flank_right_bp'],
            chm13_len=len(chm),
            n_hap32=len(h), distinct_hap32=len(set(s for _, s in h)),
            hap32_min=min(hl), hap32_med=statistics.median(hl), hap32_max=max(hl),
            hap32_spread=max(hl) - min(hl),
            n_frag=nfrag,
            n_full=len(fl), distinct_full=len(set(s for _, s in fl)),
            full_min=min(fll) if fll else '', full_med=statistics.median(fll) if fll else '',
            full_max=max(fll) if fll else '',
            full_bp_distinct=sum(len(s) for s in set(s for _, s in fl)),
            hap32_N=sum(1 for _, s in h if 'N' in s.upper()),
            vg_fp=r['vg_fp'], vg_fn=r['vg_fn'],
        )
        rows.append(row)
    with open(out, 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    strata = ['hotspot_vntr', 'hotspot_other', 'control_vntr_matched', 'control_vntr_correct',
              'control_nontr_sv']
    def fmt(xs):
        xs = [float(x) for x in xs if x != '']
        if not xs:
            return '-'
        return '%g / %g / %g / %g' % (round(min(xs)), round(q(xs, .5)), round(q(xs, .9)), round(max(xs)))
    cols = [('span_bp', 'span bp'), ('chm13_len', 'CHM13 len'), ('core_bp', 'core bp'),
            ('period', 'period'), ('copies', 'copies'),
            ('distinct_hap32', 'distinct hap32 seqs'), ('hap32_med', 'hap32 median len'),
            ('hap32_max', 'hap32 max len'), ('hap32_spread', 'hap32 max-min len'),
            ('n_full', 'full-panel records'), ('distinct_full', 'distinct full-panel seqs'),
            ('full_med', 'full median len'), ('full_max', 'full max len'),
            ('full_bp_distinct', 'distinct full-panel bp')]
    print('| stratum | n | ' + ' | '.join(c[1] for c in cols) + ' |')
    print('|' + '---|' * (len(cols) + 2))
    for s in strata + ['ALL']:
        rr = [x for x in rows if s == 'ALL' or x['stratum'] == s]
        print('| %s | %d | ' % (s, len(rr)) + ' | '.join(fmt([x[c] for x in rr]) for c, _ in cols) + ' |')
    print()
    print('cells: min / median / p90 / max over regions')
    print('total hap32 records', sum(x['n_hap32'] for x in rows),
          'total full-panel records', sum(x['n_full'] for x in rows),
          'regions with fragments', sum(1 for x in rows if x['n_frag']),
          'regions n_hap32<34', sum(1 for x in rows if x['n_hap32'] < 34),
          'regions with N in hap32', sum(1 for x in rows if x['hap32_N']))
    # period bins for VNTR strata
    for s in strata:
        ps = [float(x['period']) for x in rows if x['stratum'] == s and x['period'] != '']
        bins = {'<7': 0, '7-15': 0, '16-30': 0, '31-60': 0, '61-200': 0, '>200': 0}
        for p in ps:
            k = '<7' if p < 7 else '7-15' if p <= 15 else '16-30' if p <= 30 else '31-60' if p <= 60 else '61-200' if p <= 200 else '>200'
            bins[k] += 1
        print(s, 'period bins', bins, 'no period', sum(1 for x in rows if x['stratum'] == s and x['period'] == ''))
    big = sorted(rows, key=lambda x: -float(x['full_max'] or 0))[:12]
    print('largest full-panel alleles:', [(x['region_id'], x['stratum'], x['hap32_max'], x['full_max']) for x in big])


if __name__ == '__main__':
    main()
