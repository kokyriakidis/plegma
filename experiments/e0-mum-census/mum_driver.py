#!/usr/bin/env python3
"""Pairwise MUM anchoring experiment over the vntr-eval hap32 packages.

usage: python3 -I mum_driver.py <regions dir> <mumcheck binary> <out.tsv> [jobs] [Lmin] [kinds]
Pairs per region (distinct hap32 sequences only):
  ref: CHM13 vs every other distinct sequence
  med: the medoid (least haplotype-weighted multiset 15-mer Jaccard distance, as pgrealign) vs every other
  rnd: up to 15 random pairs of distinct non-CHM13 sequences (seed = region id)
Core zone of CHM13 = the annotated core; of any other sequence = [left flank, len - right flank) (inferred).
"""
import collections, csv, json, os, random, subprocess, sys
from multiprocessing import Pool

MAX_CELLS = 1.5e9


def read_fasta(path):
    name, seq, out = None, [], []
    with open(path) as f:
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


def kmers(s, k=15):
    c = collections.Counter()
    for i in range(len(s) - k + 1):
        w = s[i:i + k]
        if 'N' not in w:
            c[w] += 1
    return c


def medoid(seqs, weights):
    ks = [kmers(s) for s in seqs]
    tot = [sum(c.values()) for c in ks]
    n = len(seqs)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = ks[i], ks[j]
            if len(a) > len(b):
                a, b = b, a
            sh = sum(min(v, b.get(k, 0)) for k, v in a.items())
            den = tot[i] + tot[j] - sh
            d = 1 - (sh / den if den else 1)
            D[i][j] = D[j][i] = d
    best = min(range(n), key=lambda i: (sum(weights[j] * D[i][j] for j in range(n)), i))
    return best


def region_jobs(rdir, rid, Lmin, kinds):
    js = json.load(open(os.path.join(rdir, rid, 'region.json')))
    recs = read_fasta(os.path.join(rdir, rid, 'hap32.fa'))
    cnt = collections.Counter(s.upper() for _, s in recs)
    chm = [s.upper() for n, s in recs if n.startswith('CHM13')][0]
    seqs = list(cnt)
    w = [cnt[s] for s in seqs]
    fl = js['core_start'] - js['span_start']
    fr = js['span_end'] - js['core_end']
    def zone(s):
        if s == chm:
            return fl, js['core_end'] - js['span_start'] + 1
        return fl, max(fl, len(s) - fr)
    lines = []
    def add(kind, i_a, a, b):
        za, zb = zone(a), zone(b)
        dp = 1 if len(a) * len(b) <= MAX_CELLS else 0
        lines.append('%s|%s|%d|%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s' % (
            rid, kind, i_a, cnt[b], Lmin, za[0], za[1], zb[0], zb[1], dp, a, b))
    others = [s for s in seqs if s != chm]
    if 'ref' in kinds:
        for s in others:
            add('ref', 0, chm, s)
    if 'med' in kinds and len(seqs) > 1:
        m = seqs[medoid(seqs, w)]
        for s in seqs:
            if s != m:
                add('med', int(m == chm), m, s)
    if 'rnd' in kinds and len(others) > 2:
        rng = random.Random(rid)
        pairs = [(i, j) for i in range(len(others)) for j in range(i + 1, len(others))]
        rng.shuffle(pairs)
        for i, j in pairs[:15]:
            add('rnd', 0, others[i], others[j])
    return lines


def run_chunk(args):
    binary, lines = args
    p = subprocess.run([binary], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    out = p.stdout.splitlines()
    return out[1:] if out else []


def main():
    rdir, binary, out = sys.argv[1:4]
    jobs = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    Lmin = int(sys.argv[5]) if len(sys.argv) > 5 else 20
    kinds = sys.argv[6].split(',') if len(sys.argv) > 6 else ['ref', 'med', 'rnd']
    with open(os.path.join(rdir, 'regions.tsv')) as f:
        idx = list(csv.DictReader(f, delimiter='\t'))
    strat = {r['region_id']: r['stratum'] for r in idx}
    chunks = []
    for r in idx:
        lines = region_jobs(rdir, r['region_id'], Lmin, kinds)
        for i in range(0, len(lines), 4):
            chunks.append((binary, lines[i:i + 4]))
    # biggest chunks first
    chunks.sort(key=lambda c: -sum(len(l) ** 2 for l in c[1]))
    header = None
    with open(out, 'w') as fo:
        fo.write('region\tstratum\tkind\tcentre_is_chm13\tweight_b\tid\tlenA\tlenB\topt\tlcs\tn_mum\tmum_bp\tn_core\tn_flank'
                 '\tn_cross\tcore_bp\tflank_bp\tzero_core\tzero_flank\tzero_cross\texcess_core\texcess_flank'
                 '\tchain_n\tchain_bp\tchain_core_n\tchain_core_bp\tchain_cost\trec_cost\trec_anchor_bp'
                 '\trec_maxwin\trec_cells\trec_depth\trec_leaves\n')
        with Pool(jobs) as pool:
            for res in pool.imap_unordered(run_chunk, chunks):
                for l in res:
                    f = l.split('\t')
                    rid, kind, ca, wb = f[0].split('|')
                    fo.write('\t'.join([rid, strat[rid], kind, ca, wb] + f) + '\n')
                fo.flush()


if __name__ == '__main__':
    main()
