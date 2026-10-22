#!/usr/bin/env python3
"""Runs build/progressive on vntr-eval regions and scores each MSA with vntr-eval's own tools.

usage: python3 -I experiments/progressive/run.py <vntr-eval dir> <work dir> <jobs> <variant> [<variant> ...]
           [--stratum S] [--regions R1,R2,...] [--timeout seconds]

A variant is `<config>:<order>:<costs>:<partner>:<fill>[:<rounds>]` (see progressive.mojo), named in vntr-eval as
`plegma_<config>_<order>_<costs>_<partner>_<fill>`: the MSA goes to <work dir>/msa, the graph to
<vntr-eval>/candidates/<name>/ and the scores to <vntr-eval>/results/<name>/. <work dir>/progressive.tsv logs
every run: region, variant, status, aligner seconds, nodes, columns.
"""
import csv
import json
import os
import subprocess
import sys
from multiprocessing import Pool

BINARY = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'build', 'progressive')


def name_of(variant):
    return 'plegma_' + variant.replace(':', '_')


def run(task):
    vntr, work, region, variant, timeout = task
    name = name_of(variant)
    msa = os.path.join(work, 'msa', f'{region}.{name}.msa.fa')
    gfa = os.path.join(vntr, 'candidates', name, f'{region}.gfa')
    scores = os.path.join(vntr, 'results', name, f'{region}.json')
    hap32 = os.path.join(vntr, 'regions', region, 'hap32.fa')
    if os.path.exists(scores):
        return [region, variant, 'done', '', '', '']
    try:
        out = subprocess.run([BINARY, hap32, msa] + variant.split(':'), capture_output=True, text=True,
                             timeout=timeout)
    except subprocess.TimeoutExpired:
        return [region, variant, 'timeout', timeout, '', '']
    if out.returncode != 0:
        return [region, variant, 'failed: ' + (out.stderr.strip().splitlines() or ['?'])[-1][:200], '', '', '']
    fields = out.stdout.strip().split('\t')
    seconds, nodes, columns = fields[-1], fields[-4], fields[-3]
    os.makedirs(os.path.dirname(gfa), exist_ok=True)
    graph = subprocess.run(['python3', 'tools/msa_graph.py', msa, hap32, gfa, '--engine', 'native'], cwd=vntr,
                           capture_output=True, text=True)
    if graph.returncode != 0:
        return [region, variant, 'msa_graph failed: ' + graph.stderr.strip()[-200:], seconds, nodes, columns]
    score = subprocess.run(['python3', 'tools/evaluate.py', os.path.join('regions', region), gfa, '--no-truvari',
                            '--skip', 'truth'], cwd=vntr, capture_output=True, text=True)
    status = 'ok' if score.returncode == 0 else 'evaluate failed: ' + score.stderr.strip()[-200:]
    return [region, variant, status, seconds, nodes, columns]


def main():
    vntr, work, jobs = sys.argv[1], sys.argv[2], int(sys.argv[3])
    rest = sys.argv[4:]
    option = lambda flag: rest[rest.index(flag) + 1] if flag in rest else None
    stratum, regions, timeout = option('--stratum'), option('--regions'), int(option('--timeout') or 3600)
    variants = [a for i, a in enumerate(rest) if not a.startswith('--') and (i == 0 or not rest[i - 1].startswith('--'))]
    names = regions.split(',') if regions else sorted(
        d for d in os.listdir(os.path.join(vntr, 'regions')) if os.path.exists(os.path.join(vntr, 'regions', d, 'region.json')))
    if stratum:
        names = [r for r in names
                 if json.load(open(os.path.join(vntr, 'regions', r, 'region.json')))['stratum'] == stratum]
    names.sort(key=lambda r: -os.path.getsize(os.path.join(vntr, 'regions', r, 'hap32.fa')))
    os.makedirs(os.path.join(work, 'msa'), exist_ok=True)
    tasks = [(vntr, work, r, v, timeout) for v in variants for r in names]
    with Pool(jobs) as pool, open(os.path.join(work, 'progressive.tsv'), 'a') as log:
        writer = csv.writer(log, delimiter='\t')
        for row in pool.imap_unordered(run, tasks):
            writer.writerow(row)
            log.flush()
            print('\t'.join(map(str, row)), flush=True)


if __name__ == '__main__':
    main()
