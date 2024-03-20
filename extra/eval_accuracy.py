#!/usr/bin/env python3

import numpy as np
import argparse
import sys
import os
import csv
import collections
from tqdm import tqdm

import common


def load_distances(cached_dists, locus, path_fmt):
    if locus in cached_dists:
        return cached_dists[locus]

    path = path_fmt.format(locus)
    if not os.path.exists(path):
        dists = None
    with common.open(path) as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                break
        reader = csv.DictReader(f, delimiter='\t', fieldnames=fields)
        dists = collections.defaultdict(lambda: [{}, None])
        for row in reader:
            curr_dist = dists[row['target']]
            query = row['query']
            if query == '*':
                curr_dist[0] = None
                continue

            query1, query2 = query.split(',')
            div = float(row['divergence'])
            curr_dist[0][(query1, query2)] = div
            curr_dist[0][(query2, query1)] = div
            if curr_dist[1] is None and row['loo'] == 'T':
                curr_dist[1] = div

    cached_dists[locus] = dists
    return dists


def get_qv(val):
    return np.inf if val == 0 else -10 * np.log10(val)


def process_line(dists, sample, locus, genotype, loo):
    sample_dists, avail_div = dists[sample]
    if sample_dists is None:
        return 'NA\tNA\tNA\tNA' if loo else 'NA\tNA'

    assert sample_dists

    query1, query2 = genotype.split(',')
    try:
        div = sample_dists[(query1, query2)]
    except KeyError:
        sys.stderr.write(f'Cannot find genotype {genotype} for {sample} at {locus}\n')
        exit(1)

    s = f'{div:.9f}\t{get_qv(div):.9f}'
    if loo:
        if avail_div is None:
            s += '\tNA\tNA'
        else:
            s += f'\t{avail_div:.9f}\t{get_qv(avail_div):.9f}'
    return s


def main():
    parser = argparse.ArgumentParser(
        description='Evaluating genotyping accuracy',
        usage='%(prog)s -i summary.csv -d dist_path -o out.csv')
    parser.add_argument('-i', '--input', metavar='FILE',
        help='CSV summary file with Locityper genotyping results.')
    parser.add_argument('-d', '--distances', metavar='STR',
        help='Path, where locus name is replaced with `{}`. '
            'Only loci with available distances will be analyzed.')
    parser.add_argument('-o', '--output', metavar='FILE',
        help='Output CSV file.')
    parser.add_argument('--loo', action='store_true',
        help='Experiment performed in the leave-one-out setting.')
    args = parser.parse_args()

    cached_dists = {}
    with common.open(args.input) as inp, common.open(args.output, 'w') as out:
        header = None
        for line in inp:
            if not line.startswith('#'):
                header = line.strip().split('\t')
                break

        assert header is not None
        out.write(f'# {" ".join(sys.argv)}\n')
        out.write('\t'.join(header))
        out.write('\tpred_div\tpred_qv')
        if args.loo:
            out.write('\tavail_div\tavail_qv')
        out.write('\n')

        sample_col = header.index('sample')
        locus_col = header.index('locus')
        gt_col = header.index('genotype')
        for raw_line in tqdm(inp):
            raw_line = raw_line.strip()
            line = raw_line.split('\t')
            locus = line[locus_col]
            dists = load_distances(cached_dists, locus, args.distances)
            if dists is None:
                continue
            suffix = process_line(dists, line[sample_col], locus, line[gt_col], args.loo)
            out.write(f'{raw_line}\t{suffix}\n')
    missing = sum(dists is None for dists in cached_dists.values())
    if missing:
        sys.stderr.write(f'Skipping {missing} loci: no distance files found\n')


if __name__ == '__main__':
    main()
