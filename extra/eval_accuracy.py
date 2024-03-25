#!/usr/bin/env python3

import numpy as np
import argparse
import sys
import os
import csv
import collections
from tqdm import tqdm

import common
import gt_dist


def load_distances(cached_dists, locus, paf_fmt, disc_fmt):
    if locus in cached_dists:
        return cached_dists[locus]

    paf_filename = paf_fmt.format(locus)
    if not os.path.exists(paf_filename):
        dists = None
    else:
        dists = gt_dist.Distances(disc_fmt.format(locus), paf_filename)
    cached_dists[locus] = dists
    return dists


def get_qv(val):
    return np.inf if val == 0 else -10 * np.log10(val)


def process_line(distances, sample, locus, genotype, loo):
    target = distances.get_sample_haplotypes(sample)
    if len(target) != 2:
        return '\t*\tNA\tNA' * (1 + loo)

    query = genotype.split(',')
    edit, size, div = distances.calc_distance(target, query)
    s = f'\t{edit}/{size}\t{div:.9f}\t{get_qv(div):.9f}'
    if loo:
        _, edit, size, div = distances.find_closest_loo(target)
        s += f'\t{edit}/{size}\t{div:.9f}\t{get_qv(div):.9f}'
    return s


def main():
    parser = argparse.ArgumentParser(
        description='Evaluating genotyping accuracy',
        usage='%(prog)s -i summary.csv -a paf_path -d discarded_path -o out.csv [--loo]')
    parser.add_argument('-i', '--input', metavar='FILE',
        help='CSV summary file with Locityper genotyping results.')
    parser.add_argument('-a', '--alignments', metavar='STR',
        help='Path to PAF files, where locus name is replaced with `{}`. '
            'Only loci with available distances will be analyzed.')
    parser.add_argument('-d', '--discarded', metavar='STR',
        help='Path to files with discarded haplotypes, where locus name is replaced with `{}`.')
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
        out.write('\tpred_edit\tpred_div\tpred_qv')
        if args.loo:
            out.write('\tavail_edit\tavail_div\tavail_qv')
        out.write('\n')

        sample_col = header.index('sample')
        locus_col = header.index('locus')
        gt_col = header.index('genotype')
        for raw_line in tqdm(inp):
            raw_line = raw_line.strip()
            line = raw_line.split('\t')
            locus = line[locus_col]
            dists = load_distances(cached_dists, locus, args.alignments, args.discarded)
            if dists is None:
                continue
            suffix = process_line(dists, line[sample_col], locus, line[gt_col], args.loo)
            out.write(f'{raw_line}{suffix}\n')
    missing = sum(dists is None for dists in cached_dists.values())
    if missing:
        sys.stderr.write(f'Skipping {missing} loci: no distance files found\n')


if __name__ == '__main__':
    main()
