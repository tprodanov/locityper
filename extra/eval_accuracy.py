#!/usr/bin/env python3

import numpy as np
import argparse
import sys
import os
import csv
import collections
import itertools
from tqdm import tqdm

import common
import gt_dist


def load_distances(cached_dists, locus, excluded_samples, paf_fmt, db_dir):
    if locus in cached_dists:
        return cached_dists[locus]

    paf_filename = paf_fmt.format(locus)
    locus_dir = f'{db_dir}/loci/{locus}'
    excl_haps = ()
    if not os.path.exists(paf_filename):
        dists = None
    elif not os.path.exists(locus_dir):
        sys.stderr.write(f'Database directory `{locus_dir}` does not exist.\n')
        exit(1)
    else:
        dists = gt_dist.Distances(os.path.join(locus_dir, 'discarded_haplotypes.txt'), paf_filename)
        if excluded_samples:
            excl_haps = set()
            for sample in excluded_samples:
                excl_haps.update(dists.get_sample_haplotypes(sample))

    cached_dists[locus] = dists, excl_haps
    return dists, excl_haps


def get_qv(val):
    return np.inf if val == 0 else -10 * np.log10(val)


def process_line(line, distances, excl_haps, sample, locus, genotype, loo, out):
    target = distances.get_sample_haplotypes(sample)
    if len(target) != 2:
        target = tuple(target) + (None,) * (2 - len(target))

    query = genotype.split(',')
    query_dists = distances.calc_distance(target, query).iter_strs()
    ext_out = bool(loo or excl_haps)
    if ext_out:
        loo_dists = distances.find_closest(target, loo, excl_haps)[1].iter_strs()
    else:
        loo_dists = itertools.repeat(None)

    for i, (s1, s2) in enumerate(zip(query_dists, loo_dists), 1):
        curr_line = f'{line}\t{"hap" if i <= 2 else "gt"}\t{s1}'
        if ext_out:
            curr_line += f'\t{s2}\n'
        else:
            curr_line += '\n'
        out.write(curr_line)


def main():
    parser = argparse.ArgumentParser(
        description='Evaluating genotyping accuracy',
        usage='%(prog)s -i summary.csv -a paf_path -d discarded_path -o out.csv [--loo]')
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='CSV summary file with Locityper genotyping results.')
    parser.add_argument('-d', '--database', metavar='DIR', required=True,
        help='Path to the database.')
    parser.add_argument('-a', '--alignments', metavar='STR', required=True,
        help='Path to PAF files, where locus name is replaced with `{}`. '
            'Only loci with available distances will be analyzed.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    parser.add_argument('--loo', action='store_true',
        help='Experiment performed in the leave-one-out setting.')
    parser.add_argument('-e', '--excluded', metavar='FILE',
        help='File with sample names, which were excluded from the analysis.')
    args = parser.parse_args()

    if not os.path.exists(args.database):
        sys.stderr.write(f'Database `{args.database}` does not exist\n')
        exit(1)

    excluded_samples = None
    if args.excluded:
        with common.open(args.excluded) as inp:
            excluded_samples = [line.strip() for line in inp]
        sys.stderr.write(f'Exclude {len(excluded_samples)} samples\n')

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
        out.write('\tquery_type\tedit\tsize\tdiv\tqv')
        if args.loo or excluded_samples is not None:
            out.write('\tavail_edit\tavail_size\tavail_div\tavail_qv')
        out.write('\n')

        sample_col = header.index('sample')
        locus_col = header.index('locus')
        gt_col = header.index('genotype')
        for raw_line in tqdm(inp):
            raw_line = raw_line.strip()
            line = raw_line.split('\t')
            locus = line[locus_col]
            dists, excl_haps = load_distances(cached_dists, locus, excluded_samples, args.alignments, args.database)
            if dists is None:
                continue
            process_line(raw_line, dists, excl_haps, line[sample_col], locus, line[gt_col], args.loo, out)
    missing = sum(dists is None for dists in cached_dists.values())
    if missing:
        sys.stderr.write(f'Skipping {missing} loci: no distance files found\n')


if __name__ == '__main__':
    main()
