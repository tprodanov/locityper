#!/usr/bin/env python3

import sys
import os
import argparse
import warnings
from multiprocessing import Pool
import itertools
import numpy as np
import pandas as pd
import functools
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm
import re
import gzip
import json

from common import open_stream


def _recursive_find_tuples(path, matches, all_tuples, curr_tuple=(), depth=0, shift=0):
    if depth == len(matches):
        all_tuples.append(curr_tuple)
        return 0
    m = matches[depth]
    dir = path[ : m.start() + shift]
    if not os.path.isdir(dir):
        return 1

    suffix = path[m.end() + shift : ]
    subdirs = []
    with os.scandir(dir) as it:
        for entry in it:
            if not entry.name.startswith('.') and entry.is_dir():
                subdirs.append(entry)

    not_found = 0
    for entry in subdirs:
        not_found += _recursive_find_tuples(entry.path + suffix, matches, all_tuples,
            curr_tuple + (entry.name,), depth + 1, len(entry.path) - m.end())
    return not_found


def load_tags(path1, path2):
    path1 = os.path.abspath(path1)
    matches = list(re.finditer(r'\{([a-zA-Z0-9_]+)\}', path1))
    tags = [m.group(1) for m in matches]
    all_tuples = []
    not_found = _recursive_find_tuples(path1, matches, all_tuples)
    if not_found:
        sys.stderr.write(f'Skipped {not_found} directories\n')

    filt_tuples = []
    for tup in all_tuples:
        d = dict(zip(tags, tup))
        if os.path.exists(path2.format(**d)):
            filt_tuples.append(tup)

    sys.stderr.write('Found {} tag combinations. Discarded {} of them\n'.format(
        len(all_tuples), len(all_tuples) - len(filt_tuples)))
    return tags, filt_tuples


def process(res, liks, filter_liks, dist):
    liks = liks.groupby('genotype').last().reset_index()
    liks = liks.merge(dist, on='genotype').reset_index(drop=True)

    if filter_liks is None:
        filter_liks = liks
    else:
        filter_liks = filter_liks.merge(dist, on='genotype').reset_index(drop=True)
    min_dist = filter_liks.dist.min()
    closest_gts = list(filter_liks[filter_liks.dist == min_dist].genotype)

    liks.set_index('genotype', inplace=True)
    s = '{}\t{}\t{:.10f}\t{:.10f}\t'.format(filter_liks.shape[0], liks.shape[0],
        pearsonr(liks.lik, liks.dist).statistic, spearmanr(liks.lik, liks.dist).statistic)

    closest_gt_liks = [liks.loc[gt].lik if gt in liks.index else np.nan for gt in closest_gts]
    s += '{}\t{}\t{}\t'.format(min_dist, ';'.join(closest_gts), ';'.join(map('{:.3f}'.format, closest_gt_liks)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        closest_lik = np.nanmin(closest_gt_liks)
    closest_rank = liks.shape[0] if np.isnan(closest_lik) else (liks.lik > closest_lik).sum() + 1
    s += f'{closest_rank}\t'

    ml_gt = res['genotype']
    ml_dist = liks.loc[ml_gt].dist
    ml_dist_rank = (filter_liks.dist < ml_dist).sum() + 1
    s += '{}\t{}\t{}\t{}'.format(ml_gt, liks.loc[ml_gt].dist, liks.loc[ml_gt].lik, ml_dist_rank)

    weighted_dist = 0.0
    weight_sum = 0.0
    for option in res['options']:
        w = option['prob']
        weighted_dist += liks.loc[option['genotype']].dist * w
        weight_sum += w
    weighted_dist /= weight_sum

    s += '\t{:.5f}\t{:.10f}\n'.format(res['quality'], weighted_dist)
    return s


def load_and_process(tup, tags, input_fmt, dist_fmt):
    d = dict(zip(tags, tup))
    input_dir = input_fmt.format(**d)
    with gzip.open(os.path.join(input_dir, 'res.json.gz'), 'rt') as f:
        res = json.load(f)
    liks = pd.read_csv(os.path.join(input_dir, 'lik.csv.gz'), sep='\t', comment='#')
    try:
        filter_liks = pd.read_csv(os.path.join(input_dir, 'filter.csv.gz'), sep='\t', comment='#')
    except FileNotFoundError:
        filter_liks = None
    dist = pd.read_csv(dist_fmt.format(**d), sep='\t', comment='#')
    s = process(res, liks, filter_liks, dist)
    return '\t'.join(tup + (s,))


def main():
    parser = argparse.ArgumentParser(
        description='Summarize genotyping result.',
        usage='%(prog)s -i path -d path -o out.csv [-@ threads]')
    parser.add_argument('-i', '--input', metavar='STR', required=True,
        help='Path to genotyping results. Tags within `{tag}` are automatically found. '
            'Input directories must contain `lik.csv.gz` and `res.json.gz` files.')
    parser.add_argument('-d', '--distances', metavar='STR',  required=True,
        help='Path to distances. Tags within `{tag}` are automatically found.')
    parser.add_argument('-o', '--output', metavar='FILE',  required=True,
        help='Output CSV file.')
    parser.add_argument('-@', '--threads', metavar='INT', type=int, default=8,
        help='Number of threads [%(default)s].')
    args = parser.parse_args()

    tags, tag_tuples = load_tags(args.input, args.distances)

    out = open_stream(args.output, 'w')
    out.write('# {}\n'.format(' '.join(sys.argv)))
    for tag in tags:
        out.write(tag + '\t')
    out.write('total_gts\trem_gts\tpearsonr\tspearmanr\tsmallest_dist\tclosest_gts\tclosest_gt_liks\tclosest_rank\t'
        'ml_gt\tml_gt_dist\tml_gt_lik\tml_dist_rank\tgt_qual\tweighted_dist\n')

    n = len(tag_tuples)
    threads = min(n, args.threads)
    f = functools.partial(load_and_process, tags=tags, input_fmt=args.input, dist_fmt=args.distances)
    pbar = tqdm if n > 1 else lambda x, *args, **kwargs: x

    if threads > 1:
        with Pool(threads) as pool:
            for line in pbar(pool.imap(f, tag_tuples), total=n):
                out.write(line)
    else:
        for line in pbar(map(f, tag_tuples), total=n):
            out.write(line)


if __name__ == '__main__':
    main()
