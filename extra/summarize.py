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


def process(res, sol, filtering, dist):
    sol = sol.groupby('genotype').last().reset_index()
    sol = sol.merge(dist, on='genotype').reset_index(drop=True)

    if filtering is None:
        filtering = sol
    else:
        filtering = filtering.merge(dist, on='genotype').reset_index(drop=True)
    min_dist = filtering.dist.min()
    closest_gts = list(filtering[filtering.dist == min_dist].genotype)

    sol.set_index('genotype', inplace=True)
    s = '{}\t{}\t{:.10f}\t{:.10f}\t'.format(filtering.shape[0], sol.shape[0],
        pearsonr(sol.lik, sol.dist).statistic, spearmanr(sol.lik, sol.dist).statistic)

    closest_gt_liks = [sol.loc[gt].lik if gt in sol.index else np.nan for gt in closest_gts]
    s += '{}\t{}\t{}\t'.format(min_dist, ';'.join(closest_gts), ';'.join(map('{:.3f}'.format, closest_gt_liks)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        closest_lik = np.nanmin(closest_gt_liks)
    closest_rank = sol.shape[0] if np.isnan(closest_lik) else (sol.lik > closest_lik).sum() + 1
    s += f'{closest_rank}\t'

    ml_gt = res['genotype']
    ml_dist = sol.loc[ml_gt].dist
    ml_dist_rank = (filtering.dist < ml_dist).sum() + 1
    s += '{}\t{}\t{}\t{}'.format(ml_gt, sol.loc[ml_gt].dist, sol.loc[ml_gt].lik, ml_dist_rank)

    weighted_dist = 0.0
    weight_sum = 0.0
    for option in res['options']:
        w = option['prob']
        weighted_dist += sol.loc[option['genotype']].dist * w
        weight_sum += w
    weighted_dist /= weight_sum

    s += '\t{:.5f}\t{:.5}\t{:.10f}\n'.format(res['quality'], sol.loc[ml_gt].depth0_frac, weighted_dist)
    return s


def load_and_process(tup, tags, input_fmt, dist_fmt):
    d = dict(zip(tags, tup))
    input_dir = input_fmt.format(**d)
    with gzip.open(os.path.join(input_dir, 'res.json.gz'), 'rt') as f:
        res = json.load(f)
    sol = pd.read_csv(os.path.join(input_dir, 'sol.csv.gz'), sep='\t', comment='#')
    try:
        filtering = pd.read_csv(os.path.join(input_dir, 'filter.csv.gz'), sep='\t', comment='#')
    except FileNotFoundError:
        filtering = None
    dist = pd.read_csv(dist_fmt.format(**d), sep='\t', comment='#')
    s = process(res, sol, filtering, dist)
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
        'ml_gt\tml_gt_dist\tml_gt_lik\tml_dist_rank\tgt_qual\tdepth0_frac\tweighted_dist\n')

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
