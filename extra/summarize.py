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
import glob

import common


def load_tags(path):
    matches = list(re.finditer(r'\{([a-zA-Z0-9_]+)\}', path))
    tags = []
    filename_pattern = '^'
    glob_pattern = ''
    s = 0
    for m in matches:
        tags.append(m.group(1))
        filename_pattern += re.escape(path[s:m.start()]) + r'([a-zA-Z0-9_+-]+)'
        glob_pattern += path[s:m.start()] + '*'
        s = m.end()
    filename_pattern += re.escape(path[s:]) + '$'
    glob_pattern += path[s:]

    filename_pattern = re.compile(filename_pattern)
    all_tuples = []
    sys.stderr.write(f'Searching for files/directories `{glob_pattern}`\n')
    for filename in tqdm(glob.iglob(glob_pattern)):
        m = filename_pattern.match(filename)
        if m is not None:
            all_tuples.append(tuple(m.groups()))
    return tags, all_tuples


def discard_tags(tags, all_tuples, path2):
    sel_tuples = []
    disc_tuples = []
    for tup in all_tuples:
        d = dict(zip(tags, tup))
        if os.path.exists(path2.format(**d)):
            sel_tuples.append(tup)
        else:
            disc_tuples.append(tup)

    sys.stderr.write(f'Found {len(all_tuples)} tag combinations. Discarded {len(disc_tuples)} of them\n')
    return sel_tuples, disc_tuples


def process(prefix, res, sol, filt, dist):
    if sol is None:
        sol = pd.DataFrame(columns=['genotype', 'lik'])
    else:
        max_stage = sol.stage.max()
        sol = sol[sol.stage == max_stage].groupby('genotype').mean()

    if filt is not None:
        filt.set_index('genotype', inplace=True)
        filt_gts = list(filt.index)
        highest_score = filt.score.max()
    else:
        filt_gts = list(sol.index)

    true_gts = dist[dist.distance == 0].genotype
    dist.set_index('genotype', inplace=True)

    pred_gt = res['genotype']
    reached_row = dist.loc[pred_gt]
    reached_dist = reached_row.distance
    reached_div = reached_row.divergence
    weighted_dist = 0.0
    weighted_div = 0.0
    weight_sum = 0.0
    gt_probs = {}
    for option in res.get('options', ()):
        gt = option['genotype']
        w = option['prob']
        if w is None:
            weighted_dist = np.nan
            break
        row = dist.loc[gt]
        weighted_dist += row.distance * w
        weighted_div += row.divergence * w
        weight_sum += w
        gt_probs[gt] = option['log10_prob']
        if gt_probs[gt] is None:
            gt_probs[gt] = -np.inf
    if weight_sum > 0:
        weighted_dist /= weight_sum
        weighted_div /= weight_sum
    else:
        weighted_dist = np.nan
        weighted_div = np.nan

    x = sol.lik
    y = dist.loc[sol.index].distance
    try:
        pearson = pearsonr(x, y).statistic
    except ValueError:
        pearson = np.nan
    try:
        spearman = spearmanr(x, y).statistic
    except ValueError:
        spearman = np.nan

    avail_dist = np.inf
    avail_div = np.inf
    for gt in filt_gts:
        row = dist.loc[gt]
        if row.distance < avail_dist:
            avail_dist = row.distance
            avail_div = row.divergence
    weighted_dist = '{:.5f}'.format(weighted_dist).rstrip('0').rstrip('.')
    warnings = '*' if 'warnings' not in res else ','.join(sorted(res['warnings']))
    s1 = '{}{}\t{}\t{:.5f}\t{:.5f}\t{:.0f}\t{:.7f}\t{:.0f}\t{:.7f}\t{:s}\t{:.7f}\t{:.3f}\t{}\n'.format(
        prefix, filt.shape[0] if filt is not None else sol.shape[0], sol.shape[0], pearson, spearman,
        avail_dist, avail_div, reached_dist, reached_div, weighted_dist, weighted_div, float(res['quality']), warnings)

    highest_lik = sol.lik.max()
    interesting_gts = set(true_gts)
    interesting_gts.update(dist.index[dist.distance == avail_dist])
    interesting_gts.update(gt_probs.keys())
    s2 = ''

    for gt in interesting_gts:
        d = dist.loc[gt].distance
        s2 += '{}{}\t{}\t{}\t{:.10f}\t'.format(prefix, gt, d, 1 + np.sum(dist.distance < d), dist.loc[gt].divergence)
        if filt is not None and gt in filt.index:
            score = filt.loc[gt].score
            s2 += '{:.5f}\t{:.5f}\t{}\t'.format(score, score - highest_score, 1 + np.sum(filt.score > score))
        else:
            s2 += 'nan\tnan\tnan\t'

        if gt in sol.index:
            lik = sol.loc[gt].lik
            s2 += '{:.5f}\t{:.5f}\t{}\t'.format(lik, lik - highest_lik, 1 + np.sum(sol.lik > lik))
        else:
            s2 += 'nan\tnan\tnan\t'

        s2 += '{:.7g}\n'.format(gt_probs.get(gt, np.nan))
    return s1, s2


def load_and_process(tup, tags, input_fmt, dist_fmt):
    d = dict(zip(tags, tup))
    input_dir = input_fmt.format(**d)
    with gzip.open(os.path.join(input_dir, 'res.json.gz'), 'rt') as f:
        res = json.load(f)
    try:
        sol = pd.read_csv(os.path.join(input_dir, 'sol.csv.gz'), sep='\t', comment='#')
    except FileNotFoundError:
        sol = None
    try:
        filtering = pd.read_csv(os.path.join(input_dir, 'filter.csv.gz'), sep='\t', comment='#')
    except FileNotFoundError:
        filtering = None
    dist = pd.read_csv(dist_fmt.format(**d), sep='\t', comment='#')
    prefix = ''.join(map('{}\t'.format, tup))
    try:
        return process(prefix, res, sol, filtering, dist)
    except:
        sys.stderr.write(f'Encountered problem at directory {input_dir}\n')
        raise


def main():
    parser = argparse.ArgumentParser(
        description='Summarize genotyping result.',
        usage='%(prog)s -i path -d path -o out.csv [-@ threads]')
    parser.add_argument('-i', '--input', metavar='STR', required=True,
        help='Path to genotyping results. Tags within `{tag}` are automatically found. '
            'Input directories must contain `sol.csv.gz`, `filter.csv.gz` and `res.json.gz` files.')
    parser.add_argument('-d', '--distances', metavar='STR',  required=True,
        help='Path to distances. Tags within `{tag}` are automatically found.')
    parser.add_argument('-o', '--output', metavar='STR',  required=True,
        help='Output prefix.')
    parser.add_argument('-@', '--threads', metavar='INT', type=int, default=8,
        help='Number of threads [%(default)s].')
    args = parser.parse_args()

    path_prefix = args.output + ('/' if os.path.isdir(args.output) else '.')
    tags, tag_tuples = load_tags(args.input)
    tag_tuples, disc_tuples = discard_tags(tags, tag_tuples, args.distances)
    tags_prefix = ''.join(map('{}\t'.format, tags))

    if disc_tuples:
        with common.open(f'{path_prefix}missing.csv.gz', 'w') as out:
            out.write('\t'.join(tags) + '\n')
            for tup in disc_tuples:
                out.write('\t'.join(tup) + '\n')
    if not tag_tuples:
        sys.stderr.write('No input files found, stopping.\n')
        exit(1)

    out_summary = common.open(f'{path_prefix}summary.csv.gz', 'w')
    out_summary.write('# {}\n'.format(' '.join(sys.argv)))
    out_summary.write(f'{tags_prefix}avail_gts\tafter_filt_gts\tpearsonr\tspearmanr\t'
        'avail_dist\tavail_div\tachieved_dist\tachieved_div\tweighted_dist\tweighted_div\tquality\twarnings\n')

    out_gts = common.open(f'{path_prefix}gts.csv.gz', 'w')
    out_gts.write('# {}\n'.format(' '.join(sys.argv)))
    out_gts.write(f'{tags_prefix}genotype\tdist\tdist_rank\tdivergence\tfilt_score\tscore_diff\t'
        'score_rank\tlik\tlik_diff\tlik_rank\tprob\n')

    n = len(tag_tuples)
    threads = min(n, args.threads)
    f = functools.partial(load_and_process, tags=tags, input_fmt=args.input, dist_fmt=args.distances)
    pbar = tqdm if n > 1 else lambda x, *args, **kwargs: x

    if threads > 1:
        with Pool(threads) as pool:
            for s1, s2 in pbar(pool.imap(f, tag_tuples), total=n):
                out_summary.write(s1)
                out_gts.write(s2)
    else:
        for s1, s2 in pbar(map(f, tag_tuples), total=n):
            out_summary.write(s1)
            out_gts.write(s2)


if __name__ == '__main__':
    main()
