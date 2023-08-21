#!/usr/bin/env python3

import sys
import argparse
from multiprocessing import Pool
import itertools
import numpy as np
import pandas as pd
import functools
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm

from common import open_stream


def load_tags(path):
    tags = []
    tag_vals = []
    if path is not None:
        with open_stream(path) as f:
            for line in f:
                line_split = line.strip().split()
                tags.append(line_split[0].rstrip(':'))
                tag_vals.append(tuple(line_split[1:]))
                assert tags[-1] and tag_vals[-1]
    tag_tuples = list(itertools.product(*tag_vals))
    return tuple(tags), tag_tuples


def process(liks, dist):
    liks = liks.groupby('genotype').last().reset_index()
    merged = liks.merge(dist, on='genotype').reset_index(drop=True)
    total_gts = merged.shape[0]
    merged = merged[merged.stage > 0].reset_index(drop=True)
    merged = merged.sort_values(['dist', 'lik'], ascending=[True, False]).reset_index(drop=True)
    top_lik_ix = np.argmax(merged.lik)

    s = '{}\t{}\t{:.10}\t{:.10}'.format(total_gts, merged.shape[0],
        pearsonr(merged.lik, merged.dist).statistic,
        spearmanr(merged.lik, merged.dist).statistic)
    s += '\t{}\t{}\t{}\t{}'.format(merged.genotype[0], merged.dist[0], merged.lik[0],
        (merged.lik > merged.lik[0]).sum() + 1)
    i = np.argmax(merged.lik)
    s += '\t{}\t{}\t{}\t{}\n'.format(merged.genotype[i], merged.dist[i], merged.lik[i], i + 1)
    return s


def load_and_process(tup, tags, lik_fmt, dist_fmt):
    d = dict(zip(tags, tup))
    liks = pd.read_csv(lik_fmt.format(**d), sep='\t', comment='#')
    dist = pd.read_csv(dist_fmt.format(**d), sep='\t', comment='#')
    s = process(liks, dist)
    return '\t'.join(tup + (s,))


def main():
    parser = argparse.ArgumentParser(
        description='Summarize genotyping result.',
        usage='%(prog)s -i dir -d dist.csv')
    parser.add_argument('-g', '--genotyping', metavar='STR',
        help='Path to genotyping results. Parts like `{tag}` are replaced, tags are supplied separately. '
            'Input directories must contain `lik.csv.gz` and `res.json.gz` files.')
    parser.add_argument('-d', '--distances', metavar='FILE',
        help='Path to distances. Parts like `{tag}` are replaced.')
    parser.add_argument('-t', '--tags', metavar='FILE',
        help='Optional: file, where each line starts with `tag:`, and contains all possible values through whitespace.')
    parser.add_argument('-o', '--output', metavar='FILE',
        help='Output CSV file.')
    parser.add_argument('-@', '--threads', metavar='INT', type=int, default=8,
        help='Number of threads [%(default)s].')
    args = parser.parse_args()

    tags, tag_tuples = load_tags(args.tags)

    out = open_stream(args.output, 'w')
    out.write('# {}\n'.format(' '.join(sys.argv)))
    for tag in tags:
        out.write(tag + '\t')
    out.write('total_gts\test_gts\tpearsonr\tspearmanr\tclosest_gt\tclosest_gt_dist\tclosest_gt_lik\tclosest_lik_rank'
        '\tml_gt\tml_gt_dist\tml_gt_lik\tml_dist_rank\n')

    lik_fmt = args.genotyping + '/lik.csv.gz'
    dist_fmt = args.distances

    n = len(tag_tuples)
    threads = min(n, args.threads)
    f = functools.partial(load_and_process, tags=tags, lik_fmt=lik_fmt, dist_fmt=dist_fmt)

    if threads > 1:
        with Pool(threads) as pool:
            for line in tqdm(pool.imap(f, tag_tuples), total=n):
                out.write(line)
    else:
        for line in tqdm(map(f, tag_tuples), total=n):
            out.write(line)


if __name__ == '__main__':
    main()
