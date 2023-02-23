#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import itertools
import functools
import collections
import sys
from scipy.stats import pmean


def main():
    parser = argparse.ArgumentParser(
        description='Calculates distance to a target tuple of haplotype to all other tuples of haplotypes')
    parser.add_argument('dist', metavar='<file>',
        help='Input CSV file with distances pairwise distances.')
    parser.add_argument('-t', '--target', metavar='<str>', required=True,
        help='Target tuple of haplotypes (separated by comma, without any spaces).')
    parser.add_argument('-m', '--mean', type=float, metavar='<float>', default=1,
        help='Power mean (Hölder mean) parameter [default: %(default)s].')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<file>', default=sys.stdout,
        help='Output CSV file (stdout by default).')
    args = parser.parse_args()

    targets = list(map(str.strip, args.target.split(',')))
    distances = pd.read_csv(args.dist, delimiter='\t')

    power = args.mean
    if power == 1:
        mean = np.mean
    elif power == np.inf:
        mean = np.max
    elif power == -np.inf:
        mean = np.min
    else:
        mean = functools.partial(pmean, p=power)

    to_target = [distances[distances.rid == target].reset_index() for target in targets]
    nrows = [df.shape[0] for df in to_target]
    assert all(nrows)
    results = collections.defaultdict(lambda: np.inf)
    for ixs in itertools.product(*map(range, nrows)):
        query = []
        scores = []
        for df, i in zip(to_target, ixs):
            query.append(df.qid[i])
            scores.append(df.dist[i])
        score = mean(scores)
        query.sort()
        query = ','.join(sorted(query))
        results[query] = min(results[query], score)

    out = args.output
    out.write('# {}\n'.format(' '.join(sys.argv)))
    out.write('target\tquery\tscore\n')
    prefix = ','.join(targets) + '\t'
    for query, score in sorted(results.items(), key=lambda t: t[1]):
        out.write('{}{}\t{:.10f}\n'.format(prefix, query, score))


if __name__ == '__main__':
    main()
