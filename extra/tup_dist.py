#!/usr/bin/env python3

import argparse
import numpy as np
import itertools
import functools
import collections
import sys
import gzip


def _load_distances(f, tag, targets):
    prefix = f'{tag}:'
    distances = { target: { target: 0.0 } for target in targets }
    for line in f:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        contig1 = line[0]
        contig2 = line[5]
        need_contig1 = contig1 in distances
        need_contig2 = contig2 in distances
        if not need_contig1 and not need_contig2:
            continue

        dist = None
        for col in line[12:]:
            if col.startswith(prefix):
                dist = float(col.split(':', 2)[2])
                break
        else:
            raise RuntimeError(f'Cannot find tag {tag} in the alignment between {contig1} and {contig2}')

        if need_contig1:
            distances[contig1][contig2] = dist
        if need_contig2:
            distances[contig2][contig1] = dist
    # Convert Dict<str, Dict<str, float>> into List<List<(str, float)>>.
    return [list(distances[target].items()) for target in targets]


def _create_power_mean(power):
    if power == 1:
        return np.mean
    elif power == np.inf:
        return np.max
    elif power == -np.inf:
        return np.min

    import scipy
    return functools.partial(scipy.stats.pmean, p=power)


def open_stream(filename, mode='r'):
    assert mode == 'r' or mode == 'w'
    if filename is None or filename == '-':
        return sys.stdin if mode == 'r' else sys.stdout
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    else:
        return open(filename, mode)


def main():
    parser = argparse.ArgumentParser(
        description='Calculates distance between a target haplotypes tuple and all other haplotype tuples',
        usage='%(prog)s alns.paf target -o out.csv [-t tag -p power]')
    parser.add_argument('paf', metavar='FILE',
        help='Input PAF[.gz] file with pairwise distances.')
    parser.add_argument('target', metavar='STR',
        help='Target tuple of haplotypes (separated by comma, without spaces).')
    parser.add_argument('-o', '--output', metavar='FILE', required=False,
        help='Output CSV file [stdout].')
    parser.add_argument('-t', '--tag', metavar='STR', default='dv',
        help='Take haplotype distance from this field [%(default)s].')
    parser.add_argument('-p', '--power', type=float, metavar='FLOAT', default=2,
        help='Power mean (Hölder mean) parameter [%(default)s].')
    args = parser.parse_args()

    mean = _create_power_mean(args.power)
    targets = list(map(str.strip, args.target.split(',')))
    with open_stream(args.paf) as f:
        distances = _load_distances(f, args.tag, targets)
    assert min(map(len, distances)) > 1
    counts = list(map(len, distances))
    assert all(counts)

    results = collections.defaultdict(lambda: np.inf)
    for ixs in itertools.product(*map(range, counts)):
        query = []
        scores = []
        for target_distances, i in zip(distances, ixs):
            query.append(target_distances[i][0])
            scores.append(target_distances[i][1])
        score = mean(scores)
        query.sort()
        query = ','.join(sorted(query))
        results[query] = min(results[query], score)

    with open_stream(args.output, 'w') as out:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        out.write('# target: {}\n'.format(args.target))
        out.write('# tag: {}\n'.format(args.tag))
        out.write('# power: {}\n'.format(args.power))
        out.write('query\tdist\n')
        for query, score in sorted(results.items(), key=lambda t: t[1]):
            out.write('{}\t{:.10f}\n'.format(query, score))


if __name__ == '__main__':
    main()
