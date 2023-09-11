#!/usr/bin/env python3

import sys
import itertools
import collections
import operator
import argparse
import numpy as np

from common import open_stream


def averaging_function(mode):
    mode = mode.lower()
    if mode == 'min':
        return np.min
    elif mode == 'max':
        return np.max
    elif mode == 'sum':
        return np.sum
    elif mode == 'mean':
        return np.mean

    power = float(mode)
    if power == 1:
        return np.mean
    elif power == np.inf:
        return np.max
    elif power == -np.inf:
        return np.min

    import scipy
    import functools
    return functools.partial(scipy.stats.pmean, p=power)


class Contigs:
    def __init__(self):
        self.contigs = []
        self.order = {}

    def add_if_needed(self, contig):
        if contig not in self.order:
            self.order[contig] = len(self.contigs)
            self.contigs.append(contig)

    def __len__(self):
        return len(self.contigs)

    def __getitem__(self, i):
        return self.contigs[i]

    def fmt_gt(self, ixs):
        return ','.join(self.contigs[i] for i in ixs)


def _load_distances(f, get_info, targets, verbose):
    targets = { target: i for i, target in enumerate(targets) }
    contigs = Contigs()
    uns_distances = []

    if verbose:
        sys.stderr.write('Load distances\n')
    for record in f:
        contig1, contig2, dist = get_info(record)
        # Second contig is reference, therefore goes first.
        contigs.add_if_needed(contig2)
        contigs.add_if_needed(contig1)
        if contig1 in targets or contig2 in targets:
            uns_distances.append((contig1, contig2, dist))
            if verbose:
                sys.stderr.write(f'    {contig1} {contig2} -> {dist:.10f}\n')
    if not uns_distances:
        sys.stderr.write(f'    ERROR: Target genotype not found\n')
        exit(1)

    n = len(targets)
    m = len(contigs)
    distances = np.zeros((n, m))
    for contig1, contig2, dist in uns_distances:
        i = targets.get(contig1)
        j = targets.get(contig2)
        if i is not None:
            distances[i, contigs.order[contig2]] = dist
        if j is not None:
            distances[j, contigs.order[contig1]] = dist
    return contigs, distances


def _process_distances(targets, contigs, distances, mean, verbose):
    n, m = distances.shape
    best_dist = collections.defaultdict(lambda: np.inf)
    range_n = np.arange(n)
    if verbose:
        sys.stderr.write('Process distances\n')
    for ixs in itertools.product(range(m), repeat=n):
        curr_dist = distances[range_n, ixs]
        dist = mean(curr_dist)
        query = tuple(sorted(ixs))
        old_val = best_dist[query]
        if verbose:
            sys.stderr.write('    {}: {} -> {:.10f} (old {:.10f})\n'
                .format(contigs.fmt_gt(query), curr_dist, dist, old_val))
        best_dist[query] = min(old_val, dist)
    return best_dist


def _get_info_paf(tag):
    prefix = tag + ':'
    def inner(line):
        line = line.strip().split('\t')
        contig1 = line[0]
        contig2 = line[5]
        try:
            dist = float(next(col.rsplit(':', 1)[1] for col in line[12:] if col.startswith(prefix)))
        except StopIteration:
            sys.stderr.write(f'Distance not available for contigs {contig1} and {contig2}\n')
            exit(1)
        return contig1, contig2, dist
    return inner


def _get_info_bam(tag):
    def inner(record):
        return record.query_name, record.reference_name, record.get_tag(tag)
    return inner


def main():
    parser = argparse.ArgumentParser(
        description='Calculates distance between a target genotype and all other genotypes',
        usage='%(prog)s alns.paf target -o out.csv [-f field -a averaging]')
    parser.add_argument('-i', '--input', metavar='FILE',
        help='Input PAF[.gz] or BAM file with pairwise distances.')
    parser.add_argument('-g', '--genotype', metavar='STR',
        help='Target genotype (through comma).')
    parser.add_argument('-o', '--output', metavar='FILE', required=False,
        help='Output CSV file [stdout].')
    parser.add_argument('-f', '--field', metavar='STR', default='dv',
        help='Take haplotype distance from this field [%(default)s].')
    parser.add_argument('-a', '--averaging', metavar='STR', default='mean',
        help='Averaging function [%(default)s]. '
            'Possible values: `min`, `max`, `sum`, `mean`, or FLOAT: Hölder mean.')
    parser.add_argument('-v', '--verbose', action='store_true',
        help='Produce more output')
    args = parser.parse_args()

    targets = list(map(str.strip, args.genotype.split(',')))
    if args.input.endswith('.bam'):
        import pysam
        save = pysam.set_verbosity(0)
        with pysam.AlignmentFile(args.input, require_index=False) as f:
            pysam.set_verbosity(save)
            contigs, distances = _load_distances(f, _get_info_bam(args.field), targets, args.verbose)
    else:
        with open_stream(args.input) as f:
            contigs, distances = _load_distances(f, _get_info_paf(args.field), targets, args.verbose)

    aver = averaging_function(args.averaging)
    best_dist = _process_distances(targets, contigs, distances, aver, args.verbose)

    with open_stream(args.output, 'w') as out:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        out.write('# target: {}\n'.format(args.genotype))
        out.write('# field: {}\n'.format(args.field))
        out.write('# averaging: {}\n'.format(args.averaging))
        out.write('genotype\tdist\n')
        for query, dist in sorted(best_dist.items(), key=operator.itemgetter(1)):
            out.write('{}\t{:.10g}\n'.format(contigs.fmt_gt(query), dist))


if __name__ == '__main__':
    main()
