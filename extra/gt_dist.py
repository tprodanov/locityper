#!/usr/bin/env python3

import sys
import itertools
import operator
import argparse
import numpy as np
import os
import pysam

import common


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


def _load_distances(f, targets, verbose):
    targets_dict = { target: i for i, target in enumerate(targets) }
    contigs = Contigs()
    uns_distances = []

    if verbose:
        sys.stderr.write('Loading distances\n')
    for record in f:
        contig1 = record.query_name
        contig2 = record.reference_name
        dist = record.get_tag('NM')
        div = record.get_tag('dv')
        if div != 0:
            aln_len = int(round(dist / div))
        else:
            cigar = record.cigartuples
            assert len(cigar) == 1
            aln_len = cigar[0][1]

        # Second contig is reference, therefore goes first.
        contigs.add_if_needed(contig2)
        contigs.add_if_needed(contig1)
        if contig1 in targets_dict:
            uns_distances.append((contig1, contig2, dist, aln_len))
            if verbose:
                sys.stderr.write(f'    {contig1} {contig2}: NM {dist}, div {div:.6g}, aln len {aln_len} \n')
    if not uns_distances:
        sys.stderr.write('    ERROR: Target genotype {} not found\n'.format(','.join(targets)))
        exit(1)

    n = len(targets)
    m = len(contigs)
    distances = np.full((n, m), np.nan)
    aln_lens = np.full((n, m), np.nan)
    for contig1, contig2, dist, aln_len in uns_distances:
        i = targets_dict.get(contig1)
        if i is not None:
            distances[i, contigs.order[contig2]] = float(dist)
            aln_lens[i, contigs.order[contig2]] = float(aln_len)
    na_ixs = np.where(np.isnan(distances))
    if len(na_ixs[0]):
        sys.stderr.write('Distances unavailable for {} haplotype pairs, for example for {} and {}\n'.format(
            len(na_ixs[0]), targets[na_ixs[0][0]], contigs[na_ixs[1][0]]))
        exit(1)
    return contigs, distances, aln_lens


def _process_distances(ploidy, contigs, distances, aln_lens, verbose):
    n, m = distances.shape
    best_dist = {}
    range_n = np.arange(n)
    if verbose:
        sys.stderr.write('Process distances\n')
    for ixs in itertools.product(range(m), repeat=ploidy):
        curr_dist = distances[range_n, ixs]
        dist = np.sum(curr_dist)
        div = dist / np.sum(aln_lens[range_n, ixs])
        query = tuple(sorted(ixs))
        old_dist, old_div = best_dist.get(query, (np.inf, np.inf))
        if verbose:
            sys.stderr.write('    {}: {} -> {:.0f}, div {:.6g} (old {:.0f}, {:.6g})\n'
                .format(contigs.fmt_gt(query), curr_dist, dist, div, old_dist, old_div))
        if dist < old_dist:
            best_dist[query] = (dist, div)
    return best_dist


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
    parser.add_argument('-v', '--verbose', action='store_true',
        help='Produce more output.')
    parser.add_argument('-F', '--force', action='store_true',
        help='Replace output file, if present.')
    args = parser.parse_args()

    if os.path.exists(args.output):
        if args.force:
            os.remove(args.output)
        else:
            sys.stderr.write(f'Skipping output file `{args.output}` (use --force to rewrite)\n')
            return

    targets = list(map(str.strip, args.genotype.split(',')))
    ploidy = len(targets)
    unique_targets = sorted(set(targets))
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(args.input, require_index=False) as f:
        pysam.set_verbosity(save)
        contigs, distances, aln_lens = _load_distances(f, unique_targets, args.verbose)
    best_dist = _process_distances(ploidy, contigs, distances, aln_lens, args.verbose)

    tmp_filename = common.temporary_filename(args.output)
    with common.open(tmp_filename, 'w') as out:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        out.write('genotype\tdistance\tdivergence\n')
        for query, (dist, div) in sorted(best_dist.items(), key=operator.itemgetter(1)):
            out.write('{}\t{:.0f}\t{:.10f}\n'.format(contigs.fmt_gt(query), dist, div))
    os.rename(tmp_filename, args.output)


if __name__ == '__main__':
    main()
