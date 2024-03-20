#!/usr/bin/env python3

import argparse
from collections import defaultdict
import itertools
import sys
import operator

import common


def load_distances(discarded_path, paf_path):
    discarded = {}
    with common.open(discarded_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            hap, haps2 = line.strip().split('=')
            discarded[hap.strip()] = tuple(map(str.strip, haps2.split(',')))

    def group(hap):
        return (hap,) + discarded.get(hap, ())

    seq_lengths = {}
    distances = defaultdict(dict)
    with common.open(paf_path) as paf:
        for line in paf:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            hap1 = line[0]
            hap2 = line[5]

            if hap1 not in seq_lengths:
                len1 = int(line[1])
                for hap1a in group(hap1):
                    seq_lengths[hap1a] = len1
            if hap2 not in seq_lengths:
                len2 = int(line[6])
                for hap2a in group(hap2):
                    seq_lengths[hap2a] = len2
            nmatches = int(line[10])
            aln_size = int(line[11])
            assert aln_size != 0, f'Missing alignment between {hap1} and {hap2}'
            assert nmatches <= aln_size
            dist = (aln_size - nmatches, aln_size)
            for hap1a, hap2a in itertools.product(group(hap1), group(hap2)):
                distances[hap1a][hap2a] = dist
                distances[hap2a][hap1a] = dist

    for hap, length in seq_lengths.items():
        distances[hap][hap] = (0, length)
    return { hap: list(dists.items()) for hap, dists in distances.items() }


def get_genotype(s, split, sep):
    tup = tuple(map(str.strip, s.strip().split(split)))
    assert tup
    if len(tup) == 1:
        sample = tup[0]
        return (sample, (f'{sample}{sep}1', f'{sample}{sep}2'))
    gt_str = ','.join(tup)
    assert len(tup) == 2, f'Diploid genotype required ({gt_str})'
    return gt_str, tup


def load_target_genotypes(args):
    genotypes = []
    if args.genotype:
        for val in args.genotype:
            genotypes.append(get_genotype(val, ',', args.sep))
    if args.gt_file:
        with common.open(args.gt_file) as f:
            for line in f:
                genotypes.append(get_genotype(line, None, args.sep))
    return genotypes


def is_loo(target, query):
    return all(h1 != h2 for h1, h2 in itertools.product(target, query))


def calc_gt_distances(genotypes, distances, out, max_entries):
    for gt_str, genotype in genotypes:
        hap_dists = []
        for hap in genotype:
            curr_dists = distances.get(hap)
            assert curr_dists, f'Unknown haplotype {hap}'
            hap_dists.append(curr_dists)

        pred_dists = {}
        # for dist_combin in itertools.product(*hap_dists):
        for (hap1, (edit1, size1)), (hap2, (edit2, size2)) in itertools.product(*hap_dists):
            edit = edit1 + edit2
            size = size1 + size2
            div = edit / size
            query = (hap1, hap2) if hap1 <= hap2 else (hap2, hap1)
            last_pred = pred_dists.get(query)
            if last_pred is None or last_pred[0] > div:
                pred_dists[query] = (div, edit, size)

        pred_dists = sorted(pred_dists.items(), key=operator.itemgetter(1))
        for query, (div, edit, size) in pred_dists[:max_entries]:
            loo = 'T' if is_loo(genotype, query) else 'F'
            out.write(f'{gt_str}\t{",".join(query)}\t{loo}\t{edit}\t{size}\t{div:.9f}\n')


def main():
    parser = argparse.ArgumentParser(
        description='Calculates distance between a target genotype and all other genotypes',
        usage='%(prog)s -i alns.paf -d discarded.txt (-g hap1,hap2 | -G genotypes.txt) -o out.csv')
    parser.add_argument('-i', '--input', metavar='FILE',
        help='Input PAF[.gz] file with pairwise distances.')
    parser.add_argument('-d', '--discarded', metavar='FILE',
        help='File with discarded haplotypes.')
    parser.add_argument('-g', '--genotype', metavar='STR', nargs='+',
        help='One or more target genotype (haplotypes through comma) or sample name.')
    parser.add_argument('-G', '--gt-file', metavar='FILE',
        help='List of target genotypes/samples.')
    parser.add_argument('-o', '--output', metavar='FILE',
        help='Output CSV file.')
    parser.add_argument('-s', '--sep', metavar='STR', default='.',
        help='Separator between sample and haplotype [default: %(default)s].')
    parser.add_argument('-n', '--max-entries', metavar='INT',
        help='Output at most INT entries per target genotype [default: all].')
    args = parser.parse_args()

    genotypes = load_target_genotypes(args)
    distances = load_distances(args.discarded, args.input)

    max_entries = args.max_entries or sys.maxsize
    with common.open(args.output, 'w') as out:
        out.write('target\tquery\tloo\tedit_dist\taln_size\tdivergence\n')
        calc_gt_distances(genotypes, distances, out, max_entries)


if __name__ == '__main__':
    main()
