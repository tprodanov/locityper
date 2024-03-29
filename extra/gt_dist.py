#!/usr/bin/env python3

import argparse
from collections import defaultdict
import itertools
import sys
import operator
import os
import re
import numpy as np

import common


class Distances:
    def __init__(self, discarded_path, paf_path, verbose=False):
        self.paf_path = paf_path
        discarded = {}
        if os.path.exists(discarded_path):
            with common.open(discarded_path) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    hap, haps2 = line.strip().split('=')
                    discarded[hap.strip()] = tuple(map(str.strip, haps2.split(',')))
        elif verbose:
            sys.stderr.write(f'Cannot open `{discarded_path}`, assume there are no discarded haplotypes\n')

        def group(hap):
            return (hap,) + discarded.get(hap, ())

        seq_lengths = {}
        self.distances = defaultdict(dict)
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
                    self.distances[hap1a][hap2a] = dist
                    self.distances[hap2a][hap1a] = dist

        for hap, length in seq_lengths.items():
            for hap1, hap2 in itertools.product(group(hap), repeat=2):
                self.distances[hap1][hap2] = (0, length)
                self.distances[hap2][hap1] = (0, length)

        pattern = re.compile(r'[._][1-9]$')
        self.sample_haps = defaultdict(list)
        for hap in self.distances:
            m = re.search(pattern, hap)
            if m:
                self.sample_haps[hap[:m.start()]].append(hap)

    def get_sample_haplotypes(self, sample):
        return self.sample_haps.get(sample, ())

    def all_distances(self, genotype):
        hap_dists = []
        for hap in genotype:
            curr_dists = self.distances.get(hap)
            if curr_dists is None:
                return None
            hap_dists.append(curr_dists.items())

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
        return pred_dists

    def calc_distance(self, gt1, gt2):
        assert len(gt1) == len(gt2)
        best_div = np.inf
        best_distances = None

        for perm2 in itertools.permutations(gt2):
            distances = []
            sum_edit = 0
            sum_size = 0
            for hap1, hap2 in zip(gt1, perm2):
                if hap1 is None:
                    distances.append((None, None))
                    continue
                try:
                    edit, size = self.distances[hap1][hap2]
                except KeyError:
                    sys.stderr.write(f'Cannot calculate distance between {",".join(gt1)} and {",".join(gt2)}'
                        f' (missing distance {hap1} - {hap2}) (see {self.paf_path})\n')
                sum_edit += edit
                sum_size += size
                distances.append((edit, size))

            div = sum_edit / sum_size if sum_size else np.inf
            if div <= best_div:
                best_div = div
                best_distances = distances
        return GtDist(best_distances)

    def find_closest_loo(self, gt):
        loo_gt = []
        distances = []
        for hap in gt:
            if hap is None:
                loo_gt.append(None)
                distances.append((None, None))
                continue

            best_hap = None
            best_div = np.inf
            best_edit = None
            for hap2, (edit, size) in self.distances[hap].items():
                if edit / size < best_div and hap2 not in gt:
                    best_div = edit / size
                    best_edit = (edit, size)
                    best_hap = hap2
            loo_gt.append(best_hap)
            distances.append(best_edit)
        return loo_gt, GtDist(distances)

    def average_divergence(self):
        edits = []
        divs = []
        for hap1, hap_dists in self.distances.items():
            for hap2, (edit, size) in hap_dists.items():
                if hap1 < hap2:
                    edits.append(edit)
                    divs.append(edit / size)
        return np.mean(edits), np.mean(divs)


def edit_to_str(edit, size):
    div = edit / size
    qv = np.inf if div == 0 else -10 * np.log10(div)
    return f'{edit}\t{size}\t{div:.9f}\t{qv:.9f}'


class GtDist:
    def __init__(self, distances):
        self.distances = distances

    def iter_strs(self):
        sum_edit = 0
        sum_size = 0
        all_present = True
        for edit, size in self.distances:
            if edit is None:
                all_present = False
                yield 'NA\tNA\tNA\tNA'
            else:
                sum_edit += edit
                sum_size += size
                yield edit_to_str(edit, size)

        if all_present:
            yield edit_to_str(sum_edit, sum_size)
        else:
            yield 'NA\tNA\tNA\tNA'


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
        pred_dists = distances.all_distances(genotype)
        if pred_dists is None:
            out.write(f'{gt_str}\t*\tNA\tNA\tNA\tNA\n')
            continue
        pred_dists = sorted(pred_dists.items(), key=operator.itemgetter(1))
        for query, (div, edit, size) in pred_dists[:max_entries]:
            loo = 'T' if is_loo(genotype, query) else 'F'
            out.write(f'{gt_str}\t{",".join(query)}\t{loo}\t{edit}\t{size}\t{div:.9f}\n')


def main():
    parser = argparse.ArgumentParser(
        description='Calculating distances between a target genotype and all other genotypes',
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
    parser.add_argument('-v', '--verbose', action='store_true',
        help='Output more information to stderr.')
    args = parser.parse_args()

    genotypes = load_target_genotypes(args)
    distances = Distances(args.discarded, args.input)

    max_entries = args.max_entries or sys.maxsize
    with common.open(args.output, 'w') as out:
        out.write(f'# {" ".join(sys.argv)}\n')
        out.write('target\tquery\tloo\tedit_dist\taln_size\tdivergence\n')
        calc_gt_distances(genotypes, distances, out, max_entries)


if __name__ == '__main__':
    main()
