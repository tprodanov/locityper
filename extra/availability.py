#!/usr/bin/env python3

import sys

import common
from gt_dist import Distances


def load_target_genotypes(f, distances):
    target_genotypes = []
    for line in f:
        sample = line.strip()
        haps = distances.get_sample_haplotypes(sample)
        if len(haps) == 0:
            continue
        elif len(haps) == 1:
            haps.append('*')
        elif len(haps) > 2:
            sys.stderr.write(f'WARN: Sample {sample} should have <= 2 haplotypes (found {len(haps)})\n')
            continue
        target_genotypes.append(haps)
    return target_genotypes


def load_haplotypes(f, distances):
    haplotypes = set()
    for line in f:
        line = line.strip()
        # Haplotype by itself
        if line in distances.distances:
            haplotypes.add(line)
        for hap in distances.get_sample_haplotypes(line):
            haplotypes.add(hap)
    return haplotypes


def find_closest(genotypes, distances, excl_haplotypes, out):
    for genotype in genotypes:
        closest_haps, gt_dists = distances.find_closest(genotype, loo=True, excl_haps=excl_haplotypes)
        closest_haps = ['*' if hap is None else hap for hap in closest_haps]
        closest_haps.append(','.join(closest_haps))

        this_haps = (genotype[0], genotype[1], ','.join(genotype))
        for this_hap, oth_hap, s in zip(this_haps, closest_haps, gt_dists.iter_strs()):
            if this_hap != '*':
                out.write(f'{this_hap}\t{oth_hap}\t{s}\n')


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='For each sample in --S1, calculate haplotype availability assuming leave-one-out scheme. \
        In addition, availability can be calculated assuming all samples from --S2 are discarded.')
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Input PAF[.gz] file with pairwise distances.')
    parser.add_argument('-d', '--discarded', metavar='FILE', required=False,
        help='File with discarded haplotypes.')
    parser.add_argument('-S', '--samples', metavar='FILE', required=True,
        help='Calculate availability for each sample from this file.')
    parser.add_argument('--include', metavar='FILE', required=False,
        help='Assume that only these samples/haplotypes are available.')
    parser.add_argument('--exclude', metavar='FILE', required=False,
        help='Assume that all samples/haplotypes except these are available.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    args = parser.parse_args()

    distances = Distances(args.discarded, args.input)
    with common.open(args.samples) as f:
        target_genotypes = load_target_genotypes(f, distances)

    assert args.include is None or args.exclude is None, '--include and --exclude cannot go together'

    if args.include:
        with common.open(args.include) as f:
            incl_haplotypes = load_haplotypes(f, distances)
        excl_haplotypes = set(distances.distances.keys()) - incl_haplotypes
    elif args.exclude:
        with common.open(args.exclude) as f:
            excl_haplotypes = load_haplotypes(f, distances)
    else:
        excl_haplotypes = ()

    with common.open(args.output, 'w') as out:
        out.write('target\tclosest\tedit\tsize\tdiv\tqv\n')
        find_closest(target_genotypes, distances, excl_haplotypes, out)


if __name__ == '__main__':
    main()
