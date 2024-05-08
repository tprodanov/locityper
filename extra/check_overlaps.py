#!/usr/bin/env python3

import os
import sys
import collections
import argparse
import random

import common


class Locus(collections.namedtuple('Locus', 'chrom start end name')):
    def __str__(self):
        return f'{self.name} ({self.chrom}:{self.start+1:,}-{self.end:,})'

    def __lt__(self, oth):
        # Order first by chromosome, then by start, and then in descending order, by end.
        return (self.chrom, self.start, -self.end).__lt__((oth.chrom, oth.start, -oth.end))


def load_loci(indir):
    loci = []
    for subdir_name in os.listdir(indir):
        subdir = os.path.join(indir, subdir_name)
        if os.path.exists(os.path.join(subdir, 'success')):
            with open(os.path.join(subdir, 'ref.bed')) as inp:
                chrom, start, end, locus_name = next(inp).strip().split()
                start = int(start)
                end = int(end)

                if locus_name != subdir_name:
                    sys.stderr.write(f'WARN: Locus name does not match directory name at {subdir} ({locus_name})\n')
                loci.append(Locus(chrom, start, end, subdir_name))
    return loci


def find_redundant(loci):
    loci.sort()
    furthest = None
    for locus in loci:
        if furthest is None or furthest.chrom != locus.chrom or furthest.end < locus.end:
            if furthest is not None and furthest.chrom == locus.chrom and furthest.end > locus.start:
                sys.stderr.write(f'...    {str(furthest):50} overlaps {locus}\n')
            furthest = locus
        else:
            sys.stderr.write(f'!!!    {str(locus):50} contained completely in {furthest}\n')
            yield locus


def main():
    parser = argparse.ArgumentParser(
        description='Check for overlaps across target loci.',
        usage='%(prog)s db [--move]')
    parser.add_argument('input', metavar='DIR',
        help='Locityper loci database.')
    parser.add_argument('-m', '--move', action='store_true', dest='move',
        help='Move redundant loci to `--output`.')
    parser.add_argument('-o', '--output', metavar='DIR',
        help='Move discarded loci to this directory (default: <input>/redundant).')
    args = parser.parse_args()

    loci_dir = os.path.join(args.input, 'loci')
    sys.stderr.write(f'Loading coordinates from {loci_dir}/*/ref.bed\n')
    loci = load_loci(loci_dir)
    sys.stderr.write(f'Found {len(loci)} loci\n')
    if args.move:
        if args.output is None:
            args.output = os.path.join(args.input, 'redundant')
        sys.stderr.write(f'Writing redundant loci to {args.output}\n')
        common.mkdir(args.output)

    redundant = 0
    for locus in find_redundant(loci):
        redundant += 1
        if args.move:
            new_dir = os.path.join(args.output, locus.name)
            while os.path.isdir(new_dir):
                rand_dir = os.path.join(args.output, '{}-{:X}'.format(locus.name, random.randint(0, 0xffffffff)))
                sys.stderr.write(f'WARN: Output directory {new_dir} already exists, moving to {rand_dir}\n')
                new_dir = rand_dir
            os.rename(os.path.join(loci_dir, locus.name), new_dir)
    sys.stderr.write(f'{redundant} / {len(loci)} redundant loci\n')


if __name__ == '__main__':
    main()
