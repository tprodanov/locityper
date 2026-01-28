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

    def __len__(self):
        return self.end - self.start


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


def find_redundant(loci, out_csv):
    loci.sort()
    current = []
    redundant = []

    for locus in loci:
        # reversed because this way we don't care about loci, removed in the process.
        for i in reversed(range(len(current))):
            oth = current[i]
            if oth.chrom != locus.chrom or oth.end <= locus.start:
                current.pop(i)

        locus_str = str(locus)
        this_redundant = False
        for oth in current:
            assert oth.chrom == locus.chrom
            overlap = min(locus.end, oth.end) - max(locus.start, oth.start)
            assert overlap > 0
            if overlap < len(locus):
                sys.stderr.write(f'...    {locus_str:40}    overlaps     {oth}\n')
            elif locus.start == oth.start and locus.end == oth.end:
                sys.stderr.write(f'===    {locus_str:40} identical with  {oth}\n')
                this_redundant = True
            else:
                sys.stderr.write(f'!!!    {locus_str:40}  contained in   {oth}\n')
                this_redundant = True
            if out_csv:
                out_csv.write(f'{oth.name}\t{locus.name}\t{overlap}'
                    f'\t{overlap / len(oth):.6f}\t{overlap / len(locus):.6f}\n')

        if this_redundant:
            redundant.append(locus)
        current.append(locus)
    return redundant


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
    parser.add_argument('-O', '--out-csv', metavar='FILE',
        help='Optional: output CSV file with overlap information.')
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

    if args.out_csv:
        out_csv = common.open(args.out_csv, 'w')
        out_csv.write('locus1\tlocus2\toverlap\tfrac_of1\tfrac_of2\n')
    else:
        out_csv = None

    redundant = find_redundant(loci, out_csv)
    for locus in redundant:
        if args.move:
            new_dir = os.path.join(args.output, locus.name)
            while os.path.isdir(new_dir):
                rand_dir = os.path.join(args.output, '{}-{:X}'.format(locus.name, random.randint(0, 0xffffffff)))
                sys.stderr.write(f'WARN: Output directory {new_dir} already exists, moving to {rand_dir}\n')
                new_dir = rand_dir
            os.rename(os.path.join(loci_dir, locus.name), new_dir)
    redundant_str = ', '.join(locus.name for locus in redundant[:5]) + (' ...' if len(redundant) > 5 else '')
    sys.stderr.write(f'{len(redundant)} / {len(loci)} redundant loci: {redundant_str}\n')


if __name__ == '__main__':
    main()
