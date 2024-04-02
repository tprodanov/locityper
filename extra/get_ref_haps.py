#!/usr/bin/env python3

import argparse
import os
import gzip
from collections import defaultdict
import re
import sys

import common


def process_dir(dir, locus, out):
    disc_path = os.path.join(dir, 'discarded_haplotypes.txt')
    discarded = []
    if os.path.exists(disc_path):
        with common.open(disc_path) as f:
            for line in f:
                if not line.startswith('#'):
                    gt, gts2 = line.strip().split('=')
                    gt = gt.strip()
                    gts2 = list(map(str.strip, gts2.split(',')))
                    for gt2 in gts2:
                        discarded.append((gt, gt2))

    haplotypes = []
    with gzip.open(os.path.join(dir, 'haplotypes.fa.gz'), 'rt') as f:
        for line in f:
            if line.startswith('>'):
                haplotypes.append(line.split()[0][1:].strip())

    suffix = re.compile(r'[._][12]$')
    sample_to_haps = defaultdict(list)
    for hap in haplotypes:
        m = re.search(suffix, hap)
        if m:
            sample = hap[:m.start()]
            sample_to_haps[sample].append(hap)
    for hap0, hap in discarded:
        m = re.search(suffix, hap)
        if m:
            sample = hap[:m.start()]
            sample_to_haps[sample].append(hap0)

    for sample, haps in sample_to_haps.items():
        for hap in haps:
            out.write(f'{locus}\t{sample}\t{hap}\n')


def main():
    parser = argparse.ArgumentParser(
        description='Identify reference haplotypes for each sample and locus.')
    parser.add_argument('-d', '--database', metavar='DIR', required=True,
        help='Locityper database.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    args = parser.parse_args()

    indir = os.path.join(args.database, 'loci')
    with common.open(args.output, 'w') as out:
        out.write('locus\tsample\thaplotype\n')
        for locus in os.listdir(indir):
            subdir = os.path.join(indir, locus)
            if os.path.exists(os.path.join(subdir, 'success')):
                process_dir(subdir, locus, out)


if __name__ == '__main__':
    main()
