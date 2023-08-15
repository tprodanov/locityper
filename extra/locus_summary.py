#!/usr/bin/env python3

import argparse
import os
import gzip
import numpy as np
import tqdm
import sys

from common import open_stream


def get_divergencies(filename):
    with gzip.open(filename, 'rt') as inp:
        haplotypes = set()
        divergencies = []
        for line in inp:
            line = line.strip().split('\t')
            haplotypes.add(line[0])
            haplotypes.add(line[5])
            dv = next(tag for tag in line[12:] if tag.startswith('dv:'))[5:]
            divergencies.append(dv)
    return haplotypes, np.array(divergencies, dtype=np.float64)


def process_locus(dir, out):
    with open(os.path.join(dir, 'ref.bed')) as inp:
        chrom, start, end, name = next(inp).strip().split('\t')
        start = int(start)
        end = int(end)

    with gzip.open(os.path.join(dir, 'kmers.gz'), 'rt') as inp:
        next(inp)
        kmers = np.fromiter(map(int, next(inp).strip().split()), dtype=np.uint16)

    with gzip.open(os.path.join(dir, 'haplotypes.fa.gz'), 'rt') as inp:
        m = sum(line.startswith('>') for line in inp)

    haplotypes, divergencies = get_divergencies(os.path.join(dir, 'all_haplotypes.paf.gz'))
    n = len(haplotypes)
    assert n * (n - 1) / 2 == len(divergencies)
    aver_dist = np.mean(divergencies)
    out.write(f'{name}\t{chrom}:{start+1}-{end}\t{end-start}\t{n}\t{m}\t{aver_dist:.7f}')
    for i in range(1, N_CDFS + 1):
        out.write('\t{:.4f}'.format(np.mean(kmers <= i)))
    out.write('\n')
    out.flush()


N_CDFS = 5


def main():
    parser = argparse.ArgumentParser(description='Calculates average distance between gene alleles')
    parser.add_argument('input', metavar='DIR',
        help='Input directory, with multiple loci.')
    parser.add_argument('-L', '--loci', metavar='FILE',
        help='Only analyze loci with names from this list.')
    parser.add_argument('-o', '--output', metavar='FILE', required=False,
        help='Output CSV file [stdout].')
    args = parser.parse_args()

    out = open_stream(args.output, 'w')
    out.write('locus\tregion\tlength\thaplotypes\tfilt_haplotypes\taver_dist')
    for i in range(1, N_CDFS + 1):
        out.write(f'\tcdf{i}')
    out.write('\n')

    if args.loci is not None:
        with open(args.loci) as inp:
            loci = list(filter(bool, map(str.strip, inp)))
    else:
        loci = os.listdir(args.input)

    for locus in tqdm.tqdm(loci):
        dir = os.path.join(args.input, locus)
        if os.path.isfile(os.path.join(dir, 'success')):
            process_locus(dir, out)


if __name__ == '__main__':
    main()
