#!/usr/bin/env python3

import argparse
import sys
from Bio import SeqIO

import common


def process_gtf(gtf_file, hap_name, hap_len, gene):
    """
    Returns list [(start, end, ty)] where ty = 0: CDSl; ty = 1: non-CDS part of the gene; ty = 2: rest.
    """
    gene_start = None
    gene_end = None
    cds = []
    tag = f'gene_name "{gene}";'
    for line in gtf_file:
        if line.startswith('#'):
            continue
        line = line.split('\t')
        if tag not in line[8]:
            continue
        ty = line[2]
        coords = (int(line[3]), int(line[4]))
        assert coords[1] <= hap_len
        if ty == 'gene':
            assert gene_start is None
            gene_start, gene_end = coords
        elif ty == 'CDS':
            cds.append(coords)

    if gene_start is None:
        assert not cds
        sys.stderr.write(f'WARN: Gene {gene} not found on haplotype {hap_name}\n')
        return [(0, hap_len, 2)]

    out = []
    if gene_start > 0:
        out.append((0, gene_start, 2))
    last_end = gene_start
    cds.sort()
    for start, end in cds:
        if last_end < start:
            out.append((last_end, start, 1))
        start = max(start, last_end)
        end = min(end, gene_end)
        if start < end:
            out.append((start, end, 0))
        last_end = max(last_end, end)
    if last_end < gene_end:
        out.append((last_end, gene_end, 1))
    if gene_end < hap_len:
        out.append((gene_end, hap_len, 2))
    return out


def main():
    parser = argparse.ArgumentParser(description='Define subregion weights based on exons.')
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Fasta file with input haplotypes.')
    parser.add_argument('-a', '--annot', metavar='STR', required=True,
        help='Path to Immuannot GTF files with gene annotation. '
            'Haplotype in the filename should be replaced with `{}`.')
    parser.add_argument('-g', '--gene', metavar='STR', required=True,
        help='Gene name.')
    parser.add_argument('-w', '--weights', metavar='FLOAT', nargs=3, default='1.0 0.5 0.1',
        help='Weights for exons, introns and intergenic sequence [%(default)s].')
    parser.add_argument('-d', '--discarded', metavar='FILE',
        help='Optional: Text file denoting discarded haplotypes.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output BED file.')
    args = parser.parse_args()

    if isinstance(args.weights, str):
        weights = tuple(map(float, args.weights.split()))
    else:
        weights = tuple(map(float, args.weights))

    haps = []
    with common.open(args.input) as f:
        for record in SeqIO.parse(f, 'fasta'):
            haps.append((record.id, len(record.seq)))

    entries = {}
    for name, length in haps:
        with common.open(args.annot.replace('{}', name)) as gtf:
            entries[name] = process_gtf(gtf, name, length, args.gene)

    if args.discarded:
        with common.open(args.discarded) as f:
            for line in f:
                left, right = map(str.strip, line.split('='))
                for name in map(str.strip, right.split(',')):
                    entries[name] = entries[left]

    with common.open(args.output, 'w') as out:
        for name, curr_entries in entries.items():
            for start, end, ix in curr_entries:
                out.write(f'{name}\t{start}\t{end}\t{weights[ix]:.8g}\n')


if __name__ == '__main__':
    main()
