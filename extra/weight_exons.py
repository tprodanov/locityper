#!/usr/bin/env python3

import argparse
import sys
import os
from Bio import SeqIO

import common


def process_gtf(gtf_file, hap_name, hap_len, gene):
    """
    Returns list [(start, end, ty)] where ty = 0: CDSl; ty = 1: non-CDS part of the gene; ty = 2: rest.
    """
    # triples (pos, is_start? T/F, type: 0 exon, 1 gene, 2 boundaries).
    positions = [(0, True, 2), (hap_len, False, 2)]

    cds = []
    tag = f'gene_name "{gene}";'
    for line in gtf_file:
        if line.startswith('#'):
            continue
        line = line.split('\t')
        if tag not in line[8]:
            continue
        ty = line[2]

        start = int(line[3])
        end = int(line[4])
        assert end <= hap_len
        if start >= end:
            continue

        if ty == 'gene':
            positions.append((start, True, 1))
            positions.append((end, False, 1))
        elif ty == 'CDS' or ty == 'start_codon' or ty == 'stop_codon':
            positions.append((start, True, 0))
            positions.append((end, False, 0))

    enabled = [0, 0, 0]
    positions.sort()
    out = []
    last_pos = 0
    for pos, is_start, ty in positions:
        if last_pos < pos:
            curr = min(i for i, v in enumerate(enabled) if v)
            out.append((last_pos, pos, curr))
        enabled[ty] += 1 if is_start else -1
        assert enabled[ty] >= 0
        last_pos = pos
    assert not any(enabled) and last_pos == hap_len
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
    parser.add_argument('-w', '--weights', metavar='FLOAT', nargs=3, default='1.0 0.5 0.01',
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

    if args.discarded and os.path.exists(args.discarded):
        with common.open(args.discarded) as f:
            for line in f:
                left, right = map(str.strip, line.split('='))
                for name in map(str.strip, right.split(',')):
                    entries[name] = entries[left]
    elif args.discarded:
        sys.stderr.write(f'WARN: {args.discarded} does not exist\n')

    with common.open(args.output, 'w') as out:
        for name, curr_entries in entries.items():
            for start, end, ix in curr_entries:
                out.write(f'{name}\t{start}\t{end}\t{weights[ix]:.8g}\n')


if __name__ == '__main__':
    main()
