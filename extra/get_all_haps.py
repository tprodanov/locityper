#!/usr/bin/env python3

import os
import shutil
import gzip
from Bio import SeqIO

from common import open
from into_fasta import write_fasta


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='DIR',
        help='Database directory for a specific locus')
    parser.add_argument('-o', '--output', metavar='FILE',
        help='Output fasta[.gz] file with all haplotypes [default: DIR/all_haplotypes.fa.gz].')
    args = parser.parse_args()

    if args.output is None:
        args.output = os.path.join(args.input, 'all_haplotypes.fa.gz')

    in_fa = gzip.open(os.path.join(args.input, 'haplotypes.fa.gz'), 'rt')
    out_fa = open(args.output, 'w')
    discarded_filename = os.path.join(args.input, 'discarded_haplotypes.txt')
    if not os.path.exists(discarded_filename):
        shutil.copyfileobj(in_fa, out_fa)
        exit(0)

    reader = SeqIO.parse(in_fa, 'fasta')
    records = {}
    for rec in reader:
        seq = str(rec.seq)
        records[rec.id] = seq
        write_fasta(out_fa, rec.id, seq)

    with open(discarded_filename) as f:
        for line in f:
            lhs, rhs = line.split('=')
            seq = records[lhs.strip()]
            for name in rhs.split(','):
                write_fasta(out_fa, name.strip(), seq)


if __name__ == '__main__':
    main()
