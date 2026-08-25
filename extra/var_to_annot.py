#!/usr/bin/env python3

import os.path
from collections import defaultdict
import common
import pysam


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='FILE',
        help='Input VCF file.')
    parser.add_argument('-R', '--ref-hap', metavar='STR',
        help='Reference haplotype name.')
    parser.add_argument('-r', '--regions', metavar='STR', nargs='+', required=True,
        help='Fetch regions. Format: "chr:start-end" or "chr:pos@radius".')
    parser.add_argument('-d', '--min-diff', metavar='INT', type=int, default=1,
        help='Minimum difference between shortest and longest alleles [%(default)s].')
    parser.add_argument('-o', '--output', metavar='FILE',
        help='Output CSV file with per-variant information.')
    parser.add_argument('-a', '--annotation', metavar='FILE',
        help='Output annotation CSV file, compatible with `annotate.py`. '
            'If starts with "+", appends instead of rewriting.')
    parser.add_argument('-n', '--name', metavar='STR',
        help='Locus name [default: last directory name of the input VCF file].')
    # parser.add_argument('-s', '--sum', action='store_true',
    #     help='For annotation, sum allele lengths of all variants in the region.')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.input)
    out_csv = common.open(args.output, 'wt') if args.output is not None else common.Sink()
    out_csv.write('chrom\tpos\thap\tallele_diff\n')

    if args.annotation is not None and args.annotation.startswith('+'):
        annot_mode = 'at'
        args.annotation = args.annotation.lstrip('+')
        write_header = not os.path.exists(args.annotation)
    else:
        annot_mode = 'wt'
        write_header = True
    out_annot = common.open(args.annotation, annot_mode) if args.annotation is not None else common.Sink()
    if write_header:
        out_annot.write('locus\thaplotype\tannotation\ttag\n')

    locus = args.name if args.name is not None else os.path.basename(os.path.dirname(args.input))

    NAN = float('nan')
    sum_diff = defaultdict(int)
    n_vars = 0
    for reg in args.regions:
        chrom, coord = reg.split(':')
        if '@' in coord:
            pos, radius = map(int, coord.split('@'))
            start = pos - radius
            end = pos + radius - 1
        else:
            start, end = map(int, coord.split('-'))
            start -= 1

        for var in vcf.fetch(chrom, start, end):
            allele_lens = list(map(len, var.alleles))
            if max(allele_lens) - min(allele_lens) < args.min_diff:
                continue
            n_vars += 1
            ref_len = allele_lens[0]
            if args.ref_hap is not None:
                out_annot.write(f'{locus}\t{args.ref_hap}\t0\t{var.pos}\n')
            for sample, data in var.samples.items():
                gt = data['GT']
                for i, allele_ix in enumerate(gt, 1):
                    hap = sample if len(gt) == 1 else f'{sample}.{i}'
                    curr_diff = allele_lens[allele_ix] - ref_len if allele_ix is not None else NAN
                    sum_diff[hap] += curr_diff
                    out_csv.write(f'{var.chrom}\t{var.pos}\t{hap}\t{curr_diff}\n')
                    out_annot.write(f'{locus}\t{hap}\t{curr_diff}\t{var.pos}\n')

    if n_vars > 1:
        if args.ref_hap is not None:
            out_annot.write(f'{locus}\t{args.ref_hap}\t0\tsum\n')
        for hap, annot in sorted(sum_diff.items()):
            out_annot.write(f'{locus}\t{hap}\t{annot}\tsum\n')


if __name__ == '__main__':
    main()
