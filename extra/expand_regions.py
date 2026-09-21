#!/usr/bin/env python3

import pysam
import sys
from common import open


def extend_region(vcf, chrom: str, start: int, end: int, name: str | None, args):
    new_start = start
    new_end = end
    for var in vcf.fetch(chrom, start, end):
        if max(var.info['AF']) < args.af or all(len(allele) < args.length for allele in var.alleles):
            continue
        new_start = min(new_start, max(0, var.start - args.padding))
        new_end = max(new_end, var.start + len(var.ref) + args.padding)
    left_extension = start - new_start
    right_extension = new_end - end

    name_str = f' ({name})' if name else ''
    if left_extension > args.max_extension or right_extension > args.max_extension:
        sys.stderr.write(
            f'! WARN  Could not extend region {chrom}:{start+1}-{end}{name_str}: '
            f'minimal extension {left_extension:,} bp (left) and {right_extension:,} bp (right)\n')
        return None
    sys.stderr.write(f'  INFO  Extending region {chrom}:{start+1}-{end}{name_str} by '
        f'{left_extension:,} bp (left) and {right_extension:,} bp (right)\n')
    return new_start, new_end


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Expand regions such that the boundaries do not overlap variants from the VCF file.')
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Input BED file.')
    parser.add_argument('-v', '--vcf', metavar='FILE', required=True,
        help='Input VCF file.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output BED file.')
    parser.add_argument('--af', type=float, metavar='NUM', default=0.005,
        help='Skip variants with allele fraction under this value [%(default)s].')
    parser.add_argument('-l', '--length', type=int, metavar='INT', default=2,
        help='Skip variants with length under this value [%(default)s].')
    parser.add_argument('-p', '--padding', type=int, metavar='INT', default=5,
        help='Add padding after extension [%(default)s]. Padded region is not checked for variant overlaps.')
    parser.add_argument('-m', '--max-extension', type=int, metavar='INT', default=50000,
        help='Do not allow region extension beyond this size [%(default)s].')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    with open(args.input) as f, open(args.output, 'w') as out:
        for line in f:
            if line.startswith('#'):
                continue
            split = line.strip().split('\t')
            chrom = split[0]
            start = int(split[1])
            end = int(split[2])
            name = split[3] if len(split) >= 4 else None
            new_region = extend_region(vcf, chrom, start, end, name, args)
            if new_region is None:
                continue
            new_start, new_end = new_region
            out.write(f'{chrom}\t{new_start}\t{new_end}')
            if name is not None:
                out.write('\t' + '\t'.join(split[3:]))
            out.write('\n')


if __name__ == '__main__':
    main()
