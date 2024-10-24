#!/usr/bin/env python3

import argparse
import collections
import sys

import common


def main():
    parser = argparse.ArgumentParser(description='Annotate genotyping results.',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Locityper CSV file for multiple loci and samples.')
    parser.add_argument('-a', '--annotation', metavar='FILE', required=True,
        help='Tab-separated annotation for locus haplotypes.\n'
        'Columns: locus, haplotypes, annotation, tag. Tag is an optional column,\n'
        'used if the same locus has multiple annotations.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    args = parser.parse_args()

    with common.open(args.annotation) as f:
        annotation = collections.defaultdict(lambda: collections.defaultdict(dict))
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            locus = line[0]
            hap = line[1]
            annot = line[2]
            assert annot != ''
            tag = line[3] if len(line) >= 4 else '*'
            if tag == '':
                tag = '*'
            if hap in annotation[locus][tag]:
                sys.stderr.write(f'WARN: Haplotype {hap} appears twice for locus {locus}, tag {tag}\n')
            annotation[locus][tag][hap] = annot

    with common.open(args.input) as inp, common.open(args.output, 'w') as out:
        out.write(f'# {" ".join(sys.argv)}\n')
        out.write('sample\tlocus\ttag\thap1\thap2\n')
        for row in common.read_csv(inp):
            locus = row['locus']
            locus_annot = annotation.get(locus)
            if not locus_annot:
                continue

            sample = row['sample']
            gt = row['genotype']
            if gt == '*':
                for tag in locus_annot:
                    out.write(f'{sample}\t{locus}\t{tag}\t<NOCALL>\t<NOCALL>\n')
                continue

            for tag, tag_annot in locus_annot.items():
                out.write(f'{sample}\t{locus}\t{tag}')
                for hap in gt.split(','):
                    out.write('\t{}'.format(tag_annot.get(hap, '<UNKNOWN>')))
                out.write('\n')


if __name__ == '__main__':
    main()
