#!/usr/bin/env python3

import os
import sys
import argparse
import gzip
from tqdm import tqdm

import common
import summarize


def process_file(prefix, filename, out, thresholds):
    with gzip.open(filename, 'rt') as f:
        total_base = None
        total_call = None
        saved_lines = [None] * len(thresholds)

        for line in f:
            if 'total baseline variants' in line:
                total_base = int(line.split(':')[1].strip())
            elif 'total call variants' in line:
                total_call = int(line.split(':')[1].strip())
            elif not line.startswith('#'):
                line = list(map(float, line.strip().split()))
                score = line[0]
                for i, thresh in enumerate(thresholds):
                    if score >= thresh:
                        saved_lines[i] = line

        for thresh, line in zip(thresholds, saved_lines):
            out.write(f'{prefix}{total_base}\t{total_call}\t{thresh:g}\t')
            if line is not None:
                out.write('\t'.join(map('{:.10g}'.format, line)))
                if len(line) < 8:
                    out.write('\tNA' * (8 - len(line)))
                out.write('\n')
            else:
                # FN = #base.
                out.write(f'NA\t0\t0\t0\t{total_base}\tNA\tNA\tNA\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='STR', required=True,
        help='Path to vcfeval output directories. Tags within `{tag}` are automatically found.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    parser.add_argument('-t', '--thresholds', metavar='STR', default='0,10,20,30',
        help='Score thresholds via comma [%(default)s].')
    args = parser.parse_args()

    tags, tag_tuples = summarize.load_tags(args.input)
    thresholds = list(map(int, args.thresholds.split(',')))
    with common.open(args.output, 'w') as out:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        out.write(''.join(map('{}\t'.format, tags)))
        out.write('type\ttotal_base\ttotal_call\tthreshold\tscore\ttp_base\tfp\ttp_call\tfn\tprecision\trecall\tf1\n')
        for tup in tqdm(tag_tuples):
            dir = args.input.format(**dict(zip(tags, tup)))
            prefix = ''.join(map('{}\t'.format, tup))
            process_file(f'{prefix}snps\t', os.path.join(dir, 'snp_roc.tsv.gz'), out, thresholds)
            process_file(f'{prefix}indels\t', os.path.join(dir, 'non_snp_roc.tsv.gz'), out, thresholds)
            process_file(f'{prefix}all\t', os.path.join(dir, 'weighted_roc.tsv.gz'), out, thresholds)


if __name__ == '__main__':
    main()
