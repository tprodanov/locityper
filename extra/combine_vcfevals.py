#!/usr/bin/env python3

import os
import sys
import argparse
import gzip
import pysam
import numpy as np
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


def process_vcfs(prefix, base_vcf, calls_vcf, out, thresholds):
    assert len(base_vcf.samples) == len(calls_vcf.samples) == 1
    total_baseline = np.zeros(3, dtype=int)
    for rec in base_vcf:
        assert len(rec.alleles) == 2
        if rec.allele_variant_types[1] == 'SNP':
            total_baseline[[0, 2]] += 1
        else:
            total_baseline[[1, 2]] += 1

    over_thresh = np.zeros((3, len(thresholds)), dtype=int)
    for rec in calls_vcf:
        assert len(rec.alleles) == 2
        qual = rec.samples[0]['GQ']
        if rec.allele_variant_types[1] == 'SNP':
            over_thresh[[0, 2], :] += qual >= thresholds
        else:
            over_thresh[[1, 2], :] += qual >= thresholds

    for i, ty in enumerate(('snps', 'indels', 'all')):
        curr_base = total_baseline[i]
        for thresh, curr_calls in zip(thresholds, over_thresh[i]):
            out.write(f'{prefix}{ty}\t{curr_base}\t{curr_calls}\t{thresh}\t')
            assert curr_base == 0 or curr_calls == 0
            if curr_base:
                out.write(f'NA\t0\t0\t0\t{curr_base}\tNA\tNA\tNA\n')
            else:
                out.write(f'NA\t0\t{curr_call}\t0\t0\tNA\tNA\tNA\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='STR', required=True,
        help='Path to vcfeval output directories. Tags within `{tag}` are automatically found.')
    parser.add_argument('-b', '--baseline', metavar='STR', required=True,
        help='Path to baseline calls, containing the same tags as --input.')
    parser.add_argument('-c', '--calls', metavar='STR', required=True,
        help='Path to variant calls, containing the same tags as --input.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    parser.add_argument('-t', '--thresholds', metavar='STR', default='0,10,20,30',
        help='Score thresholds via comma [%(default)s].')
    args = parser.parse_args()

    tags1, tag_tuples1 = summarize.load_tags(args.input)
    tags2, tag_tuples2 = summarize.load_tags(args.baseline)
    assert set(tags1) == set(tags2)

    tags = tags1
    tag_tuples = set(tag_tuples1) | set(tag_tuples2)

    thresholds = np.array(list(map(int, args.thresholds.split(','))))
    with common.open(args.output, 'w') as out:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        out.write(''.join(map('{}\t'.format, tags)))
        out.write('type\ttotal_base\ttotal_call\tthreshold\tscore\ttp_base\tfp\ttp_call\tfn\tprecision\trecall\tf1\n')
        for tup in tqdm(tag_tuples):
            fmt = dict(zip(tags, tup))
            dir = args.input.format(**fmt)
            if os.path.isdir(dir):
                prefix = ''.join(map('{}\t'.format, tup))
                process_file(f'{prefix}snps\t', os.path.join(dir, 'snp_roc.tsv.gz'), out, thresholds)
                process_file(f'{prefix}indels\t', os.path.join(dir, 'non_snp_roc.tsv.gz'), out, thresholds)
                process_file(f'{prefix}all\t', os.path.join(dir, 'weighted_roc.tsv.gz'), out, thresholds)
            else:
                with pysam.VariantFile(args.baseline.format(**fmt)) as base_vcf, \
                    pysam.VariantFile(args.calls.format(**fmt)) as calls_vcf:
                    process_vcfs(prefix, base_vcf, calls_vcf, out, thresholds)



if __name__ == '__main__':
    main()
