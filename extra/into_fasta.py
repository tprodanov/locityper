#!/usr/bin/env python3

import numpy as np
import argparse
import sys
import os
from collections import defaultdict
import multiprocessing
from Bio import SeqIO
import gzip

import common
import into_vcf


def parse_float(val):
    return float(val) if val != 'NA' else np.nan


def load_haplotypes(path, subset_loci):
    loci_dir = os.path.join(path, 'loci')
    if not os.path.isdir(loci_dir):
        common.error(f'Directory `{path}` does not contain "loci" subdirectory')
    if subset_loci is not None:
        subset_loci = set(subset_loci)

    haplotypes = {}
    for entry in os.scandir(os.path.join(path, 'loci')):
        if not entry.is_dir():
            continue
        locus = entry.name
        if subset_loci is not None and locus not in subset_loci:
            continue
        locus_dir = entry.path
        reader = SeqIO.parse(gzip.open(os.path.join(locus_dir, 'haplotypes.fa.gz'), 'rt'), 'fasta')
        haplotypes[locus] = { rec.id: str(rec.seq) for rec in reader }

    sys.stderr.write(f'Loaded {len(haplotypes)} loci\n')
    return haplotypes


def process_locus(locus, samples, preds, haplotypes, out_dir, compress):
    extension = 'fasta.gz' if compress else 'fasta'
    with common.open(f'{out_dir}/{locus}.{extension}', 'w') as out:
        for sample in samples:
            gt, features = preds[sample]
            for i, allele in enumerate(gt, 1):
                try:
                    hap = haplotypes[allele]
                except KeyError:
                    sys.stderr.write(f'Could not find haplotype {allele} for locus {locus} and sample {sample}\n')
                    continue

                out.write(f'>{sample}.{i} {allele}')
                for key, val in features.items():
                    if key != 'GQ':
                        out.write(f' {key}={val}')
                out.write('\n')
                for i in range(0, len(hap), 120):
                    out.write(hap[i:i+120] + '\n')
    return locus


def main():
    parser = argparse.ArgumentParser(description='Convert Locityper predictions into FASTA file.',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='CSV file with Locityper predictions.')
    parser.add_argument('-d', '--database', metavar='FILE', required=True,
        help='Path to Locityper database.')
    parser.add_argument('-o', '--output', metavar='DIR', required=True,
        help='Output directory.')
    parser.add_argument('-@', '--threads', metavar='INT', type=int, default=8,
        help='Analyze loci in this many threads [%(default)s].')
    parser.add_argument('--subset-loci', metavar='STR', nargs='+',
        help='Limit the analysis to these loci.')
    parser.add_argument('--compress', choices=('yes', 'y', 'no', 'n'), default='yes',
        help='Compress output files [%(default)s].')
    args = parser.parse_args()

    common.mkdir(args.output)
    sys.stderr.write('Loading haplotypes\n')
    haplotypes = load_haplotypes(args.database, args.subset_loci)
    sys.stderr.write('Loading genotype predictions\n')
    with common.open(args.input) as f:
        preds, samples = into_vcf.load_predictions(f, '1')
    sys.stderr.write(f'Loaded {len(samples)} samples\n')
    loci_inters = set(haplotypes.keys()) & set(preds.keys())
    if len(loci_inters) < len(haplotypes):
        sys.stderr.write('WARN: Genotype predictions missing for {} loci, such as {}\n'.format(
            len(haplotypes) - len(loci_inters), list(set(haplotypes.keys()) - loci_inters)[:5].join(', ')))

    total = len(loci_inters)
    finished = 0
    def callback(locus):
        nonlocal finished
        finished += 1
        sys.stderr.write(f'Finished [{finished:3}/{total}] {locus}\n')

    compress = args.compress.startswith('y')
    with multiprocessing.Pool(args.threads) as pool:
        results = [
            pool.apply_async(process_locus,
                (locus, samples, preds[locus], haplotypes[locus], args.output, compress),
                callback=callback)
            for locus in loci_inters]
        for res in results:
            res.get()
        pool.close()
        pool.join()


if __name__ == '__main__':
    main()
