#!/usr/bin/env python3

import numpy as np
import argparse
import sys
import pysam
import json
import gzip
import os
import collections


def create_vcf_header(chrom, reference, samples):
    header = pysam.VariantHeader()
    header.add_line('##contig=<ID={},length={}>'.format(chrom, reference.get_reference_length(chrom)))
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">')
    header.add_line('##FORMAT=<ID=GQ0,Number=1,Type=Float,'
        'Description="Initial genotype quality, not accounting for warnings">')
    header.add_line('##FORMAT=<ID=WARN,Number=1,Type=String,Description="Genotype warnings">')
    for sample in samples:
        header.add_sample(sample)
    return header


def load_predictions(prediction_paths, locus, samples):
    predictions = {}
    for sample in samples:
        dirname = prediction_paths.format(locus=locus, sample=sample)
        with gzip.open(os.path.join(dirname, 'res.json.gz'), 'rt') as inp:
            res = json.load(inp)
            pred = {}
            pred['GT'] = res['genotype']
            gq = np.floor(float(res['quality']) * 10) / 10
            if 'warnings' in res:
                pred['GQ'] = 0.0
                pred['GQ0'] = gq
                pred['WARN'] = ';'.join(w.split('(', 1)[0] for w in res['warnings'])
            else:
                pred['GQ'] = gq
            predictions[sample] = pred
    return predictions


def copy_genotype(var, pred_gt):
    out_gt = []
    for allele in pred_gt.split(','):
        if '.' in allele:
            sample, hap = allele.rsplit('.', 1)
            var_gt = var.samples[sample]['GT']
            out_gt.append(None if var_gt is None else var_gt[int(hap) - 1])
        elif allele == 'GRCh38':
            out_gt.append(0)
        else:
            sample = allele
            var_gt = var.samples[sample]['GT']
            assert var_gt is None or len(var_gt) == 1
            out_gt.append(None if var_gt is None else var_gt[0])
    return out_gt


def process_locus(variants, reference, samples, locus_tup, prediction_paths, output_paths):
    chrom, start, end, locus = locus_tup
    predictions = load_predictions(prediction_paths, locus, samples)

    out_filename = output_paths.format(locus=locus)
    header = create_vcf_header(chrom, reference, samples)
    with pysam.VariantFile(out_filename, 'w', header=header) as out_vcf:
        for var in variants.fetch(chrom, start, end):
            newvar = out_vcf.new_record()
            newvar.chrom = chrom
            newvar.start = var.start
            newvar.alleles = var.alleles
            for sample in samples:
                pred = predictions[sample]
                fmt = newvar.samples[sample]
                fmt['GT'] = copy_genotype(var, pred['GT'])
                fmt.phased = True
                for key in ('GQ', 'GQ0', 'WARN'):
                    if key in pred:
                        fmt[key] = pred[key]
            out_vcf.write(newvar)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--variants', metavar='FILE', required=True,
        help='Indexed phased variants file.')
    parser.add_argument('-p', '--predictions', metavar='FILE', required=True,
        help='Paths to all predictions. Must include patterns `{locus}` and `{sample}`, '
            'which will be replaced with appropriate names.')
    parser.add_argument('-S', '--samples', metavar='FILE', required=True,
        help='List of samples.')
    parser.add_argument('-L', '--loci', metavar='FILE', required=True,
        help='BED file with loci.')
    parser.add_argument('-r', '--reference', metavar='FILE', required=True,
        help='Reference FASTA file.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Paths to output VCF files. Must include `{locus}` pattern.')
    args = parser.parse_args()

    reference = pysam.FastaFile(args.reference)
    variants = pysam.VariantFile(args.variants)
    with open(args.samples) as inp:
        samples = list(map(str.strip, inp))
    with open(args.loci) as inp:
        loci = list(map(str.strip, inp))

    for i, line in enumerate(loci, 1):
        chrom, start, end, locus = line.split('\t')
        start = int(start)
        end = int(end)
        sys.stderr.write('[{:3}/{}] {}\n'.format(i, len(loci), locus))
        process_locus(variants, reference, samples, (chrom, start, end, locus), args.predictions, args.output)


if __name__ == '__main__':
    main()
