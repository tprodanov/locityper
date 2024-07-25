#!/usr/bin/env python3

import numpy as np
import argparse
import sys
import pysam
import os
from collections import defaultdict
import multiprocessing
import operator
from simpleeval import simple_eval

import common


def parse_float(val):
    return float(val) if val != 'NA' else np.nan


def load_predictions(f, expr):
    by_locus = defaultdict(dict)
    samples = set()
    for row in common.read_csv(f):
        locus = row['locus']
        sample = row['sample']

        genotype = row['genotype'].split(',')
        features = dict(
            qual=parse_float(row['quality']),
            reads=int(row['total_reads']),
            unexpl=int(row['unexpl_reads']),
            wdist=parse_float(row['weight_dist']),
            warn=row['warnings'],
        )
        passes = bool(simple_eval(expr, names=features))
        features['GQ'] = int(passes)
        by_locus[locus][sample] = (genotype, features)
        samples.add(sample)
    return by_locus, sorted(samples)


def load_database(path, subset_loci):
    loci_dir = os.path.join(path, 'loci')
    if not os.path.isdir(loci_dir):
        common.error(f'Directory `{path}` does not contain "loci" subdirectory')
    loci = {}
    if subset_loci is not None:
        subset_loci = set(subset_loci)

    for entry in os.scandir(os.path.join(path, 'loci')):
        if not entry.is_dir():
            continue
        locus = entry.name
        if subset_loci is not None and locus not in subset_loci:
            continue

        locus_dir = entry.path
        with open(os.path.join(locus_dir, 'ref.bed')) as f:
            chrom, start, end, name = next(f).strip().split()
        if name != locus:
            common.error(f'Locus directory `{locus_dir}` contains locus `{name}` (expected `{locus}`)')
        start = int(start)
        end = int(end)
        # if chrom not in contigs or end > contigs[chrom].length:
        #     error(f'Chromosome {chrom} is missing or it is too short in the input VCF file')
        loci[locus]= (chrom, start, end)
    sys.stderr.write(f'Loaded {len(loci)} loci\n')
    return loci


def create_vcf_header(chrom_lengths, samples):
    header = pysam.VariantHeader()
    for chrom, length in chrom_lengths:
        header.add_line('##contig=<ID={},length={}>'.format(chrom, length))
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="0/1 value based on prediction filtering">')
    header.add_line('##FORMAT=<ID=qual,Number=1,Type=Float,Description="Locus genotype quality">')
    header.add_line('##FORMAT=<ID=reads,Number=1,Type=Integer,'
        'Description="Total number of reads used to predict locus genotype">')
    header.add_line('##FORMAT=<ID=unexpl,Number=1,Type=Integer,'
        'Description="Number of unexplained reads within the locus">')
    header.add_line('##FORMAT=<ID=wdist,Number=1,Type=Float,'
        'Description="Weighted distance between primary and non-primary genotype predictions">')
    header.add_line('##FORMAT=<ID=warn,Number=1,Type=String,Description="Genotype warnings">')
    for sample in samples:
        header.add_sample(sample)
    return header


def copy_genotype(pred_gt, var, genome_name):
    out_gt = []
    for allele in pred_gt:
        if len(allele) > 2 and (allele[-2] == '.' or allele[-2] == '_'):
            sample = allele[:-2]
            hap = int(allele[-1])
            var_gt = var.samples[sample]['GT']
            out_gt.append(None if var_gt is None else var_gt[int(hap) - 1])
        elif allele == genome_name:
            out_gt.append(0)
        else:
            sample = allele
            var_gt = var.samples[sample]['GT']
            assert var_gt is None or len(var_gt) == 1
            out_gt.append(None if var_gt is None else var_gt[0])
    return out_gt


def process_locus(locus, coords, samples, preds, genome_name, output):
    chrom, start, end = coords
    out_filename = os.path.join(output, f'{locus}.vcf.gz')

    global in_vcf
    chrom_len = in_vcf.header.contigs[chrom].length
    header = create_vcf_header(((chrom, chrom_len),), samples)
    with pysam.VariantFile(out_filename, 'wz', header=header) as out_vcf:
        for var in in_vcf.fetch(chrom, start, end):
            newvar = out_vcf.new_record()
            newvar.chrom = chrom
            newvar.start = var.start
            newvar.alleles = var.alleles
            newvar.id = var.id
            for i, sample in enumerate(samples):
                fmt = newvar.samples[i]
                pred = preds.get(sample)
                if pred is None:
                    continue

                fmt['GT'] = copy_genotype(pred[0], var, genome_name)
                fmt.phased = True
                for key, val in pred[1].items():
                    fmt[key] = val
            out_vcf.write(newvar)
    pysam.tabix_index(out_filename, preset='vcf', force=True)
    return locus


def create_thread(vcf_filename):
    global in_vcf
    in_vcf = pysam.VariantFile(vcf_filename)


def merge_vars(rec1, rec2, n_samples):
    assert rec1.id == rec2.id
    for sample_ix in range(n_samples):
        fmt1 = rec1.samples[sample_ix]
        fmt2 = rec2.samples[sample_ix]
        if fmt1['GT'] != fmt2['GT']:
            gq1 = (fmt1['GQ'], fmt1['qual'])
            gq2 = (fmt1['GQ'], fmt2['qual'])
            if gq1 < gq2:
                fmt1.clear()
                for key, val in fmt2.items():
                    fmt1[key] = val


def merge_vcfs(in_vcf, samples, out_dir, loci):
    chroms = set(chrom for chrom, _, _ in loci.values())
    chroms = { chrom: in_vcf.get_tid(chrom) for chrom in chroms }
    chrom_lengths = [(chrom, in_vcf.header.contigs[chrom].length)
        for chrom, _ in sorted(chroms.items(), key=operator.itemgetter(1))]

    n_samples = len(samples)
    records = {}
    loci = sorted(loci.items(), key=lambda tup: (chroms[tup[1][0]], tup[1][1], tup[1][2]))
    total = 0
    for locus in map(operator.itemgetter(0), loci):
        with pysam.VariantFile(os.path.join(out_dir, f'{locus}.vcf.gz')) as vcf:
            for record in vcf:
                key = (chroms[record.chrom], record.start, record.alleles)
                if key in records:
                    merge_vars(records[key], record, n_samples)
                    total += 1
                else:
                    records[key] = record
    sys.stderr.write(f'    Merged {total} variants\n')

    header = create_vcf_header(chrom_lengths, samples)
    out_filename = os.path.join(out_dir, 'merged.vcf.gz')
    with pysam.VariantFile(out_filename, 'wz', header=header) as out_vcf:
        for key in sorted(records.keys()):
            rec = records[key]
            rec.translate(header)
            out_vcf.write(rec)
    pysam.tabix_index(out_filename, preset='vcf', force=True)


def main():
    parser = argparse.ArgumentParser(description='Convert Locityper predictions into VCF file.',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='CSV file with Locityper predictions.')
    parser.add_argument('-d', '--database', metavar='FILE', required=True,
        help='Path to Locityper database.')
    parser.add_argument('-v', '--variants', metavar='FILE', required=True,
        help='Indexed phased variants file.\n'
            'Must have matching sample names with the loci alleles.')
    parser.add_argument('-o', '--output', metavar='DIR', required=True,
        help='Output directory.')
    parser.add_argument('-g', '--genome-name', metavar='STR', required=True,
        help='Reference genome name.')
    parser.add_argument('-@', '--threads', metavar='INT', type=int, default=8,
        help='Analyze loci in this many threads [%(default)s].')

    DEF_EXPR = 'warn == "*" and wdist < 30 and (unexpl < 1000 or unexpl < reads * 0.2)'
    parser.add_argument('-f', '--filtering', metavar='STR', default=DEF_EXPR,
        help='Simple expression to determine if the locus passes filterin. Default = `%(default)s`.')
    parser.add_argument('--subset-loci', metavar='STR', nargs='+',
        help='Limit the analysis to these loci.')
    args = parser.parse_args()

    common.mkdir(args.output)
    loci = load_database(args.database, args.subset_loci)
    with common.open(args.input) as f:
        preds, samples = load_predictions(f, args.filtering)
    sys.stderr.write(f'Loaded {len(samples)} samples\n')
    loci_inters = set(loci.keys()) & set(preds.keys())
    if len(loci_inters) < len(loci):
        sys.stderr.write('WARN: Genotype predictions missing for {} loci, such as {}\n'.format(
            len(loci) - len(loci_inters), list(set(loci.keys()) - loci_inters)[:5].join(', ')))

    total = len(loci_inters)
    finished = 0
    def callback(locus):
        nonlocal finished
        finished += 1
        sys.stderr.write(f'Finished [{finished:3}/{total}] {locus}\n')

    with multiprocessing.Pool(args.threads, initializer=create_thread, initargs=(args.variants,)) as pool:
        results = [pool.apply_async(process_locus,
            (locus, loci[locus], samples, preds[locus], args.genome_name, args.output), callback=callback)
            for locus in loci_inters]
        for res in results:
            res.get()
        pool.close()
        pool.join()

    sys.stderr.write('Merging variants\n')
    merge_vcfs(pysam.VariantFile(args.variants), samples, args.output, loci)


if __name__ == '__main__':
    main()
