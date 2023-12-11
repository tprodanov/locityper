#!/usr/bin/env python3

import numpy as np
import argparse
import sys
import pysam
import json
import gzip
import os
import collections
import glob

import common


def error(msg):
    sys.stderr.write(f'ERROR: {msg}\n')
    exit(1)


def load_input(args):
    paths = []
    if args.input is not None:
        in_dirs = args.input
        if len(in_dirs) == 1 and '*' in in_dirs[0]:
            paths = glob.glob(in_dirs[0])
        for path in in_dirs:
            components = path.split('/')
            ixs = [i for i, comp in enumerate(components) if comp == '.']
            if not ixs:
                error(f'Path `{path}` does not contain "." component')
            elif len(ixs) > 1:
                error(f'Path `{path}` contains "." component {len(ixs)} times')
            i = ixs[0]
            if i == len(components):
                error(f'Cannot get sample name from path `{path}`')
            sample = components[i + 1]
            paths.append((sample, path))
    else:
        with common.open(args.input_list) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                path, sample = line.strip().split()
                paths.append((sample, path))
    if not paths:
        error(f'Loaded zero samples')

    samples = set()
    for sample, path in paths:
        if sample in samples:
            error(f'Sample {sample} appears twice')
        elif not os.path.isdir(os.path.join(path, 'loci')):
            error(f'Directory `{path}` does not contain "loci" subdirectory')
        samples.add(sample)
    sys.stderr.write(f'Loaded {len(paths)} samples\n')
    return paths


def load_database(path, contig_lengths):
    loci_dir = os.path.join(path, 'loci')
    if not os.path.isdir(loci_dir):
        error(f'Directory `{path}` does not contain "loci" subdirectory')
    loci = []

    for entry in os.scandir(os.path.join(path, 'loci')):
        if not entry.is_dir():
            continue
        locus = entry.name
        locus_dir = entry.path

        with open(os.path.join(locus_dir, 'ref.bed')) as f:
            chrom, start, end, name = next(f).strip().split()
        if name != locus:
            error(f'Locus directory `{locus_dir}` contains locus `{name}` (expected `{locus}`)')
        start = int(start)
        end = int(end)
        if end > contig_lengths.get(chrom, end):
            error(f'Variant file does is missing chromosome {chrom}, or it is too short')
        loci.append((locus, chrom, start, end))
    sys.stderr.write(f'Loaded {len(loci)} loci\n')
    return loci


def create_vcf_header(chrom, length, input_paths):
    header = pysam.VariantHeader()
    header.add_line('##contig=<ID={},length={}>'.format(chrom, length))
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">')
    header.add_line('##FORMAT=<ID=GQ0,Number=1,Type=Float,'
        'Description="Initial genotype quality, not accounting for warnings">')
    header.add_line('##FORMAT=<ID=WARN,Number=1,Type=String,Description="Genotype warnings">')
    for sample, _ in input_paths:
        header.add_sample(sample)
    return header


def load_predictions(input_paths, locus):
    predictions = {}
    for sample, dirname in input_paths:
        with gzip.open(os.path.join(dirname, 'loci', locus, 'res.json.gz'), 'rt') as inp:
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


def copy_genotype(var, pred_gt, genome_names):
    out_gt = []
    for allele in pred_gt.split(','):
        if len(allele) > 2 and (allele[-2] == '.' or allele[-2] == '_'):
            sample = allele[:-2]
            hap = int(allele[-1])
            var_gt = var.samples[sample]['GT']
            out_gt.append(None if var_gt is None else var_gt[int(hap) - 1])
        elif allele in genome_names:
            out_gt.append(0)
        else:
            sample = allele
            var_gt = var.samples[sample]['GT']
            assert var_gt is None or len(var_gt) == 1
            out_gt.append(None if var_gt is None else var_gt[0])
    return out_gt


def process_locus(variants, contig_lengths, input_paths, locus, genome_names, output):
    locus, chrom, start, end = locus
    predictions = load_predictions(input_paths, locus)
    out_filename = os.path.join(output, f'{locus}.vcf.gz')
    header = create_vcf_header(chrom, contig_lengths[chrom], input_paths)

    with pysam.VariantFile(out_filename, 'wz', header=header) as out_vcf:
        for var in variants.fetch(chrom, start, end):
            newvar = out_vcf.new_record()
            newvar.chrom = chrom
            newvar.start = var.start
            newvar.alleles = var.alleles
            newvar.id = var.id
            for sample, pred in predictions.items():
                fmt = newvar.samples[sample]
                fmt['GT'] = copy_genotype(var, pred['GT'], genome_names)
                fmt.phased = True
                for key in ('GQ', 'GQ0', 'WARN'):
                    if key in pred:
                        fmt[key] = pred[key]
            out_vcf.write(newvar)


def main():
    parser = argparse.ArgumentParser(description='Convert Locityper predictions into VCF file.',
        formatter_class=argparse.RawTextHelpFormatter)
    in_group = parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument('-i', '--input', metavar='DIR', nargs='+',
        help='Path(s) to Locityper analyses. Each directory must contain a `loci` subdir.\n'
            'Please include "." before the sample name (..././<sample>/...).')
    in_group.add_argument('-I', '--input-list', metavar='FILE',
        help='File with two columns: path to Locityper analysis and sample name.\n'
            'Mutually exclusive with `-i/--input`.')
    parser.add_argument('-v', '--variants', metavar='FILE', required=True,
        help='Indexed phased variants file.\n'
            'Must have matching sample names with the loci alleles.')
    parser.add_argument('-d', '--database', metavar='FILE', required=True,
        help='Path to Locityper database.')
    parser.add_argument('-o', '--output', metavar='DIR', required=True,
        help='Output directory.')
    parser.add_argument('-g', '--genome', metavar='STR',
        help='Reference genome name (if the script cannot guess).')
    args = parser.parse_args()

    try:
        os.mkdir(args.output)
    except FileExistsError:
        pass

    input_paths = load_input(args)
    variants = pysam.VariantFile(args.variants)
    contig_lengths = { contig.name: contig.length for contig in variants.header.contigs.itervalues() }
    loci = load_database(args.database, contig_lengths)

    if args.genome:
        genome_names = { args.genome }
    else:
        genome_names = { 'GRCh38', 'GRCh37', 'hg38', 'hg19', 'CHM13' }

    for i, locus in enumerate(loci, 1):
        sys.stderr.write('[{:3}/{}] {}\n'.format(i, len(loci), locus[0]))
        process_locus(variants, contig_lengths, input_paths, locus, genome_names, args.output)


if __name__ == '__main__':
    main()
