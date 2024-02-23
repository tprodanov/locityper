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
import multiprocessing
import functools
import operator

import common


def error(msg):
    sys.stderr.write(f'ERROR: {msg}\n')
    exit(1)


def load_input(args):
    paths = []
    if args.input is not None:
        in_dirs = args.input
        if len(in_dirs) == 1 and '*' in in_dirs[0]:
            in_dirs = glob.glob(in_dirs[0])
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


def load_database(path, subset_loci):
    loci_dir = os.path.join(path, 'loci')
    if not os.path.isdir(loci_dir):
        error(f'Directory `{path}` does not contain "loci" subdirectory')
    loci = []
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
            error(f'Locus directory `{locus_dir}` contains locus `{name}` (expected `{locus}`)')
        start = int(start)
        end = int(end)
        # if chrom not in contigs or end > contigs[chrom].length:
        #     error(f'Chromosome {chrom} is missing or it is too short in the input VCF file')
        loci.append((locus, chrom, start, end))
    sys.stderr.write(f'Loaded {len(loci)} loci\n')
    return loci


def create_vcf_header(chrom_lengths, samples):
    header = pysam.VariantHeader()
    for chrom, length in chrom_lengths:
        header.add_line('##contig=<ID={},length={}>'.format(chrom, length))
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">')
    header.add_line('##FORMAT=<ID=GQ0,Number=1,Type=Float,'
        'Description="Initial genotype quality, not accounting for warnings">')
    header.add_line('##FORMAT=<ID=WARN,Number=1,Type=String,Description="Genotype warnings">')
    for sample in samples:
        header.add_sample(sample)
    return header


def load_predictions(input_paths, locus):
    predictions = {}
    for sample, dirname in input_paths:
        json_filename = os.path.join(dirname, 'loci', locus, 'res.json.gz')
        with gzip.open(json_filename, 'rt') as inp:
            try:
                res = json.load(inp)
            except json.decoder.JSONDecodeError:
                sys.stderr.write(f'Cannot parse json from {json_filename}\n')
                raise
            pred = {}
            if 'genotype' not in res:
                pred = None
                continue

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


def copy_genotype(var, pred_gt, genome_name):
    out_gt = []
    for allele in pred_gt.split(','):
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


def process_locus(locus, input_paths, genome_name, output):
    locus, chrom, start, end = locus
    predictions = load_predictions(input_paths, locus)
    out_filename = os.path.join(output, f'{locus}.vcf.gz')

    global in_vcf
    header = create_vcf_header(((chrom, in_vcf.header.contigs[chrom].length),),
        map(operator.itemgetter(0), input_paths))
    with pysam.VariantFile(out_filename, 'wz', header=header) as out_vcf:
        for var in in_vcf.fetch(chrom, start, end):
            newvar = out_vcf.new_record()
            newvar.chrom = chrom
            newvar.start = var.start
            newvar.alleles = var.alleles
            newvar.id = var.id
            for sample, pred in predictions.items():
                fmt = newvar.samples[sample]
                fmt['GT'] = copy_genotype(var, pred['GT'], genome_name)
                fmt.phased = True
                for key in ('GQ', 'GQ0', 'WARN'):
                    if key in pred:
                        fmt[key] = pred[key]
            out_vcf.write(newvar)
    pysam.tabix_index(out_filename, preset='vcf')
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
            gq1 = (fmt1['GQ'], fmt1.get('GQ0', 0))
            gq2 = (fmt2['GQ'], fmt2.get('GQ0', 0))
            if gq1 < gq2:
                fmt1['GT'] = fmt2['GT']


def merge_vcfs(in_vcf, input_paths, out_dir, loci):
    chroms = set(map(operator.itemgetter(1), loci))
    chroms = { chrom: in_vcf.get_tid(chrom) for chrom in chroms }
    chrom_lengths = [(chrom, in_vcf.header.contigs[chrom].length)
        for chrom, _ in sorted(chroms.items(), key=operator.itemgetter(1))]

    n_samples = len(input_paths)
    records = {}
    loci.sort(key=lambda tup: (chroms[tup[1]], tup[2], tup[3]))
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

    header = create_vcf_header(chrom_lengths, map(operator.itemgetter(0), input_paths))
    out_filename = os.path.join(out_dir, 'merged.vcf.gz')
    with pysam.VariantFile(out_filename, 'wz', header=header) as out_vcf:
        for key in sorted(records.keys()):
            rec = records[key]
            rec.translate(header)
            out_vcf.write(rec)
    pysam.tabix_index(out_filename, preset='vcf')


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
    parser.add_argument('-d', '--database', metavar='FILE', required=True,
        help='Path to Locityper database.')
    parser.add_argument('-v', '--variants', metavar='FILE', required=True,
        help='Indexed phased variants file.\n'
            'Must have matching sample names with the loci alleles.')
    parser.add_argument('-o', '--output', metavar='DIR', required=True,
        help='Output directory.')
    parser.add_argument('-g', '--genome', metavar='STR', required=True,
        help='Reference genome name.')
    parser.add_argument('-@', '--threads', metavar='INT', type=int, default=8,
        help='Analyze loci in this many threads [%(default)s].')
    parser.add_argument('--subset-loci', metavar='STR', nargs='+',
        help='Limit the analysis to these loci.')
    args = parser.parse_args()

    common.mkdir(args.output)
    input_paths = load_input(args)
    loci = load_database(args.database, args.subset_loci)

    total = len(loci)
    finished = 0
    def callback(locus):
        nonlocal finished
        finished += 1
        sys.stderr.write(f'Finished [{finished:3}/{total}] {locus}\n')

    with multiprocessing.Pool(args.threads, initializer=create_thread, initargs=(args.variants,)) as pool:
        results = [pool.apply_async(process_locus, (locus, input_paths, args.genome, args.output), callback=callback)
            for locus in loci]
        for res in results:
            res.get()
        pool.close()
        pool.join()

    sys.stderr.write('Merging variants\n')
    merge_vcfs(pysam.VariantFile(args.variants), input_paths, args.output, loci)


if __name__ == '__main__':
    main()
