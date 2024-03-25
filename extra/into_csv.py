#!/usr/bin/env python3

import argparse
import os
import sys
import gzip
import json
import multiprocessing
import math

import common
import into_vcf


def process_sample(sample, sample_dir):
    s = ''
    for entry in os.scandir(os.path.join(sample_dir, 'loci')):
        if not entry.is_dir():
            continue
        locus = entry.name

        json_filename = os.path.join(entry.path, 'res.json.gz')
        with gzip.open(json_filename, 'rt') as inp:
            try:
                res = json.load(inp)
            except json.decoder.JSONDecodeError:
                sys.stderr.write(f'Cannot parse json from {json_filename}\n')
                raise
            s += f'{sample}\t{locus}\t'
            if 'genotype' not in res:
                s += '*\tNA\t*\n'
                continue

            gt = res['genotype']
            qual = math.floor(10 * float(res['quality'])) * 0.1
            weight_dist = res.get('weight_dist', np.nan)
            warnings = ';'.join(res.get('warnings', '*'))
            s += f'{gt}\t{qual:.1f}\t{weight_dist:.5f}\t{warnings}\n'
    return s


def main():
    parser = argparse.ArgumentParser(description='Convert Locityper predictions into CSV file.',
        formatter_class=argparse.RawTextHelpFormatter)
    in_group = parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument('-i', '--input', metavar='DIR', nargs='+',
        help='Path(s) to Locityper analyses. Each directory must contain a `loci` subdir.\n'
            'Please include "." before the sample name (..././<sample>/...).')
    in_group.add_argument('-I', '--input-list', metavar='FILE',
        help='File with two columns: path to Locityper analysis and sample name.\n'
            'Mutually exclusive with `-i/--input`.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    parser.add_argument('-@', '--threads', metavar='INT', type=int, default=8,
        help='Number of threads [%(default)s].')
    args = parser.parse_args()

    input_paths = into_vcf.load_input(args)
    total = len(input_paths)
    finished = 0

    out = common.open(args.output, 'w')
    out.write('# {}\n'.format(' '.join(sys.argv)))
    out.write('sample\tlocus\tgenotype\tquality\tweight_dist\twarnings\n')

    def callback(s):
        nonlocal finished
        finished += 1
        out.write(s)
        if finished % 10 == 0 or finished == total:
            sys.stderr.write(f'Finished [{finished:4} / {total}]\n')
            out.flush()

    with multiprocessing.Pool(args.threads) as pool:
        results = [pool.apply_async(process_sample, tup, callback=callback) for tup in input_paths]
        for res in results:
            res.get()
        pool.close()
        pool.join()


if __name__ == '__main__':
    main()
