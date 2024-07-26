#!/usr/bin/env python3

import argparse
import os
import sys
import gzip
import glob
import json
import multiprocessing
import math
import numpy as np

import common


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
                common.error(f'Path `{path}` does not contain "." component')
            elif len(ixs) > 1:
                common.error(f'Path `{path}` contains "." component {len(ixs)} times')
            i = ixs[0]
            if i == len(components):
                common.error(f'Cannot get sample name from path `{path}`')
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
        common.error(f'Loaded zero samples')

    samples = set()
    for sample, path in paths:
        if sample in samples:
            common.error(f'Sample {sample} appears twice')
        elif not os.path.isdir(os.path.join(path, 'loci')):
            common.error(f'Directory `{path}` does not contain "loci" subdirectory')
        samples.add(sample)
    sys.stderr.write(f'Loaded {len(paths)} samples\n')
    return paths


def get_or_nan(res, key):
    val = res.get(key)
    return np.nan if val is None else val


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
                s += '*\n'
                continue

            try:
                gt = res['genotype']
                qual = math.floor(10 * float(res['quality'])) * 0.1
                total_reads = get_or_nan(res, 'total_reads')
                weight_dist = get_or_nan(res, 'weight_dist')
                unexpl_reads = get_or_nan(res, 'unexpl_reads')
                warnings = ';'.join(res.get('warnings', '*'))
                s += f'{gt}\t{qual:.1f}\t{total_reads}\t{unexpl_reads}\t{weight_dist:.5f}\t{warnings}\n'
            except:
                sys.stderr.write(f'Error analysing {json_filename}\n')
                raise
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

    input_paths = load_input(args)
    total = len(input_paths)
    finished = 0

    out = common.open(args.output, 'w')
    out.write('# {}\n'.format(' '.join(sys.argv)))
    out.write('sample\tlocus\tgenotype\tquality\ttotal_reads\tunexpl_reads\tweight_dist\twarnings\n')

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
