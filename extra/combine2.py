#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.special import logsumexp
import gzip

import common


def load_sol(path):
    sol = pd.read_csv(path, sep='\t', comment='#')
    max_stage = sol.stage.max()
    sol = sol[sol.stage == max_stage].groupby('genotype').mean().reset_index()
    return sol[['genotype', 'lik']]


def process_locus(dir1, dir2, out_dir):
    sol1 = load_sol(os.path.join(dir1, 'sol.csv.gz'))
    sol2 = load_sol(os.path.join(dir2, 'sol.csv.gz'))
    solj = pd.merge(sol1, sol2, how='inner', on='genotype')
    solj['lik'] = solj.lik_x + solj.lik_y
    solj.sort_values(by='lik', inplace=True, ascending=False)
    solj.to_csv(os.path.join(out_dir, 'sol.csv.gz'), sep='\t', index=False)

    liks = np.array(solj.lik) * np.log(10)
    oth_prob = logsumexp(liks[1:])
    oth_prob = oth_prob - logsumexp((liks[0], oth_prob))
    phred = max(-10 * oth_prob / np.log(10), 0.01)
    with gzip.open(os.path.join(out_dir, 'res.json.gz'), 'wt') as out:
        out.write('{\n')
        out.write('    "locus": "{}"\n'.format(os.path.basename(out_dir)))
        out.write('    "genotype": "{}"\n'.format(solj.genotype[0]))
        out.write('    "quality": {}\n'.format(phred))
        out.write('}\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--input-1', metavar='DIR', required=True,
        help='First input directory.')
    parser.add_argument('-2', '--input-2', metavar='DIR', required=True,
        help='Second input directory.')
    parser.add_argument('-o', '--output', metavar='DIR', required=True,
        help='Output directory.')
    args = parser.parse_args()

    loci1 = set(os.listdir(os.path.join(args.input_1, 'loci')))
    loci2 = set(os.listdir(os.path.join(args.input_2, 'loci')))
    common.mkdir(args.output)
    common.mkdir(os.path.join(args.output, 'loci'))
    loci = sorted(loci1 & loci2)
    for i, locus in enumerate(loci, 1):
        sys.stderr.write('[{:3}/{}] {}\n'.format(i, len(loci), locus))
        out_dir = os.path.join(args.output, 'loci', locus)
        common.mkdir(out_dir)
        process_locus(
            os.path.join(args.input_1, 'loci', locus),
            os.path.join(args.input_2, 'loci', locus),
            out_dir)


if __name__ == '__main__':
    main()
