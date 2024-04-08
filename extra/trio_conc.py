#!/usr/bin/env python3

import sys
import os.path
import argparse
from collections import defaultdict
import tqdm
import numpy as np

import common
import gt_dist


def load_pedigree(f):
    trios = []
    for row in common.read_csv(f):
        indiv = row['Individual ID']
        mother = row['Maternal ID']
        father = row['Paternal ID']
        if mother != '0' and father != '0':
            trios.append((indiv, mother, father))
    return trios


def load_predictions(f):
    reader = common.read_csv(f)
    preds = defaultdict(dict)
    for row in reader:
        sample = row['sample']
        locus = row['locus']
        preds[locus][sample] = row
    return preds


def format_dist(dist, size):
    div = dist / size
    qv = -10 * np.log10(div) if div > 0 else np.inf
    return f'{dist}\t{size}\t{div:.9f}\t{qv:.9f}'


def calc_concordance(gt_indiv, gt_mother, gt_father, distances):
    haps_i = gt_indiv.split(',')
    haps_m = gt_mother.split(',')
    haps_f = gt_father.split(',')
    assert len(haps_m) == len(haps_f) == 2

    best_div = np.inf
    best_edit = None
    best_corr = None
    for i in range(8):
        j1 = bool(i & 1)
        j2 = bool(i & 2)
        j3 = bool(i & 4)
        hap_i1 = haps_i[j1]
        hap_i2 = haps_i[1 - j1]
        hap_m = haps_m[j2]
        hap_f = haps_f[j3]
        dist1, size1 = distances.distances[hap_i1][hap_m]
        dist2, size2 = distances.distances[hap_i2][hap_f]
        dist = dist1 + dist2
        size = size1 + size2
        div = dist / size
        if div < best_div:
            best_div = div
            best_edit = (dist1, size1, dist2, size2)
            best_corr = (hap_i1, hap_m, hap_i2, hap_f)

    dist1, size1, dist2, size2 = best_edit
    hap_i1, hap_m, hap_i2, hap_f = best_corr
    return (
        f'hap\t{format_dist(dist1, size1)}\t{hap_i1}{"~="[hap_i1 == hap_m]}{hap_m}',
        f'hap\t{format_dist(dist2, size2)}\t{hap_i2}{"~="[hap_i2 == hap_f]}{hap_f}',
        f'gt\t{format_dist(dist1 + dist2, size1 + size2)}\t*'
    )


def process(preds, trios, alns_fmt, out):
    skipped_loci = 0
    skipped_trios = 0
    for locus, locus_preds in tqdm.tqdm(preds.items(), total=len(preds)):
        alns_filename = alns_fmt.format(locus)
        if not os.path.exists(alns_filename):
            skipped_loci += 1
            continue

        distances = gt_dist.Distances(None, alns_filename)
        for (indiv, mother, father) in trios:
            pred_indiv = locus_preds.get(indiv)
            pred_mother = locus_preds.get(mother)
            pred_father = locus_preds.get(father)
            if pred_indiv is None or pred_mother is None or pred_father is None:
                skipped_trios += 1
                continue

            gt_indiv = pred_indiv['genotype']
            gt_mother = pred_mother['genotype']
            gt_father = pred_father['genotype']
            if gt_indiv == '*' or gt_mother == '*' or gt_father == '*':
                skipped_trios += 1
                continue

            prefix = f'{locus}\t{indiv}\t{mother}\t{father}\t'
            for s in calc_concordance(gt_indiv, gt_mother, gt_father, distances):
                out.write(f'{prefix}{s}\n')

    if skipped_loci:
        sys.stderr.write(f'Skipped {skipped_loci} loci\n')
    if skipped_trios:
        sys.stderr.write(f'On average, skipped {skipped_trios / len(preds):.3f} trios per locus\n')


def main():
    parser = argparse.ArgumentParser(description='Calculate trio concordance.')
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Input CSV file with Locityper predictions.')
    parser.add_argument('-p', '--pedigree', metavar='FILE', required=True,
        help='CSV file with pedigree information. Required columns `Individual ID`, `Maternal ID` and `Paternal ID`.')
    parser.add_argument('-a', '--alns', metavar='STR', required=True,
        help='Path to the PAF file, where locus name is replaced with `{}`.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    args = parser.parse_args()

    with common.open(args.pedigree) as f:
        trios = load_pedigree(f)
    with common.open(args.input) as f:
        preds = load_predictions(f)
    with common.open(args.output, 'w') as out:
        out.write(f'# {" ".join(sys.argv)}\n')
        out.write('locus\tindiv\tmother\tfather\ttype\tedit\tsize\tdiv\tqv\texpl\n')
        process(preds, trios, args.alns, out)


if __name__ == '__main__':
    main()
