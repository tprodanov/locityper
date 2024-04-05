#!/usr/bin/env python3

import argparse
import operator
import sys
import os
import gzip
from Bio import SeqIO
import collections
import numpy as np
import intervaltree
import tqdm
import scipy.optimize

import common
from gt_dist import Distances


def load_segdups(bed):
    trees = collections.defaultdict(intervaltree.IntervalTree)
    header = next(bed)
    assert header.startswith('#')
    fields = header.lstrip('#').strip().split('\t')
    chrom_col = fields.index('chrom')
    start_col = fields.index('chromStart')
    end_col = fields.index('chromEnd')
    chrom2_col = fields.index('otherChrom')
    start2_col = fields.index('otherStart')
    end2_col = fields.index('otherEnd')

    for line in bed:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            chrom = line[chrom_col]
            start = int(line[start_col])
            end = int(line[end_col])

            chrom2 = line[chrom2_col]
            start2 = int(line[start2_col])
            end2 = int(line[end2_col])
            trees[chrom].addi(start, end, (chrom2, start2, end2))
    return trees


def kmer_multiplicity(seqs, k):
    sum_multipl = 0.0
    div_multipl = 0
    for seq, count in seqs:
        counter = collections.Counter()
        for i in range(len(seq) - k + 1):
            counter[seq[i : i + k]] += 1
        mult = sum(counter.values()) / len(counter)
        sum_multipl += mult * count
        div_multipl += count
    return sum_multipl / div_multipl


def ling_compl(seq, divisor):
    u1 = len(set(seq))
    u2 = len(set(zip(seq, seq[1:])))
    u3 = len(set(zip(seq, seq[1:], seq[2:])))
    return u1 * u2 * u3 / divisor


def average_ling_compl(seqs, w):
    sum_complexities = 0.0
    div_complexities = 0.0
    sum_low_complex = 0
    div_low_complex = 0
    divisor = min(4, w) * min(16, w - 1) * min(64, w - 2)

    for seq, count in seqs:
        n = len(seq)
        win_count = n // w
        m = win_count * w
        shift = (n - m) // 2
        curr_low_complex = 0
        for i in range(win_count):
            window = seq[shift + i * w : shift + i * w + w]
            assert len(window) == w
            compl = ling_compl(window, divisor)
            sum_complexities += compl * count
            curr_low_complex += compl < 0.5

        sum_low_complex += curr_low_complex * count
        div_low_complex += count
        div_complexities += win_count * count

    return sum_complexities / div_complexities, sum_low_complex / div_low_complex


def calc_dupl_frac(trees, chrom, start, end):
    sum_overl = 0
    for dup_start, dup_end, (dup_chrom2, dup_start2, dup_end2) in trees[chrom].overlap(start, end):
        size = min(end, dup_end) - max(start, dup_start)
        if dup_chrom2 == chrom:
            size -= max(0, min(end, dup_end2) - max(start, dup_start2))
        sum_overl += max(0, size)
    return sum_overl / (end - start)


def fit_heaps(curves):
    n_perm, n_haps = curves.shape
    best_params = None
    best_fit = np.inf

    def sq_err(params):
        k, gamma = params
        if k <= 0.0 or gamma <= -1.0 or gamma >= 1.0:
            return 1e30
        pred = k * np.arange(1, n_haps + 1) ** gamma
        fit = np.sum((curves - pred) ** 2)
        nonlocal best_fit, best_params
        if fit < best_fit:
            best_fit = fit
            best_params = (k, gamma)
        return fit

    x00 = np.mean(curves[:, 0])
    scipy.optimize.minimize(sq_err, x0=(x00, 0.2), method='Nelder-Mead')
    assert best_params is not None
    return best_params


def write_curves(curves, k, gamma, prefix, out):
    n, m = curves.shape
    for i in range(n):
        for j in range(m):
            out.write(f'{prefix}\t{i+1}\t{j + 1}\t{curves[i, j]}\n')

    medians = np.median(curves, axis=0)
    pred = k * np.arange(1, m + 1) ** gamma
    for j in range(m):
        out.write(f'{prefix}\tmedian\t{j+1}\t{medians[j]:.3f}\n')
        out.write(f'{prefix}\theaps\t{j+1}\t{pred[j]:.3f}\n')


def calc_saturation1(prefix, distances, n_perm, out):
    haps = list(distances.lengths.keys())
    n_haps = len(haps)
    curves = np.zeros((n_perm, n_haps))
    for i in range(n_perm):
        perm = np.random.permutation(haps)
        cum_size = 0
        for j, hap in enumerate(perm):
            cum_size += min(distances.dir_distances[hap2][hap] for hap2 in perm[:j]) if j else distances.lengths[hap]
            curves[i, j] = cum_size
    k, gamma = fit_heaps(curves)
    if out:
        write_curves(curves, k, gamma, prefix, out)
    return k, gamma


def calc_saturation2(prefix, seqs, kmer, n_perm, out):
    kmer_sets = []
    for seq, count in seqs:
        curr_kmers = set(hash(seq[i:i+kmer]) for i in range(len(seq) - kmer + 1))
        kmer_sets.extend((curr_kmers,) * count)
    n_haps = len(seqs)
    ixs = np.arange(n_haps)

    curves = np.zeros((n_perm, n_haps))
    for i in range(n_perm):
        perm = np.random.permutation(ixs)
        comb_set = set()
        for j, w in enumerate(perm):
            comb_set |= kmer_sets[w]
            curves[i, j] = len(comb_set)
    k, gamma = fit_heaps(curves)
    if out:
        write_curves(curves, k, gamma, prefix, out)
    return k, gamma


def process_locus(locus_dir, locus, out, sat_out, trees, args):
    reader = SeqIO.parse(gzip.open(os.path.join(locus_dir, 'haplotypes.fa.gz'), 'rt'), 'fasta')
    distances = Distances(
        os.path.join(locus_dir, 'discarded_haplotypes.txt'),
        os.path.join(locus_dir, 'haplotypes.paf.gz'),
        dir_dist=True)

    seqs = [(str(rec.seq), distances.group_size(rec.id)) for rec in reader]
    total_seqs = sum(map(operator.itemgetter(1), seqs))
    l = sum(len(seq) * count for seq, count in seqs) / total_seqs
    mult = kmer_multiplicity(seqs, args.kmer)
    compl, n_low_complex = average_ling_compl(seqs, args.window)
    mean_edit, mean_div = distances.average_divergence()

    with open(os.path.join(locus_dir, 'ref.bed')) as f:
        chrom, start, end, _ = next(f).split('\t')
        overl = calc_dupl_frac(trees, chrom, int(start), int(end))

    sat_k = sat_gamma = kmer_sat_k = kmer_sat_gamma = np.nan
    if args.saturation[0] != 0:
        sat_k, sat_gamma = calc_saturation1(f'{locus}\taln', distances, args.saturation[0], sat_out)
    if args.saturation[1] != 0:
        kmer_sat_k, kmer_sat_gamma = calc_saturation2(f'{locus}\tkmers', seqs,
            args.sat_kmer, args.saturation[1], sat_out)

    out.write(f'{locus}\t{total_seqs}\t{len(seqs)}\t{l:.0f}'
        f'\t{mult:.9f}\t{compl:.9f}\t{n_low_complex:.5f}\t{overl:.9f}\t{mean_edit:.9f}\t{mean_div:.9f}'
        f'\t{sat_k:.5f}\t{1 - sat_gamma:.9f}\t{kmer_sat_k:.5f}\t{1 - kmer_sat_gamma:.9f}\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', metavar='DIR', required=True,
        help='Locityper database. Each locus directory must contain `haplotypes.paf.gz`.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    parser.add_argument('-k', '--kmer', metavar='INT', type=int, default=15,
        help='k-mer size for multiplicity calculation [%(default)s].')
    parser.add_argument('-w', '--window', metavar='INT', type=int, default=100,
        help='Calculate linguistic complexity on windows of this size [%(default)s].')
    parser.add_argument('-s', '--segdup', metavar='FILE', required=True,
        help='BED file with segmental duplications.')
    parser.add_argument('-S', '--saturation', metavar='INT', type=int, nargs=2, default=(200, 50),
        help='Calculate saturation curves using INT permutations [%(default)s]. Use 0 to skip.'
            ' Two values required: for sequence & k-mers saturation calculation.')
    parser.add_argument('--sat-kmer', metavar='INT', type=int, default=25,
        help='k-mer size for saturation calculation [%(default)s].')
    parser.add_argument('-O', '--sat-out', metavar='FILE', required=False,
        help='Optional: output saturation curves into this file.')
    args = parser.parse_args()

    with common.open(args.segdup) as f:
        trees = load_segdups(f)

    indir = os.path.join(args.database, 'loci')

    with common.open(args.output, 'w') as out:
        out.write(f'# {" ".join(sys.argv)}\n')
        out.write('locus\tn_seqs\tunique_seqs\tmean_len\tmean_multiplicity\tmean_complexity\tn_low_complex\tsegdup_frac'
            '\tmean_edit\tmean_div\tsat_k\tsat_alpha\tkmer_sat_k\tkmer_sat_alpha\n')

        if args.sat_out and max(args.saturation) > 0:
            sat_out = common.open(args.sat_out, 'w')
            sat_out.write('locus\ttype\titeration\tn_haps\tsize\n')
        else:
            sat_out = None

        for locus in tqdm.tqdm(os.listdir(indir)):
            subdir = os.path.join(indir, locus)
            if not os.path.exists(os.path.join(subdir, 'success')):
                continue
            process_locus(subdir, locus, out, sat_out, trees, args)


if __name__ == '__main__':
    main()
