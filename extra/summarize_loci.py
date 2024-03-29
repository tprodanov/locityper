#!/usr/bin/env python3

import argparse
import sys
import os
import gzip
from Bio import SeqIO
import collections
import numpy as np
import intervaltree
import tqdm

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
    multiplicities = []
    for seq in seqs:
        counter = collections.Counter()
        for i in range(len(seq) - k + 1):
            counter[seq[i : i + k]] += 1
        mult = sum(counter.values()) / len(counter)
        multiplicities.append(mult)
    return np.mean(multiplicities)


def ling_compl(seq, divisor):
    u1 = len(set(seq))
    u2 = len(set(zip(seq, seq[1:])))
    u3 = len(set(zip(seq, seq[1:], seq[2:])))
    return u1 * u2 * u3 / divisor


def average_ling_compl(seqs, w):
    complexities = []
    divisor = min(4, w) * min(16, w - 1) * min(64, w - 2)
    for seq in seqs:
        n = len(seq)
        count = n // w
        m = count * w
        shift = (n - m) // 2
        for i in range(count):
            window = seq[shift + i * w : shift + i * w + w]
            assert len(window) == w
            complexities.append(ling_compl(window, divisor))
    return np.mean(complexities)


def calc_dupl_frac(trees, chrom, start, end):
    sum_overl = 0
    for dup_start, dup_end, (dup_chrom2, dup_start2, dup_end2) in trees[chrom].overlap(start, end):
        size = min(end, dup_end) - max(start, dup_start)
        if dup_chrom2 == chrom:
            size -= max(0, min(end, dup_end2) - max(start, dup_start2))
        sum_overl += max(0, size)
    return sum_overl / (end - start)


def process_locus(locus_dir, locus, out, trees, args):
    reader = SeqIO.parse(gzip.open(os.path.join(locus_dir, 'haplotypes.fa.gz'), 'rt'), 'fasta')
    seqs = [str(rec.seq) for rec in reader]
    l = sum(len(seq) for seq in seqs) / len(seqs)
    mult = kmer_multiplicity(seqs, args.kmer)
    compl = average_ling_compl(seqs, args.window)

    with open(os.path.join(locus_dir, 'ref.bed')) as f:
        chrom, start, end, _ = next(f).split('\t')
        overl = calc_dupl_frac(trees, chrom, int(start), int(end))

    distances = Distances(os.path.join(locus_dir, 'discarded_haplotypes.txt'),
        os.path.join(locus_dir, 'haplotypes.paf.gz'))
    mean_edit, mean_div = distances.average_divergence()
    out.write(f'{locus}\t{len(seqs)}\t{l:.0f}\t{mult:.9f}\t{compl:.9f}\t{overl:.9f}'
        f'\t{mean_edit:.9f}\t{mean_div:.9f}\n')


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
    args = parser.parse_args()

    with common.open(args.segdup) as f:
        trees = load_segdups(f)

    indir = os.path.join(args.database, 'loci')
    with common.open(args.output, 'w') as out:
        out.write(f'# {" ".join(sys.argv)}\n')
        out.write('locus\tn_seqs\tmean_len\tmean_multiplicity\tmean_complexity\tsegdup_frac\tmean_edit\tmean_div\n')
        for locus in tqdm.tqdm(os.listdir(indir)):
            subdir = os.path.join(indir, locus)
            if not os.path.exists(os.path.join(subdir, 'success')):
                continue
            process_locus(subdir, locus, out, trees, args)


if __name__ == '__main__':
    main()
