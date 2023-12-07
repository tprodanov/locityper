#!/usr/bin/env python3

import argparse
import re
import numpy as np
from tqdm import tqdm

import common


def get_count(size, window):
    """
    Based on the region size and the window size, returns number of windows in the region:
    (number of overlapping windows, number of non-overlapping windows).
    """
    return max(0, size - window + 1), size // window


def jaccard_dist(overlap, size1, size2):
    return 1.0 - overlap / (size1 + size2 - overlap)


PATTERN = re.compile(r'[0-9]+=')

def count_matching_windows(cigar, windows):
    """
    Based on the CIGAR string and the window size `w`, returns the number of matching overlapping/non-overlapping
    windows of each size in `windows`.
    """
    counts = np.zeros((len(windows), 2), dtype=int)
    for entry in PATTERN.findall(cigar):
        length = int(entry[:-1])
        for i, w in enumerate(windows):
            counts[i] += get_count(length, w)
    return counts


def add_dist(line, windows):
    if line.startswith('#'):
        return line
    line = line.strip()
    spl_line = line.split('\t')

    len1 = int(spl_line[1])
    len2 = int(spl_line[6])
    cigar = next(s for s in spl_line if s.startswith('cg:'))[5:]
    matching_counts = count_matching_windows(cigar, windows)

    for window, (ov_m, nov_m) in zip(windows, matching_counts):
        ov1, nov1 = get_count(len1, window)
        ov2, nov2 = get_count(len2, window)
        line += f'\tm{window}:i:{ov_m}\td{window}:f:{jaccard_dist(ov_m, ov1, ov2):.8f}'
        line += f'\tM{window}:i:{nov_m}\tD{window}:f:{jaccard_dist(nov_m, nov1, nov2):.8f}'
    line += '\n'
    return line


def main():
    parser = argparse.ArgumentParser(
        description='Calculate alignment distance as the number/fraction of matching windows.')
    parser.add_argument('paf', metavar='FILE',
        help='Input PAF[.gz] file with pairwise distances.')
    parser.add_argument('-o', '--output', metavar='FILE', required=False,
        help='Output PAF[.gz] file [default: stdout].')
    parser.add_argument('-w', '--windows', metavar='INT', type=int, nargs='+',
        help='Window sizes.')
    args = parser.parse_args()

    with common.open(args.paf) as inp, common.open(args.output, 'w') as out:
        for line in tqdm(inp):
            out.write(add_dist(line, args.windows))


if __name__ == '__main__':
    main()
