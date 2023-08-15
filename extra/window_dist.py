#!/usr/bin/env python3

import argparse
import re
import numpy as np
from tqdm import tqdm

from common import open_stream


PATTERN = re.compile(r'[0-9]+=')

def count_matching_windows(cigar, windows):
    """
    Based on the CIGAR string and the window size `w`.
    """
    counts = [0] * len(windows)
    for entry in PATTERN.findall(cigar):
        length = int(entry[:-1])
        for i, w in enumerate(windows):
            counts[i] += max(0, length - w + 1)
    return counts


def add_dist(line, windows):
    if line.startswith('#'):
        return line
    line = line.strip()
    spl_line = line.split('\t')

    len1 = int(spl_line[1])
    len2 = int(spl_line[6])
    cigar = next(s for s in spl_line if s.startswith('cg:'))[5:]
    overls = count_matching_windows(cigar, windows)

    for window, overl in zip(windows, overls):
        size1 = len1 - window + 1
        size2 = len2 - window + 1
        jaccard_dist = 1.0 - overl / (size1 + size2 - overl)
        line += f'\tm{window}:i:{overl}\td{window}:f:{jaccard_dist:.8f}'
    line += '\n'
    return line


def main():
    parser = argparse.ArgumentParser(
        description='Calculate alignment distance as the number/fraction of matching windows.')
    parser.add_argument('paf', metavar='FILE',
        help='Input PAF[.gz] file with pairwise distances.')
    parser.add_argument('-o', '--output', metavar='FILE', required=False,
        help='Output PAF[.gz] file [default: stdout].')
    parser.add_argument('-w', '--window', metavar='INT', type=int, nargs='+',
        help='Window sizes.')
    args = parser.parse_args()

    with open_stream(args.paf) as inp, open_stream(args.output, 'w') as out:
        for line in tqdm(inp):
            out.write(add_dist(line, args.window))


if __name__ == '__main__':
    main()
