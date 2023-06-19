#!/usr/bin/env python3

import argparse
import re
import numpy as np
from tqdm import tqdm

from tup_dist import open_stream


def define_harmonic_sums(power, max_size):
    s = 0.0
    sums = [s]
    for i in range(1, max_size + 1):
        s += 1 / i ** power
        sums.append(s)
    return np.array(sums)


PATTERN = re.compile(r'[0-9]+[A-Z=]')

def error_array(cigar):
    """
    Based on the CIGAR string, returns True/False list,
    where True represents positions of errors.
    Additionally, padd array with Falses.
    """
    arr = []
    for entry in PATTERN.findall(cigar):
        length = int(entry[:-1])
        op = entry[-1]
        arr.extend((op != '=',) * length)
    return arr


def calc_dist(arr, window, hsums):
    mult = 1 / window
    dist = 0.0
    s = sum(arr[:window])
    dist += mult * hsums[s]
    for i in range(len(arr) - window):
        s = s + arr[i + window] - arr[i]
        dist += mult * hsums[s]
    return dist


def add_dist(line, window, hsums, padding, max_dist_const):
    if line.startswith('#'):
        return line
    line = line.strip()
    cigar = next(s for s in line.split('\t') if s.startswith('cg:'))[5:]

    arr = error_array(cigar)
    l = len(arr)
    if l <= window:
        dist = hsums[sum(arr)]
        max_dist = hsums[len(arr)]
    else:
        max_dist = max_dist_const + (len(arr) - window + 1) * hsums[-1] / window
        if padding:
            arr = [False] * padding + arr + [False] * padding
        dist = calc_dist(arr, window, hsums)
    diverg = dist / max_dist
    return f'{line}\thd:f:{dist:.8f}\thf:f:{diverg:.8f}\n'


def main():
    parser = argparse.ArgumentParser(
        description='Recalculate alignment distance by penalizing close errors less.')
    parser.add_argument('paf', metavar='FILE',
        help='Input PAF[.gz] file with pairwise distances.')
    parser.add_argument('-o', '--output', metavar='FILE', required=False,
        help='Output PAF[.gz] file [default: stdout].')
    parser.add_argument('-w', '--window', metavar='INT', type=int, default=150,
        help='Window size [default: %(default)s].')
    parser.add_argument('-p', '--power', metavar='FLOAT', type=float, default=1,
        help='Harmonic sum power [default: %(default)s].')
    parser.add_argument('-P', '--padding', metavar='INT', type=int,
        help='Padding size [default: window-1].')
    args = parser.parse_args()

    hsums = define_harmonic_sums(args.power, args.window)
    if args.padding is None:
        args.padding = args.window - 1
    args.padding = min(args.padding, args.window - 1)
    max_dist_const = 2.0 * hsums[:args.padding].sum() / args.window

    with open_stream(args.paf) as inp, open_stream(args.output, 'w') as out:
        for line in tqdm(inp):
            out.write(add_dist(line, args.window, hsums, args.padding, max_dist_const))


if __name__ == '__main__':
    main()
