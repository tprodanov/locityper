#!/usr/bin/env python3

import argparse
import os
import sys

import common


def load_haplogroups(f):
    fields = next(f).strip().split()[1:]
    groups = {}
    for line in f:
        line = line.strip().split()
        groups[line[0]] = tuple(line[1:])
    return fields, groups


def get_gt_groups(gt, groups):
    hap1, hap2 = gt.split(',')
    return tuple(map(','.join, map(sorted, zip(groups[hap1], groups[hap2]))))


def add_columns(inp, out, group_fields, groups):
    for line in inp:
        if line.startswith('#'):
            out.write(line)
        else:
            s = line.strip()
            fields = s.split('\t')
            true_field = fields.index('closest_gts')
            pred_field = fields.index('ml_gt')
            s += '\t{}\t{}\n'.format(
                '\t'.join(map('true_{}'.format, group_fields)),
                '\t'.join(map('pred_{}'.format, group_fields)))
            out.write(s)
            break

    n = len(group_fields)
    for line in inp:
        s = line.strip()
        split = s.split('\t')
        true_gts = split[true_field].split(';')
        pred_gt = split[pred_field]

        true_groups = [set() for _ in range(n)]
        for gt in true_gts:
            for true_group, curr_group in zip(true_groups, get_gt_groups(gt, groups)):
                true_group.add(curr_group)
        for true_group in true_groups:
            s += '\t{}'.format(';'.join(true_group))
        s += '\t' + '\t'.join(get_gt_groups(pred_gt, groups))
        out.write(s + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Input CSV file, to which additional columns will be added.')
    parser.add_argument('-H', '--haplogroups', metavar='FILE', required=True,
        help='CSV file with haplogroups')
    args = parser.parse_args()

    tmp_filename = os.path.join(os.path.dirname(args.input), 'tmp-{}'.format(os.path.basename(args.input)))
    with common.open_stream(args.haplogroups) as inp:
        fields, groups = load_haplogroups(inp)
    with common.open_stream(args.input) as inp, common.open_stream(tmp_filename, 'w') as out:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        add_columns(inp, out, fields, groups)
    os.rename(tmp_filename, args.input)


if __name__ == '__main__':
    main()
