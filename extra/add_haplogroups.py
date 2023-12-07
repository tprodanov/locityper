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
            gt_field = fields.index('genotype')
            s += '\t{}\n'.format('\t'.join(map('{}'.format, group_fields)))
            out.write(s)
            break

    n = len(group_fields)
    for line in inp:
        s = line.strip()
        split = s.split('\t')
        s += '\t' + '\t'.join(get_gt_groups(split[gt_field], groups))
        out.write(s + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Input CSV file, to which additional columns will be added.')
    parser.add_argument('-H', '--haplogroups', metavar='FILE', required=True,
        help='CSV file with haplogroups')
    args = parser.parse_args()

    tmp_filename = common.temporary_filename(args.input)
    with common.open(args.haplogroups) as inp:
        fields, groups = load_haplogroups(inp)
    with common.open(args.input) as inp, common.open(tmp_filename, 'w') as out:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        add_columns(inp, out, fields, groups)
    os.rename(tmp_filename, args.input)


if __name__ == '__main__':
    main()
