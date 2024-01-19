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


def get_gt_groups(gt, groups, size):
    hap1, hap2 = gt.split(',')
    groups1 = groups.get(hap1) or ('?',) * size
    groups2 = groups.get(hap2) or ('?',) * size
    return tuple(map(','.join, map(sorted, zip(groups1, groups2))))


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

    size = len(group_fields)
    for line in inp:
        s = line.strip()
        split = s.split('\t')
        s += '\t' + '\t'.join(get_gt_groups(split[gt_field], groups, size))
        out.write(s + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar='FILE', required=True,
        help='Input CSV file, to which additional columns will be added.')
    parser.add_argument('-H', '--haplogroups', metavar='FILE', required=True,
        help='CSV file with haplogroups')
    parser.add_argument('-o', '--output', metavar='FILE', required=False,
        help='Output CSV file [default: overwrite input].')
    args = parser.parse_args()

    tmp_filename = common.temporary_filename(args.input)
    with common.open(args.haplogroups) as inp:
        fields, groups = load_haplogroups(inp)
    with common.open(args.input) as inp, common.open(tmp_filename, 'w') as out:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        add_columns(inp, out, fields, groups)
    os.rename(tmp_filename, args.output or args.input)


if __name__ == '__main__':
    main()
