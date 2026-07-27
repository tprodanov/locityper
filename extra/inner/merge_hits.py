#!/usr/bin/env python3

import sys
import gzip
from collections import defaultdict
from dataclasses import dataclass
import itertools
import typing
from contextlib import contextmanager


@dataclass
class Mapping:
    chrom: str
    start: int
    end: int
    # + or -
    strand: str
    similarity: float

    def __init__(self, cols: list[str]):
        self.chrom = cols[5]
        self.start = int(cols[7])
        self.end = int(cols[8])
        self.strand = cols[4]
        self.similarity = float(cols[9]) / float(cols[10])

    def __lt__(self, oth):
        return (self.chrom, self.start, self.end).__lt__((oth.chrom, oth.start, oth.end))


@dataclass
class Subregion:
    start: int
    end: int
    # +, -, ? (if both are present), _ (if uncovered)
    strand: str
    similarity: float

    def __len__(self) -> int:
        return self.end - self.start


def _strand_to_int(strand: str) -> int:
    """
    Returns +1 for positive strand, -1 for negative, and 0 for unknown strands.
    """
    match strand:
        case '+':
            return 1
        case '-':
            return -1
        case _:
            return 0


class MergedMapping:
    def __init__(self, target: str, chrom: str, subregions: list[Mapping]):
        self.target: str = target
        self.chrom: str = chrom
        self.start: int = subregions[0].start
        self.end: int = subregions[-1].end
        self.subregions: list[Subregion] = subregions
        self.av_similarity: float = sum(s.similarity * len(s) for s in subregions) / sum(map(len, subregions))

    def __len__(self) -> int:
        return self.end - self.start

    def consensus_strand(self, genome: str):
        first: Subregion = self.subregions[0]
        last: Subregion = self.subregions[-1]
        boundary_strand: int = _strand_to_int(first.strand) * len(first) + _strand_to_int(last.strand) * len(last)
        if first.strand != last.strand or _strand_to_int(first.strand) == 0:
            sys.stderr.write(f'{self.target}\t{genome}\tConflicting strands at the boundaries\n')
        if boundary_strand != 0:
            return '+' if boundary_strand > 0 else '-'
        total_strand = sum(_strand_to_int(subregion.strand) * len(subregion) for subregion in self.subregions)
        if total_strand == 0:
            sys.stderr.write(f'{self.target}\t{genome}\t'
                f'Cannot identify strand for the hit at {self.chrom}:{self.start+1}-{self.end}\n')
        return '+' if total_strand >= 0 else '-'

    def format_strands(self) -> str:
        curr_strand = None
        curr_len = 0
        s = ''
        for subr in self.subregions:
            if curr_len > 0 and curr_strand != subr.strand:
                if s:
                    s += ','
                s += f'{curr_len}{curr_strand}'
            curr_strand = subr.strand
            curr_len += len(subr)
        if s:
            s += ','
        s += f'{curr_len}{curr_strand}'
        return s

    def to_string(self, genome: str, target_len: int) -> str:
        return f'{self.chrom}\t{self.start}\t{self.end}\t{self.consensus_strand(genome)}\t{self.target}' \
            f'\t{len(self) / target_len:.6f}\t{self.av_similarity:.6f}'


def merge(target: str, mappings: list[Mapping]) -> MergedMapping:
    """
    Merge mappings that are already on the same chromosome and within max gap distance.
    """
    # All starts and ends, deduplicated and sorted.
    endpoints: list[int] = sorted(set(itertools.chain((m.start for m in mappings), (m.end for m in mappings))))
    # Quadratic runtime, but should be fine since the number of mappings is probably small.
    subregions: list[Subregion] = []
    for start, end in zip(endpoints, endpoints[1:]):
        strand = '_'
        similarity = 0.0
        for m in mappings:
            if start < m.end and m.start < end:
                strand = m.strand if strand == '_' else '?' if strand != m.strand else strand
                similarity = max(similarity, m.similarity)
        subregions.append(Subregion(start, end, strand, similarity))
    return MergedMapping(target, mappings[0].chrom, subregions)


def process_mappings(
    target: str,
    mappings: list[Mapping],
    max_distance: int,
) -> list[MergedMapping]:
    """
    Process mappings (unsorted and ungrouped) that belong to the same target.
    """
    mappings.sort()
    curr = []
    res = []
    for m in mappings:
        if curr and (curr[-1].chrom != m.chrom or curr[-1].end + max_distance < m.start):
            res.append(merge(target, curr))
            curr.clear()
        curr.append(m)
    if curr:
        res.append(merge(target, curr))
    return res


@contextmanager
def open_file(filename: str) -> typing.IO:
    if filename.endswith(".gz"):
        with gzip.open(filename, 'rt') as f:
            yield f
    else:
        with open(filename) as f:
            yield f


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='FILE',
        help='Alignments of the target to an assembly in PAF format.')
    parser.add_argument('-b', '--targets-bed', metavar='FILE',
        help='Input BED file with target coordinates in the reference genome.')
    parser.add_argument('-g', '--genome', metavar='STR',
        help='Genome/assembly name.')
    parser.add_argument('-d', '--distance', metavar='INT', type=int, required=True,
        help='Maximal distance between hits that will be merged.')
    parser.add_argument('-l', '--min-len', metavar='FLOAT', type=float, default=0.05,
        help='Minimal length as a fraction of the target length [%(default)s].')
    parser.add_argument('-s', '--min-simil', metavar='FLOAT', type=float, default=0.3,
        help='Minimal sequence similarity [%(default)s].')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output BED.gz file. Will contain columns: '
             'chrom, start, end, strand, target name, length fraction, similarity.')
    parser.add_argument('-c', '--cn-out', metavar='FILE', required=True,
        help='Output target copy number (after all filters) to this CSV file.')
    args = parser.parse_args()

    # target -> list of Mapping
    mappings: dict[str, list[Mapping]] = defaultdict(list)
    target_len: dict[str, int] = {}
    with gzip.open(args.input, 'rt') as f:
        for line in f:
            line = line.strip().split('\t')
            target = line[0]
            target_len[target] = int(line[1])
            mappings[target].append(Mapping(line))

    with open_file(args.targets_bed) as f:
        targets = [line.strip().split('\t')[3] for line in f]

    genome = args.genome
    with gzip.open(args.output, 'wt') as out, gzip.open(args.cn_out, 'wt') as cn_out:
        for target in targets:
            target_mappings = mappings.get(target)
            if target_mappings is None:
                cn_out.write(f'{target}\t{genome}\t0\n')
                continue
            merged_mappings = process_mappings(target, target_mappings, args.distance)
            merged_mappings.sort(key=len, reverse=True)
            curr_target_len: int = target_len[target]
            copy_num = 0
            for mapping in merged_mappings:
                if len(mapping) / curr_target_len >= args.min_len and mapping.av_similarity >= args.min_simil:
                    out.write(mapping.to_string(genome, curr_target_len) + '\n')
                    copy_num += 1
            cn_out.write(f'{target}\t{genome}\t{copy_num}\n')


if __name__ == '__main__':
    main()
