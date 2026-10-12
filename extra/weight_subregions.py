#!/usr/bin/env python3

from collections import defaultdict, namedtuple
import re
import os
from sortedcontainers import SortedList
from intervaltree import IntervalTree
from common import open


_SAMPLE_SEP = re.compile('[_./-]')

class Entry:
    def __init__(self, haplotype, contig, start, end, strand):
        self.haplotype = haplotype
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.sample = _SAMPLE_SEP.split(self.haplotype, 1)[0]


class Feature:
    def __init__(self, pattern, weight):
        assert weight >= 0.0
        self.pattern = re.compile(f'^({pattern})$') if pattern is not None else None
        self.weight = weight

    @classmethod
    def from_str(cls, s):
        pattern, weight = s.rsplit(':', 1)
        return cls(pattern, float(weight))


Endpoint = namedtuple('Endpoint', 'pos is_start feature')


def create_entries(args):
    if args.fasta and args.direct:
        with open(args.fasta) as f:
            haplotype = None
            curr_len = 0
            for line in f:
                if line.startswith('>'):
                    if haplotype is not None:
                        yield Entry(haplotype, haplotype, 0, curr_len, '+')
                    haplotype = line[1:].strip().split()[0]
                    curr_len = 0
                else:
                    curr_len += len(line.strip())
            if haplotype is not None:
                yield Entry(haplotype, haplotype, 0, curr_len, '+')

    elif args.fasta:
        with open(args.fasta) as f:
            for line in f:
                if not line.startswith('>'):
                    continue
                split = line[1:].strip().split()
                contig, coord, strand = split[1].split(':')
                start, end = map(int, coord.split('-'))
                yield Entry(split[0], contig, start-1, end, strand)

    elif args.regions:
        with open(args.regions) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                split = line.strip().split('\t')
                yield Entry(split[3], split[0], int(split[1]), int(split[2]), split[5])
    else:
        raise ValueError('Either -f/--fasta or -r/--regions must be specified')


def load_discarded_haps(filename):
    discarded_haps = set()
    if filename is not None and os.path.exists(filename):
        with open(filename) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                split = line.split(' ')
                assert split[1] == '=' or split[1] == '~'
                discarded_haps.update(split[2:])
    return discarded_haps


def process_entries(entries, discarded_haps):
    samples = set()
    lengths = {}
    trees = defaultdict(IntervalTree)
    for entry in entries:
        if entry.haplotype in discarded_haps:
            continue
        assert entry.haplotype not in lengths
        assert not trees[(entry.sample, entry.contig)].overlaps(entry.start, entry.end)
        assert entry.strand == '+' or entry.strand == '-'
        samples.add(entry.sample)
        lengths[entry.haplotype] = entry.end - entry.start
        trees[(entry.sample, entry.contig)].addi(entry.start, entry.end, (entry.haplotype, entry.strand))
    return samples, lengths, trees


# Adapted from https://github.com/HTGenomeAnalysisUnit/gtf-parser
class GffEntry:
    def is_gff3(attr_string: str) -> str:
        """Detect whether an attribute string is GTF or GFF3 format."""
        return "=" in attr_string and not re.search(r'"\s*;', attr_string)

    @staticmethod
    def _parse_gtf_attributes(raw: str) -> dict[str, str]:
        """Parse GTF-style attributes: key "value"; key2 "value2";"""
        attrs: dict[str, str] = {}
        for match in re.finditer(r'(\w+)\s+"([^"]*)"', raw):
            attrs[match.group(1)] = match.group(2)
        return attrs

    @staticmethod
    def _parse_gff3_attributes(raw: str) -> dict[str, str]:
        """Parse GFF3-style attributes: key=value;key2=value2"""
        attrs: dict[str, str] = {}
        for part in raw.split(";"):
            part = part.strip()
            if "=" in part:
                key, _, value = part.partition("=")
                attrs[key.strip()] = value.strip()
        return attrs

    def __init__(self, split, is_gff3):
        self.contig = split[0]
        self.source = split[1]
        self.feature = split[2]
        self.start = int(split[3])
        self.end = int(split[4])
        self.score = split[5]
        self.strand = split[6]
        self.frame = split[7]
        self.attributes = GffEntry._parse_gff3_attributes(split[8]) if is_gff3 \
            else GffEntry._parse_gtf_attributes(split[8])

    def __len__(self):
        return self.end - self.start


def parse_annotation(f):
    is_gff3 = None
    for line in f:
        if line.startswith('#'):
            pass
        split = line.strip().split('\t')
        if is_gff3 is None:
            is_gff3 = GffEntry.is_gff3(split[8])
        yield GffEntry(split, is_gff3)


def process_annotation(annotation, sample, gene, trees, features, endpoints):
    for entry in annotation:
        gene_name = entry.attributes['gene_name']
        if gene_name != gene:
            continue
        overlap = trees[(sample, entry.contig)].overlap(entry.start, entry.end)
        if not overlap:
            continue
        feature_ix = next(i for i, feature in enumerate(features) if feature.pattern is None or feature.pattern.match(entry.feature))
        for ovlp in overlap:
            start = max(entry.start - ovlp.begin, 0)
            end = min(entry.end, ovlp.end) - ovlp.begin
            haplotype, strand = ovlp.data
            if strand == '-':
                length = ovlp.end - ovlp.begin
                start, end = length - end, length - start
            else:
                assert strand == '+'
            endpoints[haplotype].append(Endpoint(start, True, feature_ix))
            endpoints[haplotype].append(Endpoint(end, False, feature_ix))


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Assign weights to haplotype subregions based '
        'on annotation in the GTF format. Haplotype = sample + suffix, all haplotypes are unique and only contain '
        'sequence around the region of interest. Contig names can repeat, but are unique within a given sample.')
    parser.add_argument('-a', '--annotation', metavar='FMT', required=True,
        help='Formatting string for GTF/GFF annotation files, where sample name is replaced with {}.')
    parser.add_argument('-g', '--gene', metavar='STR', required=True,
        help='Gene name.')
    parser.add_argument('-f', '--fasta', metavar='FILE',
        help='FASTA file with extracted haplotypes (mutually exclusive with -r). '
            'For each each haplotype line should have format '
            '">haplotype_name contig:start-end:strand ..." unless --direct. '
            'Sample name = haplotype name until dot/dash/underscore.')
    parser.add_argument('--direct', action='store_true',
        help='Annotations were created directly on the local haplotypes: '
            'GTF contigs = haplotype names and offset is 0.')
    parser.add_argument('-r', '--regions', metavar='FILE',
        help='BED file with columns contig, start, end, haplotype_name, _, strand (mutually exclusive with -f).')
    parser.add_argument('-d', '--discarded', metavar='FILE',
        help='Optional: file with discarded haplotypes (DB/loci/LOCUS/discarded_haplotypes.txt)')

    default_features = ['CDS|start_codon|end_codon:1', 'gene:0.5']
    default_features_str = ' '.join(map('"{}"'.format, default_features))
    parser.add_argument('-F', '--features', nargs='+', metavar='PATTERN:NUM', default=default_features,
        help='Feature and and corresponding weight, separated with ":". '
            'If multiple feature patterns match, first matched applies. '
            f'Default: -F {default_features_str}.')
    parser.add_argument('-b', '--background', metavar='NUM', default=0.01,
        help='Weight of the background sequence that does not overlap any features [%(default)s].')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output BED file.')
    args = parser.parse_args()

    entries = create_entries(args)
    discarded_haplotypes = load_discarded_haps(args.discarded)
    samples, lengths, trees = process_entries(entries, discarded_haplotypes)

    features = list(map(Feature.from_str, args.features)) + [Feature(None, args.background)]
    bg_ix = len(features) - 1
    endpoints = { haplotype: [Endpoint(0, True, bg_ix), Endpoint(length, False, bg_ix)]
        for haplotype, length in lengths.items() }
    for sample in samples:
        filename = args.annotation.replace('{}', sample)
        with open(filename) as f:
            process_annotation(parse_annotation(f), sample, args.gene, trees, features, endpoints)

    enabled = SortedList()
    with open(args.output, 'wt') as out:
        for haplotype in sorted(endpoints.keys()):
            hap_endpoints = endpoints[haplotype]
            hap_endpoints.sort()
            last_pos = 0
            enabled.clear()
            for pos, is_start, feature_ix in hap_endpoints:
                if last_pos < pos:
                    # Smallest index of an enabled feature
                    out_feature_ix = enabled[0]
                    out.write(f'{haplotype}\t{last_pos}\t{pos}\t{features[out_feature_ix].weight:.8g}\n')
                if is_start:
                    enabled.add(feature_ix)
                else:
                    enabled.remove(feature_ix)
                last_pos = pos
            assert not enabled and last_pos == lengths[haplotype]


if __name__ == '__main__':
    main()
