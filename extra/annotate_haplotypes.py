#!/usr/bin/env python3

from common import open
import weight_subregions
from sortedcontainers import SortedSet


class Feature:
    def __init__(self, s):
        """
        Parse string of format "feature:attribute[:tag=XXX][:backup=XXX][:missing=XXX]"
        """
        split = s.split(':')
        assert len(split) >= 2, f'At least two entries required in "{s}" (feature:attribute).'
        self.feature = split[0]
        self.attributes = [split[1]]
        assert '=' not in self.feature, f'First entry in "{s}" must be the feature name (no = allowed)'
        assert '=' not in self.attributes[0], f'Second entry in "{s}" must be the attribute key (no = allowed)'
        self.tag = self.feature
        self.missing = '<missing>'
        for kw in split[2:]:
            key, val = kw.split('=')
            match key:
                case 'tag':
                    self.tag = val
                case 'backup':
                    self.attributes.append(val)
                case 'missing':
                    self.missing = val


def process_annotation(annotation, sample, gene, trees, features, features_set, hap_annotations):
    for entry in annotation:
        gene_name = entry.attributes['gene_name']
        if gene_name != gene or entry.feature not in features_set:
            continue
        overlap = trees[(sample, entry.contig)].overlap(entry.start, entry.end)
        if not overlap:
            continue
        for ovlp in overlap:
            haplotype = ovlp.data[0]
            for feature in features:
                if entry.feature != feature.feature:
                    continue
                val = feature.missing
                for attr in feature.attributes:
                    curr_val = entry.attributes.get(attr)
                    if curr_val is not None:
                        val = curr_val
                        break
                hap_annotations[haplotype][feature.tag].append(val)


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Annotate haplotypes based on GTF/GFF3 file.')
    parser.add_argument('-a', '--annotation', metavar='FMT', required=True,
        help='Formatting string for GTF/GFF annotation files, where sample name is replaced with {}.')
    parser.add_argument('-g', '--gene', metavar='STR', required=True,
        help='Gene name.')
    parser.add_argument('-l', '--locus', metavar='STR',
        help='Locus name for the output CSV file [default: same as --gene].')
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
        help='Optional: file with discarded haplotypes (DB/loci/LOCUS/discarded_haplotypes.txt).')
    parser.add_argument('-F', '--features', metavar='STR:STR', nargs='+', required=True,
        help='Features of interest of format "feature:attribute[:tag=XXX][:backup=XXX][:missing=XXX]". '
            'Feature name must match the 3rd GFF column, and attribute is a key from the 9th column. '
            'Backup keys can supplied via "backup", they are only checked if the primary attribute is not found. '
            'If no key/backup keys were found, missing value will be used (<missing> by default). '
            'If missing value is empty (missing=) missing values will not be reported. '
            'By default, tag will be equal to the feature, but this can be overridden with tag=XXX.')
    parser.add_argument('-s', '--separator', metavar='STR', default=';',
        help='Separator for outputting repeated annotations [%(default)s].')
    parser.add_argument('-c', '--collapse', action='store_true',
        help='Collapse repeated values in haplotype annotations.')
    parser.add_argument('--write-tag', choices=('yes', 'no', 'auto'), default='auto',
        help='Output tag to the fourht column of the output file? '
            'Options: yes, no, auto (only if there is more than one tag, default).')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    args = parser.parse_args()

    entries = weight_subregions.parse_haplotypes(args)
    discarded_haplotypes = weight_subregions.load_discarded_haps(args.discarded)
    samples, lengths, trees = weight_subregions.process_entries(entries, discarded_haplotypes)

    features = [Feature(s) for s in args.features]
    features_set = { feature.feature for feature in features }
    tags = sorted({ feature.tag for feature in features })
    hap_annotations = { haplotype: { tag: [] for tag in tags } for haplotype in lengths.keys() }

    for sample in samples:
        filename = args.annotation.replace('{}', sample)
        with open(filename) as f:
            gff_entries = weight_subregions.parse_annotation(f)
            process_annotation(gff_entries, sample, args.gene, trees, features, features_set, hap_annotations)

    locus = args.gene if args.locus is None else args.locus
    match args.write_tag:
        case 'yes':
            write_tag = True
        case 'no':
            write_tag = False
        case 'auto':
            write_tag = len(tags) > 1
    with open(args.output, 'wt') as out:
        out.write(f'locus\thaplotype\tannotation{"\ttag" if write_tag else ""}\n')
        for haplotype in sorted(hap_annotations.keys()):
            curr_hap_annotations = hap_annotations[haplotype]
            for tag in tags:
                hap_tag_annots = curr_hap_annotations[tag]
                if args.collapse:
                    hap_tag_annots = SortedSet(hap_tag_annots)
                hap_tag_annots_str = args.separator.join(hap_tag_annots) if hap_tag_annots else '*'
                out.write(f'{locus}\t{haplotype}\t{hap_tag_annots_str}{"\t" + tag if write_tag else ""}\n')


if __name__ == '__main__':
    main()
