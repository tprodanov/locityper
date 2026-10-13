#!/usr/bin/env python3

from common import open
import weight_subregions


def process_annotation(annot, sample, gene, trees, features, features_set, hap_annotations):
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
                if entry.feature != feature[0]:
                    continue
                tag = feature[0] if len(feature) < 3 else feature[2]
                hap_annotations[haplotype][tag].append(entry.attributes[feature[1]])


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
    parser.add_argument('-f' '--features', metavar='STR:STR', nargs='+',
        help='Features of interest. Each feature (3rd GFF column) must be followed by a colon (:) '
            'and an attribute key (9th column). Optionally, you can set tag name after yet another colon.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    args = parser.parse_args()

    entries = weight_subregions.parse_haplotypes(args)
    discarded_haplotypes = weight_subregions.load_discarded_haps(args.discarded)
    samples, lengths, trees = weight_subregions.process_entries(entries, discarded_haplotypes)

    features = [feature.split(':') for feature in args.features]
    assert all(2 <= len(feature) <= 3 for feature in features)
    features_set = { feature[0] for feature in features }
    tags = sorted({ feature[0] if len(feature) < 3 else feature[2] for feature in features })
    hap_annotations = { haplotype: { tag: [] for tag in tags } for haplotype in lengths.keys() }

    for sample in samples:
        filename = args.annotation.replace('{}', sample)
        with open(filename) as f:
            process_annotation(parse_annotation(f), sample, args.gene, trees, features, features_set, hap_annotations)

    locus = args.gene if args.locus is None else args.locus
    with open(args.output, 'wt') as out:
        out.write('locus\thaplotype\tannotation\ttag\n')
        for haplotype in sorted(hap_annotations.keys()):
            curr_hap_annotations = hap_annotations[haplotype]
            for tag in tags:
                hap_tag_annots = curr_hap_annotations[tag]
                hap_tag_annots_str = ','.join(hap_tag_annots) if hap_tag_annots else '*'
                out.write(f'{locus}\t{haplotype}\t{hap_tag_annots_str}\t{tag}\n')




if __name__ == '__main__':
    main()
