- Filter set of haplotypes based on gaps in read depth.
- Remove `bio` dependency.
- `add`: allow FASTA input.
- WFA: Stop very bad alignments.
- `draw_depth`: handle repeated genotypes.
- Read depth estimation: underestimating variance at low subsampling rates.
- Count edit distance only within the region of interest!
- What to do with read pairs where one of the mates is unmapped?
- Check if alt alleles have Ns. Discard missing haplotypes

Improve `locityper add`:
* Add ability to provide sequences directly,
* Calculate sequence divergencies faster using distances between alleles,
* Do not require to check multiple `extend` parameters,
* Automatically stop extending if there are Ns,
* Reduce window size,
* Write locus summary in a text file.
* Improve Greedy solver.

In advance, check if at least 1/3 of the reads have good alignment scores. How?
Just select best 1/3 and rank haplotypes.

Make read filtering stricter

- Use `itertools` for zip3, and for combinations, combinations with replacement.
- Use separate `Genotype` struct to store tuple of contigs, including parsing and formatting.