//! Paths to various directories and files within the database and output directories.

/// k-mer counts are stored in this file.
pub(super) const KMERS: &'static str = "kmers.bin.lz4";

/// Information about a locus is stored in `database/LOCI_DIR/<locus_name>`.
pub(super) const LOCI_DIR: &'static str = "loci";
/// Reference locus location is stored in `database/loci/<locus_name>/LOCUS_BED`.
pub(super) const LOCUS_BED: &'static str = "ref.bed";
/// Unfiltered haplotypes for the locus are stored in `database/loci/<locus_name>/LOCUS_FASTA_ALL`.
pub(super) const LOCUS_FASTA_ALL: &'static str = "all_haplotypes.fa.gz";
/// Filtered haplotypes for the locus are stored in `database/loci/<locus_name>/LOCUS_FASTA`.
pub(super) const LOCUS_FASTA: &'static str = "haplotypes.fa.gz";
/// Dendrogram is stored in `database/loci/<locus_name>/LOCUS_DENDROGRAM`.
pub(super) const LOCUS_DENDROGRAM: &'static str = "all_haplotypes.nwk";

/// Background distributions are stored in `output/bg/BG_DISTR`.
pub(super) const BG_DISTR: &'static str = "distr.gz";
/// Store alignments in `output/ALNS_DIR/GENOTYPE.bam`.
pub const ALNS_DIR: &'static str = "alns";

/// File, created on the successful completion.
pub(super) const SUCCESS: &'static str = "success";
