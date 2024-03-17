//! Paths to various directories and files within the database and output directories.

/// k-mer counts are stored in this file.
pub(super) const KMERS: &'static str = "kmers.bin.lz4";

/// Information about a locus is stored in `database/LOCI_DIR/<locus_name>`.
pub(super) const LOCI_DIR: &'static str = "loci";
/// Reference locus location is stored in `database/loci/<locus_name>/LOCUS_BED`.
pub(super) const LOCUS_BED: &'static str = "ref.bed";
/// Path to all non-identical locus haplotypes.
pub(super) const LOCUS_FASTA: &'static str = "haplotypes.fa.gz";
/// Path to pairwise distances between haplotypes.
pub(super) const DISTANCES: &'static str = "dist.bin";

/// Background distributions are stored in `output/BG_DISTR`.
pub(super) const BG_DISTR: &'static str = "distr.gz";
/// Save results to `output/RES_JSON`.
pub(super) const RES_JSON: &'static str = "res.json.gz";
/// Store alignments in `output/ALNS_DIR/GENOTYPE.bam`.
pub const ALNS_DIR: &'static str = "alns";

/// File, created on the successful completion.
pub(super) const SUCCESS: &'static str = "success";
