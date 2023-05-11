//! Paths to various directories and files within the database and output directories.

/// k-mer counts are stored in files named `kmers.gz`.
pub const KMERS: &'static str = "kmers.gz";
/// Background region is stored in a FASTA file named `database/bg/bg.fa.gz`.
pub const BG_FASTA: &'static str = "bg.fa.gz";
/// Background information is stored in a directory named `database/bg` or `output/bg`.
pub const BG_DIR: &'static str = "bg";
/// Jellyfish data is stored in `database/jf/`.
pub const JF_DIR: &'static str = "jf";

/// Information about a locus is stored in `database/loci/<locus_name>`.
pub const LOCI_DIR: &'static str = "loci";
/// Reference locus location is stored in `database/loci/<locus_name>/ref.bed`.
pub const LOCUS_BED: &'static str = "ref.bed";
/// Unfiltered haplotypes for the locus are stored in `database/loci/<locus_name>/all_haplotypes.fa.gz`.
pub const LOCUS_FASTA_ALL: &'static str = "all_haplotypes.fa.gz";
/// Filtered haplotypes for the locus are stored in `database/loci/<locus_name>/haplotypes.fa.gz`.
pub const LOCUS_FASTA: &'static str = "haplotypes.fa.gz";
/// Pairwise haplotype alignments are stored in `database/loci/<locus_name>/all_haplotypes.paf.gz`.
pub const LOCUS_PAF: &'static str = "all_haplotypes.paf.gz";
/// Dendrogram is stored in `database/loci/<locus_name>/all_haplotypes.nwk`.
pub const LOCUS_DENDROGRAM: &'static str = "all_haplotypes.nwk";

/// Preprocessing parameters are stored in `output/bg/params.json`.
pub const PREPROC_PARAMS: &'static str = "params.json";
/// Background distributions are stored in `output/bg/distr.gz`.
pub const BG_DISTR: &'static str = "distr.gz";
