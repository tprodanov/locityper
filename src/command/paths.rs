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
/// Haplotypes for the locus are stored in `database/loci/<locus_name>/haplotypes.fa.gz`.
pub const LOCUS_FASTA: &'static str = "haplotypes.fa.gz";
/// Pairwise hapltype alignments are stored in `database/loci/<locus_name>/haplotypes.paf.gz`.
pub const LOCUS_PAF: &'static str = "haplotypes.paf.gz";

/// Sample parameters, estimated from the background regions, are stored in `output/bg/params.gz`.
pub const SAMPLE_PARAMS: &'static str = "params.gz";
