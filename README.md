Locityper
=========

Locityper is designed to **genotype complex loci** from *short-* and *long-read* whole genome sequencing.

Locityper is created by
[Timofey Prodanov](https://marschall-lab.github.io/people/tprodanov/) `timofey.prodanov[at]hhu.de` and
[Tobias Marschall](https://marschall-lab.github.io/people/tmarschall/) `tobias.marschall[at]hhu.de`.

Table of contents
-----------------
* [Citing Locityper](#citing-locityper)
* [Installation](#installation)
* [General usage](#general-usage)
* [Output files](#output-files)
* [See also](#see-also)

Citing Locityper
----------------

Publication in progress, please check later.

Installation
------------

For manual installation, you need the following software:
* [Rust](https://www.rust-lang.org/),
* [GCC ≥ 4.9.0](https://gcc.gnu.org/),
* [Samtools](https://www.htslib.org/),
* [Jellyfish](https://github.com/gmarcais/Jellyfish/),
* [Minimap2](https://github.com/lh3/minimap2) if you plan to analyze long-read sequencing,
* [Strobealign](https://github.com/ksahlin/strobealign) if you plan to analyze short-read sequencing.

Next, please follow the next steps:
```bash
git clone https://github.com/tprodanov/locityper
cd locityper
git clone https://github.com/smarco/WFA2-lib WFA2
cargo build --release
cargo install
```

### Compiling without WFA and GCC

It is possible to reduce the number of dependencies by sacrificing some functionality.
In particular, one can compile `locityper` without `WFA` library and appropriate `GCC` compiler with
```bash
git clone https://github.com/tprodanov/locityper
cd locityper
cargo build --release --no-default-features
cargo install
```
However, in such configuration it is impossible to align locus haplotypes and discard similar alleles,
only identical alleles will be removed.

General usage
-------------

Note that many of all commands below allow to specify the number of threads (`-@ N`, 8 by default).

### Creating a database

First, we create a database, that would later contain information about complex loci.
```bash
locityper create -d db -r reference.fasta
```
During database creation, Locityper collects *k*-mer counts across the reference genome using `jellyfish`.

If you already have a database, you can create a copy of it without any of the complex loci by running
```bash
rsync -vaP db/jf new_db
```

### Adding loci to the database

There are two ways to add loci to the database. First, locus alleles can be directly provided using a FASTA file.
```bash
locityper add -d db -s alleles.fasta=name # where `name` is the name of the locus.
```

Alternatively, you can use pangenome VCF file:
```bash
locityper add -d db -r reference.fasta -v pangenome.vcf.gz -l chr:start-end=name
```
**Pangenome VCF file** can be downloaded from
[here](https://github.com/human-pangenomics/hpp_pangenome_resources#minigraph-cactus) (see *Raw VCF*).
Next, VCF file needs to be transformed such that no variants overlap each other:
```bash
vcfbub -l 0 -i hprc-v1.1-mc-grch38.raw.vcf.gz | bgzip > hprc-v1.1-grch38.vcf.gz
tabix -p vcf hprc-v1.1-grch38.vcf.gz
```

### Preprocessing WGS dataset

Before locus genotyping can be performed, WGS dataset needs to be preprocessed.
For that, please use
```bash
locityper preproc -i reads1.fastq [reads2.fastq] -d db -r reference.fasta -o analysis
```
Use can also provide input files in FASTA format.
Locityper accepts input files compressed with `gzip`, `bgzip` or `lz4`.

During sample preprocessing Locityper examines read alignments to a long *simple* region in the reference genome
without significant duplications or other structural variantions.
Locityper attempts to automatically identify the genome version, however, if this does not happen,
please use `-b/--bg-region` to provide such region (preferable ≥3 Mb).
By default, the following regions are used: `chr17:72950001-77450000` *(CHM13)*,
`chr17:72062001-76562000` *(GRCh38)* and `chr17:70060001-74560000` *(GRCh37)*.

There are two other ways to preprocess WGS dataset:
If you already have read mappings to the whole genome, you can use them via
```bash
locityper preproc -a aligned.bam -d db -r reference.fasta -o analysis
```
However, this calculation may be less accurate than from the reads depending on the read mapping.

Additionally, you can estimate WGS characteristics using an already preprocessed file.
> [!WARNING]
> Please use this feature with care if you are certain that the two datasets are similar
> (produced with the same sequencing technology and similar library preparation).
> Preprocessing using existing dataset is significantly faster, but may produce incorrect results
> if datasets have noticeably different characteristics.
```bash
locityper preproc -i reads1.fastq [reads2.fastq] -o analysis -~ other_analysis
```

Finally, if you want to analyze a new set of loci, and do not want to keep the same output folder, you
can copy preprocessed data to a new location with
```bash
rsync -vaP analysis1/bg analysis2
```

### Genotyping WGS dataset

In order to genotype a dataset, please run
```bash
locityper genotype -i reads1.fastq [reads2.fastq] -d db -o analysis
```
