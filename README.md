Locityper
=========

Locityper is a tool, designed to genotype complex loci from short- and long-read whole genome sequencing.

Locityper is created by Timofey Prodanov `timofey.prodanov[at]hhu.de` and Tobias Marschall `tobias.marschall[at]hhu.de`.

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
* [Rust](https://www.rust-lang.org/)
* [Samtools](https://www.htslib.org/)
* [Jellyfish](https://github.com/gmarcais/Jellyfish/)
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

General usage
-------------

### Creating a database

First, we create a database, that would later contain information about complex loci.
```bash
locityper create -d db -r reference.fasta
```
Locityper would attempt to automatically identify the genome version, however, if this does not happen,
please use `-b/--bg-region` to provide a long region (>3Mb) without significant duplications.
This region is later used to evaluate WGS dataset characteristics, such as read depth and error profiles.
By default, the following regions are used: `chr17:72950001-77450000` *(CHM13)*,
`chr17:72062001-76562000` *(GRCh38)* and `chr17:70060001-74560000` *(GRCh37)*.

If you already have a database, you can create a copy of it without any of the complex loci by running command
```bash
rsync -vaP db/{bg,jf} new_db
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
Use can also use `fasta` input files, as well as gzipped input files.
