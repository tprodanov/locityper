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

### Preparing input data

#### *k*-mer counts across the reference genome

Several Locityper commands need reference *k*-mer counts,
calculated using [Jellyfish](https://github.com/gmarcais/Jellyfish/).
You can use the following code to obtain them (for example for *k = 25*):
```bash
jellyfish count --canonical --lower-count 2 --out-counter-len 2 --mer-len 25 \
    --threads 8 --size 3G --output counts.jf genome.fa
```

#### Pangenome VCF file

Locityper can automatically extract locus alleles from the pangenome reference,
stored in a VCF file with non-overlapping variants.
Pangenome VCF file can be downloaded from
[here](https://github.com/human-pangenomics/hpp_pangenome_resources#minigraph-cactus) (see *Raw VCF*).
Next, VCF file needs to be transformed such that no variants overlap each other:
```bash
vcfbub -l 0 -i hprc-v1.1-mc-grch38.raw.vcf.gz | bgzip > hprc-v1.1-grch38.vcf.gz
tabix -p vcf hprc-v1.1-grch38.vcf.gz
```

### Creating database with target loci

Database with target loci can be constructed in two ways:
first, locus alleles can be extracted **from a pangenome VCF file**:
```bash
locityper add -d db -r reference.fa -j counts.jf \
    -v pangenome.vcf.gz -l locus_name chr:start-end
```
Please use the same reference genome that was used for Jellyfish *k*-mer counting.
Regions can also be supplied using a four-column BED file using `-L` argument.
Furthermore, `-l/-L` arguments can be used together and can be repeated multiple times:
```bash
locityper add -d db -r reference.fa -j counts.jf \
    -v pangenome.vcf.gz -l locus1 chr:start-end \
    -L loci2.bed -L loci3.bed -l locus4 chr:start-end
```

Alternatively, locus alleles can be **directly provided using a FASTA file**.
```bash
locityper add -d db -r ref.fa -j counts.jf \
    -l locus_name chr:start-end alleles.fa
```
Note, that you still need to provide region coordinates in the reference genome.
This is needed to evaluate off-target *k*-mer counts.
Path to alleles can also be provided in a fifth column of the input BED file (`-L`).

You can freely add more loci to the database using the same commands.

During locus preprocessing, Locityper calculates all pairwise alignment between locus alleles.
Alignment accuracy is controlled by `-a` argument, with `-a 9` producing the most accurate alignments
at slow speed, and `-a 1` quickly producing very inaccurate alignments.
Additionally, you can use `-a 0` to skip alignment step completely.

### Preprocessing WGS dataset

Before locus genotyping can be performed, WGS dataset needs to be preprocessed.
For that, please use
```bash
locityper preproc -i reads1.fastq [reads2.fastq] -j counts.jf \
    -r reference.fa -o preproc_out
```
Input files can have FASTA or FASTQ format, and can be uncompressed or compressed with `gzip`, `bgzip` or `lz4`.

During sample preprocessing Locityper examines read alignments to a long *simple* region in the reference genome
without significant duplications or other structural variantions.
Locityper attempts to automatically identify the genome version, however, if this does not happen,
please use `-b/--bg-region` to provide such region (preferable ≥3 Mb).
By default, the following regions are used: `chr17:72950001-77450000` *(CHM13)*,
`chr17:72062001-76562000` *(GRCh38)* and `chr17:70060001-74560000` *(GRCh37)*.

If you already have read mappings to the whole genome, you can use them via
```bash
locityper preproc -a aligned.bam -d db -r reference.fa -o preproc_out
```
However, this calculation may be less accurate depending on the read mapping.

Additionally, you can estimate WGS characteristics using an already preprocessed file.
> [!WARNING]
> Please use this feature with care if you are certain that the two datasets are similar
> (produced with the same sequencing technology and similar library preparation).
> Preprocessing using existing dataset is significantly faster, but may produce incorrect results
> if datasets have noticeably different characteristics.
```bash
locityper preproc -i reads1.fastq [reads2.fastq] -o preproc_out -~ other_preproc
```

### Genotyping WGS dataset

In order to genotype a dataset, please run
```bash
locityper genotype -i reads1.fastq [reads2.fastq] -p preproc -d db -o analysis
```
