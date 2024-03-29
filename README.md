<h1 align="center">Locityper</h1>

<div align="center">

<a href="">[![Last updated](https://img.shields.io/github/last-commit/tprodanov/locityper/main.svg?label=Last%20updated&color=blue&style=flat-square)](https://github.com/tprodanov/locityper/commits/main/)</a>
<a href="">[![Last release](https://img.shields.io/github/v/tag/tprodanov/locityper.svg?label=GitHub&color=blueviolet&style=flat-square)](https://github.com/tprodanov/locityper/releases)</a>
<a href="">[![Bioconda](https://img.shields.io/conda/v/bioconda/locityper.svg?label=Bioconda&color=blue&style=flat-square)](https://anaconda.org/bioconda/locityper)</a>

</div>

Locityper is targeted genotyper for polymorphic genes based on short- and long-read whole genome sequencing.

Locityper is created by
[Timofey Prodanov](https://marschall-lab.github.io/people/tprodanov/) `timofey.prodanov[at]hhu.de` and
[Tobias Marschall](https://marschall-lab.github.io/people/tmarschall/) `tobias.marschall[at]hhu.de`
at the [Marschall Lab](https://marschall-lab.github.io/), [Heinrich Heine University Düsseldorf](hhu.de).

# Table of contents
* [Citing Locityper](#citing-locityper)
* [Installation](#installation)
* [General usage](#general-usage)
* [Output files](#output-files)
* [See also](#see-also)

# Citing Locityper

Publication in progress, please check later.

<br>

# Installation

## Conda

You can use `conda`/`mamba` to easily install `locityper`:
```bash
conda install -c bioconda locityper
# OR
mamba install -c bioconda locityper
```

## Singularity

It is possible to compile `locityper` using
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) container platform.
For that, simply run
```bash
singularity build --fakeroot locityper.sif containers/locityper.def
# OR
sudo singularity build locityper.sif containers/locityper.def
```
and then execute with
```bash
./locityper.sif ...
# Execute dependencies with
./locityper.sif jellyfish/minimap2/strobealign ...
# Alternatively, use
singularity exec locityper.sif ...
```

## Docker

Alternatively, you can compile `locityper` via [Docker](https://www.docker.com/):
```bash
sudo docker build -t locityper containers
# OR
sudo docker buildx build -t locityper containers
```
Then, you can execute with
```bash
sudo docker run locityper locityper ...
```

## Manual installation

For manual installation, you need the following software:
* [Rust](https://www.rust-lang.org/),
* [GCC ≥ 4.9.0](https://gcc.gnu.org/),
* [Samtools](https://www.htslib.org/),
* [Jellyfish](https://github.com/gmarcais/Jellyfish/),
* [Minimap2](https://github.com/lh3/minimap2) if you plan to analyze long-read sequencing,
* [Strobealign](https://github.com/ksahlin/strobealign) if you plan to analyze short-read sequencing.

Next, please follow the next steps:
<pre><code>git clone <ins>https://github.com/tprodanov/locityper</ins>
cd locityper
git clone <ins>https://github.com/smarco/WFA2-lib</ins> WFA2
cargo build --release
cargo install
</code></pre>

### Compiling without WFA and GCC

It is possible to reduce the number of dependencies by sacrificing some functionality.
In particular, one can compile `locityper` without the `WFA` library and the `GCC` compiler with

<pre><code>git clone <ins>https://github.com/tprodanov/locityper</ins>
cd locityper
cargo build --release <b>--no-default-features</b>
cargo install
</code></pre>
However, in such configuration Locityper will not calculate pairwise haplotype alignments
in order to discard similar alleles, only identical alleles will be removed.

<br>

# General usage

Note that many of all commands below allow to specify the number of threads using `-@` argument.

## Prerequisites

### *k*-mer counts across the reference genome

Locityper utilizes *k*-mer counts across the refernece genome,
calculated using [Jellyfish](https://github.com/gmarcais/Jellyfish/).
You can use the following code to obtain them (for example for *k* = 25):
```bash
jellyfish count --canonical --lower-count 2 --out-counter-len 2 --mer-len 25 \
    --threads 8 --size 3G --output counts.jf genome.fa
```

### Pangenome VCF file

Locityper can automatically extract locus alleles from the pangenome reference,
stored in a VCF file with non-overlapping variants.
Pangenome VCF file can be found
[here](https://github.com/human-pangenomics/hpp_pangenome_resources#minigraph-cactus) (see *Raw VCF*).
Next, VCF file needs to be transformed such that no variants overlap each other
using [vcfbub](https://github.com/pangenome/vcfbub):
```bash
vcfbub -l 0 -i hprc-v1.1-mc-grch38.raw.vcf.gz | bgzip > hprc-v1.1-grch38.vcf.gz
tabix -p vcf hprc-v1.1-grch38.vcf.gz
```

## Preparing a database with target loci

Database with target loci can be constructed using `locityper add` command.
First of all, locus alleles can be extracted **from a pangenome VCF file**:

<pre><code>locityper add< -d db -r reference.fa -j counts.jf \
    -v <b>pangenome.vcf.gz</b> -l <b>locus_name chr:start-end</b>
</code></pre>

Please use the same reference genome as for Jellyfish *k*-mer counting.
Regions can also be supplied with a four-column BED file using `-L` argument.
Furthermore, `-l/-L` arguments can be used together and can be repeated multiple times:

<pre><code>locityper add -d db -r reference.fa -j counts.jf \
    -v pangenome.vcf.gz <b>-l locus1 chr:start-end \
    -L loci2.bed -L loci3.bed -l locus4 chr:start-end</b>
</code></pre>

Alternatively, locus alleles can be **directly provided using a FASTA file**.

<pre><code>locityper add -d db -r ref.fa -j counts.jf \
    -l locus_name chr:start-end <b>alleles.fa</b>
</code></pre>

Note, that you still need to provide region coordinates in the reference genome.
This is needed to evaluate off-target *k*-mer counts.
Path to alleles can also be provided in a fifth column of the input BED file (`-L`).

If you don't know exact region coordinates, you can align locus alleles to the reference genome with
```bash
minimap2 -cx asm20 genome.fa alleles.fa | \
    awk '{ printf("%s:%d-%d\n",$6,$8+1,$9) }' | \
    sort | uniq -c | sort -k1,1nr | head
```
This will produce a list of possible reference coordinates ranked by the number of alleles.
You can then select the most common coordinates / largest region at your discretion.

You can freely add more loci to an existing database using the same commands.

### Aligning locus alleles

During locus preprocessing, Locityper calculates all pairwise alignment between locus alleles,
which can be found in the output directory `db/loci/<locus>/all_haplotypes.bam`.
Alignment accuracy is controlled by `-a` argument, with `-a 9` producing the most accurate alignments
at slow speed, and `-a 1` quickly producing very inaccurate alignments.
Additionally, you can use `-a 0` to skip alignment step completely.

## Preprocessing WGS dataset

Before locus genotyping can be performed, WGS dataset needs to be preprocessed.
For that, please use
<pre><code><b>locityper preproc</b> -i reads1.fastq [reads2.fastq] \
    -j counts.jf -r reference.fa -o preproc_out
</code></pre>

Input files can have FASTA or FASTQ format, and can be uncompressed or compressed with `gzip`, `bgzip` or `lz4`.

During sample preprocessing Locityper examines read alignments to a long *simple* region in the reference genome
without significant duplications or other structural variantions.
Locityper attempts to automatically identify the genome version, however, if this does not happen,
please use `-b/--bg-region` to provide such region (preferable ≥3 Mb).
By default, the following regions are used: `chr17:72950001-77450000` *(CHM13)*,
`chr17:72062001-76562000` *(GRCh38)* and `chr17:70060001-74560000` *(GRCh37)*.

If you already have read mappings to the whole genome, you can use them via
<pre><code>locityper preproc <b>-a aligned.bam</b> \
    -j counts.jf -r reference.fa -o preproc_out
</code></pre>

In some cases, reads are stored in an unmapped BAM/CRAM format.
In cases like that you can use arguments `-a reads.bam --no-index`.
However, this option is only allowed for single-end or interleaved paired-end (`-^`) data.

Additionally, you can estimate WGS characteristics using an already preprocessed file.
> [!WARNING]
> Please use this feature with care if you are certain that the two datasets are similar
> (produced with the same sequencing technology and similar library preparation).
> Preprocessing using existing dataset is significantly faster, but may produce incorrect results
> if datasets have noticeably different characteristics.
<pre><code>locityper preproc -i/-a input <b>-~ other</b> \
    -j counts.jf -r reference.fa -o preproc_out
</code></pre>

## Genotyping WGS dataset

In order to genotype a dataset, please run
<pre><code><b>locityper genotype</b> -i reads1.fastq [reads2.fastq] -p preproc -d db -o analysis
</code></pre>

Simlarly to the preprocessing step, you can analyse existing read mappings using
`-a aligned.bam [--no-index]`.

Additionally, you can specify only some loci from the database using `--subset-loci` argument
and specify genotype priors with `--priors`.

## Input files for preprocessing/genotyping WGS dataset

Locityper `preproc` and `genotype` commands allow multiple input files for a single dataset.
You can specify them by repeating `-i/-a` arguments:
```bash
locityper preproc -i readsA.fq.gz -i readsB.fq.gz ...
locityper preproc -i readsA1.fq.gz readsA2.fq.gz \
    -i readsB1.fq.gz readsB2.fq.gz ...
locityper preproc -a readsA.bam -a readsB.bam --no-index ...
```
> [!WARNING]
> All of the files must correspond to the same sample, same sequencing technology and should have roughly the same
> characteristics (read length, error rates).
> Additionally, please make sure to use the same input files for genotyping as for preprocessing.

Alternatively, you can specify an input list of files with `-I list.txt` where each line is
`<flag> <file> [<file2>]`.
Specifically, lines can be
- Two paired-end files: `p reads1.fq.gz reads2.fq.gz` or `p reads*.fq.gz`,
- Interleaved paired-end file: `pi reads.fq.gz` (same as `-i reads.fq.gz --interleaved`),
- Single-end FASTA/Q file: `s reads.fq.gz`,
- Alignment and mapped BAM/CRAM file: `a alns.bam`,
- Unmapped BAM/CRAM file: `u reads.bam` (same as `-a reads.bam --no-index`).
- Unmapped interleaved BAM/CRAM file: `ui reads.bam` (same as `-a reads.bam --no-index --interleaved`).

Multiple lines are allowed, but all must have the same flag.
