Locityper
---------

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

First, we create a database, that would later contain information about complex loci.
```bash
locityper create -d db -r reference.fasta
```


To analyze a locus/loci, you need

```bash
# Create database directory.
locityper create -

```
