FROM rust:1.92-slim-trixie

LABEL org.opencontainers.image.authors="timofey.prodanov@hhu.de"
LABEL version="1.3.5"
LABEL description="Targeted genotyper for complex polymorphic loci from short and long read WGS."

ARG PREFIX=/usr/local
ARG SAMTOOLS_VER=1.19
ARG MM2_VER=2.26
# Strobealign 0.18.0-alpha
ARG STROBEALIGN_SHA1=21c05402bfbcfdd3e2bf59a53677a2b726584464

RUN apt-get update \
    &&  apt-get install --assume-yes --no-install-recommends \
        # To download and extract sources
        curl git bzip2 \
        # To build samtools
        make cmake g++ libncurses-dev zlib1g-dev libbz2-dev liblzma-dev \
        # To build strobealign
        nasm \
        # To use built-in time command, which can measure memory usage
        zsh \
        # Locityper dependency
        jellyfish clang libclang-dev \
    #
    # Building Minimap2
    &&  mkdir /repos \
    &&  cd /repos \
    &&  curl -L https://github.com/lh3/minimap2/releases/download/v${MM2_VER}/minimap2-${MM2_VER}_x64-linux.tar.bz2 | \
        tar --no-same-owner -jxvf - \
    &&  cp minimap2-${MM2_VER}_x64-linux/minimap2    ${PREFIX}/bin \
    &&  cp minimap2-${MM2_VER}_x64-linux/k8          ${PREFIX}/bin \
    &&  cp minimap2-${MM2_VER}_x64-linux/paftools.js ${PREFIX}/bin \
    #
    # Building Samtools
    &&  cd /repos \
    &&  curl -L https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 | \
        tar --no-same-owner -jxvf - \
    &&  cd /repos/samtools-${SAMTOOLS_VER} \
    &&  ./configure \
    &&  make \
    &&  make install prefix=${PREFIX} \
    #
    # Building Strobealign
    &&  cd /repos \
    &&  git clone https://github.com/ksahlin/strobealign \
    &&  cd /repos/strobealign \
    &&  git reset --hard ${STROBEALIGN_SHA1} \
    &&  cargo build --verbose --release --bins \
    &&  cp target/release/strobealign ${PREFIX}/bin \
    #
    # Building Locityper
    &&  cd /repos \
    &&  git clone https://github.com/tprodanov/locityper \
    &&  git clone https://github.com/smarco/WFA2-lib locityper/WFA2 \
    &&  cargo install --no-track --verbose --root "${PREFIX}" --path locityper \
    #
    # Cleaning
    &&  apt-get purge --assume-yes --auto-remove \
        curl git bzip2 make cmake g++ nasm clang \
    &&  cd / \
    &&  rm -rf /repos "${PREFIX}/cargo" "${PREFIX}/rustup"
