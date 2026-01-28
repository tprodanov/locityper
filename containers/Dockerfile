FROM rust:1.92-slim-trixie

LABEL org.opencontainers.image.authors="timofey.prodanov@hhu.de"
LABEL version="1.3.3"
LABEL description="Targeted genotyper for complex polymorphic loci from short and long read WGS."

ARG PREFIX=/usr/local
ARG SAMTOOLS_VER=1.19
ARG MM2_VER=2.26
ARG STROBE_VER=0.17.0

RUN apt-get update && \
    apt-get install --assume-yes --no-install-recommends \
        curl clang libclang-dev libncurses-dev zlib1g-dev \
        bzip2 libbz2-dev liblzma-dev make cmake git g++ \
        jellyfish zsh pkgconf libisal-dev

WORKDIR /repos
RUN curl -L https://github.com/lh3/minimap2/releases/download/v${MM2_VER}/minimap2-${MM2_VER}_x64-linux.tar.bz2 | \
    tar --no-same-owner -jxvf - && \
    cp minimap2-${MM2_VER}_x64-linux/minimap2    ${PREFIX}/bin && \
    cp minimap2-${MM2_VER}_x64-linux/k8          ${PREFIX}/bin && \
    cp minimap2-${MM2_VER}_x64-linux/paftools.js ${PREFIX}/bin

WORKDIR /repos
RUN curl -L https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 | \
    tar --no-same-owner -jxvf - && \
    cd /repos/samtools-${SAMTOOLS_VER} && \
    ./configure && \
    make && \
    make install prefix=${PREFIX}

WORKDIR /repos
RUN curl -L https://github.com/ksahlin/strobealign/archive/refs/tags/v${STROBE_VER}.tar.gz | \
    tar --no-same-owner -zxvf - && \
    cd /repos/strobealign-${STROBE_VER} && \
    cmake -B build -DCMAKE_C_FLAGS="-march=native" -DCMAKE_CXX_FLAGS="-march=native" && \
    make -j -C build && \
    cp build/strobealign ${PREFIX}/bin

WORKDIR /repos
RUN git clone https://github.com/tprodanov/locityper && \
    git clone https://github.com/smarco/WFA2-lib locityper/WFA2 && \
    cargo install --no-track --verbose --root "${PREFIX}" --path locityper

WORKDIR /
RUN rm -rf /repos
