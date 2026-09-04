FROM rust:1.92-alpine3.23 AS build

ARG PREFIX=/usr/local
ARG SAMTOOLS_VER=1.19
ARG MM2_VER=2.31
ARG JF_VERSION=2.3.1
# Strobealign 0.18.0-alpha
ARG STROBEALIGN_SHA1=8ba69020e86905b249b51ede94ca6f992c555bc5

RUN apk update \
    &&  apk upgrade --no-cache \
    &&  apk add --no-cache \
        # To download and extract sources
        curl git bzip2 \
        # To build samtools
        make cmake g++ ncurses-dev zlib-dev bzip2-dev xz-dev \
        # To build strobealign
        nasm \
        # Locityper build dependency
        perl clang clang-dev \
    #
    # Copying Minimap2 binaries
    &&  mkdir /repos \
    &&  cd /repos \
    &&  curl -L https://github.com/lh3/minimap2/archive/refs/tags/v${MM2_VER}.tar.gz | \
        tar --no-same-owner -zxvf - \
    &&  cd /repos/minimap2-${MM2_VER} \
    &&  make \
    &&  cp minimap2 $PREFIX/bin \
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
    # Building Jellyfish
    &&  cd /repos \
    &&  curl -L https://github.com/gmarcais/Jellyfish/releases/download/v${JF_VERSION}/jellyfish-${JF_VERSION}.tar.gz | \
        tar --no-same-owner -zxvf - \
    &&  cd /repos/jellyfish-${JF_VERSION} \
    &&  ./configure --prefix=${PREFIX} \
    &&  make \
    &&  make install \
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
    # &&  git clone --single-branch --branch beta https://github.com/tprodanov/locityper \
    &&  git clone --depth 1 https://github.com/tprodanov/locityper \
    &&  git clone --depth 1 https://github.com/smarco/WFA2-lib locityper/WFA2 \
    &&  RUSTFLAGS="-C target-feature=-crt-static" \
        cargo install --no-track --verbose --root "${PREFIX}" --path /repos/locityper \
    #
    # Cleaning
    &&  apk del --purge \
        curl git bzip2 make cmake g++ nasm perl clang \
    &&  cd / \
    &&  rm -rf /repos "${PREFIX}/cargo" "${PREFIX}/rustup"

FROM alpine:3.23

LABEL org.opencontainers.image.authors="timofey.prodanov@hhu.de"
LABEL version="1.7.2"
LABEL description="Targeted genotyper for complex polymorphic loci from short and long read WGS."

ARG PREFIX=/usr/local

RUN apk update \
    &&  apk upgrade --no-cache \
    &&  apk add --no-cache libstdc++-dev libbz2 xz-libs zsh

COPY --from=build ${PREFIX}/lib/libjellyfish*.so* ${PREFIX}/lib
COPY --from=build \
    ${PREFIX}/bin/minimap2 \
    ${PREFIX}/bin/strobealign \
    ${PREFIX}/bin/samtools \
    ${PREFIX}/bin/jellyfish \
    ${PREFIX}/bin/locityper \
    ${PREFIX}/bin/