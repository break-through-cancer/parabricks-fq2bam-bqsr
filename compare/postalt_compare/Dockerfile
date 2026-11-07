FROM rust:1.98.1-slim-bookworm AS rsbuild
RUN apt-get update && apt-get install -y --no-install-recommends build-essential pkg-config liblzma-dev git ca-certificates && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/bsbludwig/postalt_rs /src && cd /src && git checkout b48e8b374cd4fe632f8de8932d1f71c531532d6b
WORKDIR /src
RUN cargo build --release --locked

FROM debian:bookworm-slim AS build

ARG HTSLIB_VERSION=1.22.1
ARG HTSLIB_SHA256=3dfa6eeb71db719907fe3ef7c72cb2ec9965b20b58036547c858c89b58c342f7
ARG SAMTOOLS_SHA256=02aa5cd0ba52e06c2080054e059d7d77a885dfe9717c31cd89dfe7a4047eda0e
ARG K8_VERSION=1.2
ARG K8_SHA256=a86b160a82f3233a21235d21170f3719a600dbd96bf1ec705a6eb57d770953c9
ARG POSTALT_COMMIT=27dd1da7020151be9358fce67cbe2c10801c348d
ARG POSTALT_SHA256=5648b6830f88f9a04a1682110e9e16ac95dd093d7aa8d070aa9686955a634fc1

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ca-certificates curl bzip2 \
        zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
RUN curl -fsSL -o htslib.tar.bz2 "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
    && echo "${HTSLIB_SHA256}  htslib.tar.bz2" | sha256sum -c - \
    && curl -fsSL -o samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/${HTSLIB_VERSION}/samtools-${HTSLIB_VERSION}.tar.bz2" \
    && echo "${SAMTOOLS_SHA256}  samtools.tar.bz2" | sha256sum -c - \
    && tar -xjf htslib.tar.bz2 && tar -xjf samtools.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure --prefix=/opt/samtools --disable-libcurl --disable-gcs --disable-s3 \
    && make -j"$(nproc)" && make install \
    && cd ../samtools-${HTSLIB_VERSION} \
    && ./configure --prefix=/opt/samtools --with-htslib=/opt/samtools --without-curses \
    && make -j"$(nproc)" && make install

RUN curl -fsSL -o k8.tar.bz2 "https://github.com/attractivechaos/k8/releases/download/v${K8_VERSION}/k8-${K8_VERSION}.tar.bz2" \
    && echo "${K8_SHA256}  k8.tar.bz2" | sha256sum -c - \
    && tar -xjf k8.tar.bz2 \
    && install -m 0755 k8-${K8_VERSION}/k8-x86_64-Linux /opt/k8 \
    && curl -fsSL -o /opt/bwa-postalt.js "https://raw.githubusercontent.com/lh3/bwa/${POSTALT_COMMIT}/bwakit/bwa-postalt.js" \
    && echo "${POSTALT_SHA256}  /opt/bwa-postalt.js" | sha256sum -c -

FROM debian:bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/break-through-cancer/parabricks-fq2bam-bqsr" \
      org.opencontainers.image.description="Comparison of bwa-postalt.js (k8) and postalt_rs; contains GPL-3.0 postalt_rs" \
      org.opencontainers.image.version="0.1.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps python3 time zlib1g libbz2-1.0 liblzma5 libdeflate0 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/samtools/ /usr/local/
COPY --from=build /opt/k8 /usr/local/bin/k8
COPY --from=build /opt/bwa-postalt.js /opt/postalt/bwa-postalt.js
COPY --from=rsbuild /src/target/release/postalt_rs /usr/local/bin/postalt_rs
COPY postalt_tools.py /opt/postalt/postalt_tools.py
COPY compare_engines.sh /usr/local/bin/compare_engines.sh
COPY compare_records.py /opt/postalt/compare_records.py
RUN chmod +x /usr/local/bin/compare_engines.sh && ldconfig && samtools --version | head -1 && k8 -v | head -1 && postalt_rs --version
