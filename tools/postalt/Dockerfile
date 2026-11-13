FROM rust:1.98.1-slim-bookworm AS rsbuild

ARG POSTALT_RS_COMMIT=b48e8b374cd4fe632f8de8932d1f71c531532d6b

RUN apt-get update && apt-get install -y --no-install-recommends build-essential pkg-config liblzma-dev git ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/bsbludwig/postalt_rs /src && cd /src && git checkout "${POSTALT_RS_COMMIT}"
WORKDIR /src
RUN cargo build --release --locked

FROM debian:bookworm-slim AS build

ARG HTSLIB_VERSION=1.22.1
ARG HTSLIB_SHA256=3dfa6eeb71db719907fe3ef7c72cb2ec9965b20b58036547c858c89b58c342f7
ARG SAMTOOLS_SHA256=02aa5cd0ba52e06c2080054e059d7d77a885dfe9717c31cd89dfe7a4047eda0e

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
    && make -j"$(nproc)" && make install \
    && mkdir -p /opt/samtools/share/doc/samtools /opt/samtools/share/doc/htslib \
    && cp /tmp/samtools-${HTSLIB_VERSION}/LICENSE /opt/samtools/share/doc/samtools/LICENSE \
    && cp /tmp/htslib-${HTSLIB_VERSION}/LICENSE /opt/samtools/share/doc/htslib/LICENSE

FROM debian:bookworm-slim

ARG POSTALT_RS_COMMIT=b48e8b374cd4fe632f8de8932d1f71c531532d6b

LABEL org.opencontainers.image.source="https://github.com/break-through-cancer/parabricks-fq2bam-bqsr" \
      org.opencontainers.image.description="postalt_apply: alt-aware postalt_rs (a bwa-postalt.js r985 port, GPL-3.0) on the reads it can change (samtools 1.22.1)" \
      org.opencontainers.image.licenses="GPL-3.0 AND MIT" \
      org.opencontainers.image.version="0.2.6" \
      io.github.break-through-cancer.postalt_rs.source="https://github.com/bsbludwig/postalt_rs/tree/${POSTALT_RS_COMMIT}"

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps python3 zlib1g libbz2-1.0 liblzma5 libdeflate0 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/samtools/ /usr/local/
COPY --from=rsbuild /src/target/release/postalt_rs /usr/local/bin/postalt_rs
COPY --from=rsbuild /src/LICENSE /src/NOTICE.md /usr/share/doc/postalt_rs/
COPY lib/postalt_tools.py /opt/postalt/postalt_tools.py
COPY bin/postalt_apply /usr/local/bin/postalt_apply
ENV POSTALT_TOOLS=/opt/postalt/postalt_tools.py
RUN ldconfig && samtools --version | head -1 && postalt_rs --version && postalt_apply 2>&1 | head -1 || true
