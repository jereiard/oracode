# Ion Torrent toolchain (tmap/tvc) provider. Override to point at your own
# Torrent Suite image.
ARG TS_IMAGE=jereiard/tsuite:5.18.1
FROM ${TS_IMAGE} AS tsuite

# ---------------------------------------------------------------------------
# Base: Ubuntu 20.04 + Miniforge (Python 3.10) + samtools/bcftools + BAMSurgeon
# + Picard, all from public sources. Building this here (rather than pulling a
# committed image) keeps the repository buildable by anyone.
# ---------------------------------------------------------------------------
FROM ubuntu:20.04 AS base

# Proxy support for builds behind a corporate proxy. Build with
#   --build-arg HTTP_PROXY=... --build-arg HTTPS_PROXY=...
# (dist.sh forwards these automatically). They are build-only; the final image
# does not rely on them.
ARG HTTP_PROXY
ARG HTTPS_PROXY
ARG http_proxy
ARG https_proxy
ENV HTTP_PROXY=${HTTP_PROXY} \
    HTTPS_PROXY=${HTTPS_PROXY} \
    http_proxy=${http_proxy} \
    https_proxy=${https_proxy}

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    CONDA_DIR=/opt/conda \
    PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       wget bzip2 ca-certificates tini curl unzip patch \
       build-essential \
    && rm -rf /var/lib/apt/lists/*

# Miniforge, pinned.
ARG MINIFORGE_VERSION=24.3.0-0
RUN wget --no-hsts -q \
      "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh" \
      -O /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p "${CONDA_DIR}" \
    && rm /tmp/miniforge.sh \
    && conda clean -afy

# Bioinformatics tools plus pip-installed pysam, from the pinned list captured
# from a working build.
COPY conda-explicit.txt /tmp/conda-explicit.txt
RUN conda install -y --file /tmp/conda-explicit.txt \
    && conda clean -afy \
    && rm /tmp/conda-explicit.txt \
    && python -m pip install --no-cache-dir pysam==0.22.1

# BAMSurgeon at the pinned commit, and Picard 1.131.
# Fetched over HTTPS with curl: some corporate proxies reset plain git protocol
# transfers, while curl (already used for IGV) works.
ARG BAMSURGEON_REPO=jereiard/bamsurgeon
ARG BAMSURGEON_COMMIT=50b14964eb2205fb5e3838d366174103024be050
ARG PICARD_VERSION=1.131
RUN mkdir -p /opt/reflower/deps \
    && curl -kfsSL --retry 3 \
         "https://codeload.github.com/${BAMSURGEON_REPO}/tar.gz/${BAMSURGEON_COMMIT}" \
         -o /tmp/bamsurgeon.tar.gz \
    && mkdir -p /tmp/bamsurgeon \
    && tar -xzf /tmp/bamsurgeon.tar.gz -C /tmp/bamsurgeon --strip-components=1 \
    && mv /tmp/bamsurgeon /opt/reflower/deps/bamsurgeon \
    && rm -rf /tmp/bamsurgeon /tmp/bamsurgeon.tar.gz \
    && curl -kfsSL --retry 3 \
         "https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip" \
         -o /tmp/picard.zip \
    && unzip -q /tmp/picard.zip -d /tmp/picard \
    && mv "/tmp/picard/picard-tools-${PICARD_VERSION}" /opt/reflower/deps/picard-tools-${PICARD_VERSION} \
    && rm -rf /tmp/picard /tmp/picard.zip

RUN apt-get update \
    && apt-get install -y --no-install-recommends redis-server \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Oracode image
# ---------------------------------------------------------------------------
FROM base

ARG IGV_VERSION=2.19.8
ARG IGV_SHA256=754ab357c9dfe207e663c5363741550dd9ed71b3caaaba6f4f6d0edfe66d8791

# The supplied Illumina panel BAM contains paired reads whose mate record is
# absent. Picard SamToFastq otherwise discards those selected spike-in reads,
# which systematically lowers the post-realignment VAF. Preserve and remap the
# orphan FASTQ alongside complete pairs.
COPY patches/bamsurgeon-preserve-orphans.patch /tmp/bamsurgeon-preserve-orphans.patch
RUN cd /opt/reflower/deps/bamsurgeon \
    && patch -p1 < /tmp/bamsurgeon-preserve-orphans.patch \
    && rm /tmp/bamsurgeon-preserve-orphans.patch

# The TS image is built on Ubuntu 20.04 in tsbuild/ so its binaries are ABI
# compatible with this base image. TVC additionally needs these runtime libs.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       libarmadillo9 \
       libarpack2 \
       libsuperlu5 \
       xvfb \
       xauth \
       fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

# The Broad download endpoint can present a proxy-injected certificate during
# corporate builds. Authenticity is enforced by the pinned archive SHA-256.
RUN curl -kfsSL --retry 3 \
      "https://data.broadinstitute.org/igv/projects/downloads/2.19/IGV_Linux_${IGV_VERSION}_WithJava.zip" \
      -o /tmp/igv.zip \
    && echo "${IGV_SHA256}  /tmp/igv.zip" | sha256sum -c - \
    && unzip -q /tmp/igv.zip -d /opt \
    && mv "/opt/IGV_Linux_${IGV_VERSION}" /opt/igv \
    && ln -s /opt/igv/igv.sh /usr/local/bin/igv.sh \
    && rm /tmp/igv.zip
# Keep IGV off the network for its Google OAuth lookup (air-gapped hosts get a
# proxy 407 / timeout that raises a modal dialog and stalls headless batch runs).
# run_igv passes --igvDirectory /opt/igv/offline, whose oauth-config.json makes
# OAuthUtils register a provider without contacting igv.org.
COPY oauth-config.json /opt/igv/offline/oauth-config.json
COPY igv-headless.sh /usr/local/bin/igv-headless

COPY --from=tsuite /usr/local/bin/tmap /usr/local/bin/tmap
COPY --from=tsuite /usr/local/bin/tvc /usr/local/bin/tvc
COPY --from=tsuite /usr/local/bin/tvcutils /usr/local/bin/tvcutils
COPY --from=tsuite /usr/local/bin/tvcassembly /usr/local/bin/tvcassembly
COPY --from=tsuite /usr/local/share/variantCaller/ /usr/local/share/variantCaller/

ENV TVC_PARAMETER_DIR=/usr/local/share/variantCaller/parameter_sets
RUN python -m pip install --no-cache-dir textual==8.2.7

# Pipeline code. The base image ships an older copy under /opt/reflower; the
# build always overlays the sources from this repository so the image matches
# the tree it was built from.
COPY bespoke.py curses_utils.py educe_variant.py extract_roi.py igvbatch.py \
     logfile_manager.py oracode.py proficiency.py reflow.py remede.py \
     reporting.py textual_ui.py /opt/reflower/
COPY entrypoint.sh /opt/reflower/entrypoint.sh
RUN chmod +x /opt/reflower/entrypoint.sh
RUN mkdir -p /root/miniconda3/bin
RUN ln -sfn /usr/local/bin/tmap /root/miniconda3/bin/tmap-ion
RUN chmod 777 /root/miniconda3/bin/tmap-ion
RUN tmap --version 2>&1 | grep -q '5.18.6' \
    && tvc --version 2>&1 | grep -q '5.18-6' \
    && test -f "$TVC_PARAMETER_DIR/parameter_sets.json" \
    && test -x /opt/igv/igv.sh \
    && test -x /usr/local/bin/igv-headless \
    && /opt/igv/jdk-21/bin/java -version 2>&1 | grep -q '21\.'
COPY patches/bamsurgeon-mark-orphans.patch /tmp/bamsurgeon-mark-orphans.patch
COPY patches/bamsurgeon-balance-single-end-strands.patch /tmp/bamsurgeon-balance-single-end-strands.patch
RUN cd /opt/reflower/deps/bamsurgeon \
    && patch -p1 < /tmp/bamsurgeon-mark-orphans.patch \
    && patch -p1 < /tmp/bamsurgeon-balance-single-end-strands.patch \
    && rm /tmp/bamsurgeon-mark-orphans.patch \
       /tmp/bamsurgeon-balance-single-end-strands.patch \
    && python -m py_compile \
       bin/addsnv.py \
       bin/addindel.py \
       bin/bamsurgeon/common.py \
       bin/bamsurgeon/aligners.py \
       bin/bamsurgeon/mutation.py \
       bin/bamsurgeon/replace_reads.py
WORKDIR /opt/reflower
ENTRYPOINT ["/opt/reflower/entrypoint.sh"]
CMD ["--help"]
