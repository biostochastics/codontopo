# codon-topo v0.6.0 reference build environment
#
# Purpose: reproducible container for regenerating every JSON artifact,
# CSV table, figure, and PDF from the codon-topo v0.6.0 release. The
# resulting numerical outputs should match the published values to at
# least three significant figures on any host that runs this image.
#
# Build:
#     docker build -t codontopo:v0.6.0 .
#
# Run (regenerate outputs in a mounted volume):
#     docker run --rm -v "$PWD/output:/repo/output" codontopo:v0.6.0 \
#         codon-topo all --output-dir=/repo/output --seed=135325
#
# Verify (against the published SHA-256 manifest):
#     shasum -a 256 -c output/checksums.sha256
#
# NOTE: PDF rendering (Typst) and figure rendering (R) are separate
# steps. Rebuild them with:
#     typst compile output/manuscript.typ output/manuscript.pdf
#     typst compile output/supplement.typ output/supplement.pdf
#     Rscript src/codon_topo/visualization/R/strengthened_figures.R

FROM python:3.11.14-slim

# System build tools required by scipy/numpy wheels on some platforms,
# plus curl for downloading Typst.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential curl ca-certificates git \
    && rm -rf /var/lib/apt/lists/*

# Typst 0.14.2 for manuscript/supplement PDF rebuilds.
RUN curl -fsSL https://github.com/typst/typst/releases/download/v0.14.2/typst-x86_64-unknown-linux-musl.tar.xz \
        | tar -xJ -C /usr/local/bin --strip-components=1 typst-x86_64-unknown-linux-musl/typst \
    && typst --version

WORKDIR /repo

# Copy the pinned lockfile first for cacheable pip install.
COPY requirements.lock ./
RUN pip install --no-cache-dir -r requirements.lock

# Install codon-topo itself (source copied last so pip layer caches).
COPY pyproject.toml README.md ./
COPY src ./src
COPY tests ./tests
COPY data ./data
COPY scripts ./scripts
RUN pip install --no-cache-dir -e ".[dev]"

# tRNAscan-SE and R are NOT installed in this image — they are optional
# for regenerating tRNA counts and figures respectively. Rerunning
# `codon-topo all` uses the cached tRNAscan output shipped in data/.
# To install tRNAscan-SE 2.0.12 or R 4.5, extend this image.

# Deterministic seed already set as default in cli.
ENV PYTHONHASHSEED=0
CMD ["codon-topo", "--help"]
