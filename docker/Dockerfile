# RiboPilot pipeline image: R + the RiboPilot package + the full nrDNA tool chain
# on a single PATH. The Nextflow processes call every tool directly (no
# `conda run -n <env>`), so one environment holds everything, unlike MitoPilot's
# multi-env image. Built and pushed with docker/deploy-*.sh; referenced from a
# project's .config via `process.container` (see inst/config.docker).

FROM condaforge/mambaforge:24.9.2-0

ARG GIT_VERSION_TAG=unspecified
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        jq \
        parallel \
        default-jre \
        build-essential \
        procps && \
    rm -rf /var/lib/apt/lists/*

RUN conda config --set channel_priority strict && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# --- Bioinformatics tool chain (mirrors environment.yml; nextflow itself runs on
#     the host and is not needed inside the image). sqlite provides the sqlite3
#     CLI the shell/R-script processes call directly. blast>=2.12 has TLS for the
#     HTTPS-only NCBI remote BLAST service.
RUN mamba install -y \
        'sqlite' \
        'fastp=0.23.4' \
        'fastqc' \
        'getorganelle=1.7.7.1' \
        'spades=4.1.0' \
        'blast>=2.12' \
        'bowtie2=2.5.4' \
        'minimap2=2.28' \
        'samtools=1.21' \
        'seqtk' \
        'seqkit' \
        'pybarrnap' \
        'barrnap' \
        'itsx' \
        'hmmer' \
        'entrez-direct' \
        'biopython' \
        'taxonkit' \
        'ncbi-datasets-cli' && \
    mamba clean -a -y

# --- R + the RiboPilot package dependencies (installed from conda so the package
#     itself installs with dependencies = FALSE below). Biostrings is the one
#     Bioconductor dep. Kept in its own layer so tool-chain edits don't rebuild R.
RUN mamba install -y \
        'r-base=4.4.*' \
        'r-remotes' \
        'bioconductor-biostrings' \
        'bioconductor-pwalign' \
        'bioconductor-decipher' \
        'r-cli' \
        'r-config' \
        'r-crayon' \
        'r-dbi' \
        'r-dbplyr' \
        'r-dplyr' \
        'r-gggenes' \
        'r-ggplot2' \
        'r-glue' \
        'r-golem' \
        'r-jsonlite' \
        'r-later' \
        'r-msar' \
        'r-processx' \
        'r-purrr' \
        'r-reactable' \
        'r-readr' \
        'r-rsqlite' \
        'r-shiny' \
        'r-shinyjs' \
        'r-shinywidgets' \
        'r-stringr' && \
    mamba clean -a -y

# --- Install RiboPilot last so R-code edits only invalidate this cheap layer.
#     Dependencies are already provided by conda above, so skip re-resolving them.
COPY docker/RiboPilot_*.tar.gz /pkg.tar.gz
RUN R -e 'remotes::install_local("/pkg.tar.gz", dependencies = FALSE, upgrade = "never")' && \
    R -e 'library(RiboPilot)' && \
    rm /pkg.tar.gz

# Reference databases are resolved at run time (shipped in the installed package,
# a user-built directory, or GitHub URLs), so none are baked into the image.
