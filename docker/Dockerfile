# RiboPilot pipeline image: R + the RiboPilot package + the full nrDNA tool chain
# on a single PATH. The Nextflow processes call every tool directly (no
# `conda run -n <env>`), so one environment holds everything, unlike MitoPilot's
# multi-env image. Built and pushed with docker/deploy-*.sh; referenced from a
# project's .config via `process.container` (see inst/config.docker).

# 1 bakes the rDNA panel in, 0 builds without it (deploy-*.sh --no-db). Declared
# before the first FROM because only a global ARG can be used in a FROM.
ARG RP_DB=1

FROM condaforge/mambaforge:24.9.2-0 AS base

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        parallel \
        default-jre \
        build-essential \
        procps && \
    rm -rf /var/lib/apt/lists/*

RUN conda config --set channel_priority strict && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# --- Bioinformatics tool chain, pinned; environment.yml carries the same pins.
#     nextflow runs on the host and is not needed inside the image. sqlite
#     provides the sqlite3 CLI the shell/R-script processes call directly.
RUN mamba install -y \
        'python=3.12' \
        'sqlite=3.53.4' \
        'jq=1.8.2' \
        'fastp=0.23.4' \
        'fastqc=0.12.1' \
        'getorganelle=1.7.7.1' \
        'spades=4.1.0' \
        'blast=2.17.0' \
        'bowtie2=2.5.4' \
        'minimap2=2.28' \
        'samtools=1.21' \
        'seqtk=1.5' \
        'seqkit=2.13.0' \
        'pybarrnap=0.5.2' \
        'barrnap=1.10.6' \
        'itsx=1.1.3' \
        'hmmer=3.4' \
        'entrez-direct=26.0' \
        'biopython=1.88' \
        'ncbi-datasets-cli=18.36.0' && \
    mamba clean -a -y

# --- ribovore is not on bioconda (checked: no ribovore or ribotyper package),
#     so install its ribovore-1.0.5 release straight from GitHub instead.
#     Only ribotyper is used here, which needs infernal but not blastn/ribosensor;
#     the full "linux" install pulls those extra tools too but they sit unused.
RUN mkdir -p /opt/ribovore && cd /opt/ribovore && \
    wget -q https://raw.githubusercontent.com/ncbi/ribovore/ribovore-1.0.5/install.sh -O install.sh && \
    echo "a23e58b4a30306dbc4d5e8bba0127c20b016fb1bd0d7921221ffbb25feb97bde  install.sh" > install.sh.sha256 && \
    sha256sum -c install.sh.sha256 && \
    sh install.sh linux && \
    rm install.sh install.sh.sha256
ENV RIBOINSTALLDIR=/opt/ribovore
ENV RIBOSCRIPTSDIR=$RIBOINSTALLDIR/ribovore
ENV RIBOINFERNALDIR=$RIBOINSTALLDIR/infernal/binaries
ENV RIBOEASELDIR=$RIBOINSTALLDIR/infernal/binaries
ENV RIBOSEQUIPDIR=$RIBOINSTALLDIR/sequip
ENV RIBOTIMEDIR=/usr/bin
ENV RIBOBLASTDIR=$RIBOINSTALLDIR/ncbi-blast/bin
ENV RRNASENSORDIR=$RIBOINSTALLDIR/rRNA_sensor
ENV VECPLUSDIR=$RIBOINSTALLDIR/vecscreen_plus_taxonomy
ENV PERL5LIB=$RIBOSCRIPTSDIR:$RIBOSEQUIPDIR:$VECPLUSDIR
ENV PATH=$RIBOSCRIPTSDIR:$PATH

# --- R + the RiboPilot package dependencies (installed from conda so the package
#     itself installs with dependencies = FALSE below). Biostrings is the one
#     Bioconductor dep. Kept in its own layer so tool-chain edits don't rebuild R.
RUN mamba install -y \
        'r-base=4.4.3' \
        'r-remotes=2.5.0' \
        'bioconductor-biostrings=2.74.0' \
        'bioconductor-pwalign=1.2.0' \
        'bioconductor-decipher=3.2.0' \
        'r-cli=3.6.6' \
        'r-dbi=1.3.0' \
        'r-dbplyr=2.6.0' \
        'r-dplyr=1.2.1' \
        'r-ggplot2=4.0.3' \
        'r-glue=1.8.1' \
        'r-htmltools=0.5.9' \
        'r-htmlwidgets=1.6.4' \
        'r-httr2=1.2.3' \
        'r-jsonlite=2.0.0' \
        'r-later=1.4.8' \
        'r-msar=0.6.0' \
        'r-processx=3.9.0' \
        'r-ps=1.9.3' \
        'r-reactable=0.4.5' \
        'r-readr=2.2.0' \
        'r-rsqlite=3.53.3' \
        'r-shiny=1.14.0' \
        'r-shinyjs=2.1.1' \
        'r-shinywidgets=0.9.1' \
        'r-xml2=1.5.0' && \
    mamba clean -a -y

# The panel and its decoys are versioned quarterly while the image versions on
# every package bump, so they sit below the package layer and stay cached.
# ADD extracts a local tarball in place and never stores the archive itself.
#
# Two stages, selected by RP_DB, so an image can be built before the panel
# release exists. Only the selected one is built, so db-0 needs no tarballs in
# the context. db-0 deliberately leaves RIBOPILOT_DB_DIR unset: an image with no
# panel must not advertise one, and discovery then finds nothing and goes remote.
FROM base AS db-0

FROM base AS db-1
ARG RP_DB_VERSION=2026.08.18
ADD docker/ribopilot_rrna-${RP_DB_VERSION}.tar.gz /opt/ribopilot/blastdb/${RP_DB_VERSION}/
ADD docker/decoys-${RP_DB_VERSION}.tar.gz         /opt/ribopilot/blastdb/${RP_DB_VERSION}/
# ADD does not fail on a corrupt-but-gzip-shaped archive; it falls back to a
# plain copy, so an unusable panel would otherwise bake in and only surface in
# somebody's run. The decoys are checked too: without them a search runs
# panel-only and a non-metazoan contaminant comes back unmatched instead of
# identified, which is the exact failure the decoys exist to prevent. Named one
# by one, never a shell loop: an undefined $var in a RUN expands to empty at
# build time, so a loop variable would check nothing.
# No pipe: a piped command returns 141 here even on success.
RUN blastdbcmd -db /opt/ribopilot/blastdb/${RP_DB_VERSION}/ribopilot_rrna -info > /dev/null && \
    test -s /opt/ribopilot/blastdb/${RP_DB_VERSION}/manifest.json && \
    blastdbcmd -db /opt/ribopilot/blastdb/${RP_DB_VERSION}/16S_ribosomal_RNA -info > /dev/null && \
    blastdbcmd -db /opt/ribopilot/blastdb/${RP_DB_VERSION}/SSU_eukaryote_rRNA -info > /dev/null && \
    blastdbcmd -db /opt/ribopilot/blastdb/${RP_DB_VERSION}/LSU_eukaryote_rRNA -info > /dev/null
ENV RIBOPILOT_DB_DIR=/opt/ribopilot/blastdb
# BLAST resolves taxonomy4blast.sqlite3 through $BLASTDB or the WORKING
# DIRECTORY, never through the -db path. Without this, -taxids fails from every
# cwd except the panel directory -- and per spec section 4 it does not fail
# loudly, it returns UNFILTERED rows. Nextflow runs blastn from a task work
# directory, so the panel dir is exactly where it never is.
ENV BLASTDB=/opt/ribopilot/blastdb/${RP_DB_VERSION}
# Verified from / and under the env the image actually ships with, because the
# defect this catches is "all three taxonomy files present, BLASTDB unset" --
# which every check above passes. A check that cannot fail the way production
# fails is not a check. txid33208 is the panel's own root, so it cannot be
# legitimately empty. Not folded into the RUN above: that one is extracted and
# replayed against a synthetic panel by the test suite, which has no taxonomy
# files. Does not start with `RUN blastdbcmd` for the same reason.
RUN cd / && \
    blastdbcmd -db /opt/ribopilot/blastdb/${RP_DB_VERSION}/ribopilot_rrna \
      -taxids 33208 -outfmt "%a" > /tmp/taxids.check && \
    test -s /tmp/taxids.check && \
    rm -f /tmp/taxids.check

FROM db-${RP_DB} AS final

# --- Install RiboPilot last so R-code edits only invalidate this cheap layer.
#     Dependencies are already provided by conda above, so skip re-resolving them.
#     GIT_VERSION_TAG moves here, immediately above the layer it labels: every
#     layer below an ARG that changes on each publish rebuilds on each publish,
#     and the database layers above must not pay that cost.
ARG GIT_VERSION_TAG=unspecified
COPY docker/RiboPilot_*.tar.gz /pkg.tar.gz
RUN R -e 'remotes::install_local("/pkg.tar.gz", dependencies = FALSE, upgrade = "never")' && \
    R -e 'library(RiboPilot)' && \
    rm /pkg.tar.gz

# With RP_DB=1 (the default) the rDNA panel and decoy databases are baked in
# above at RIBOPILOT_DB_DIR, so a default install works with zero configuration;
# RIBOPILOT_DB_DIR can still be overridden to point at a user-built directory or
# a different version. With RP_DB=0 nothing is baked in and every search is
# remote unless the user configures or mounts a panel.
