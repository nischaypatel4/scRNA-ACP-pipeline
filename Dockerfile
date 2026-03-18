FROM rocker/r-ver:4.4.0

LABEL maintainer="Nischay Patel <nischaypatel4@gmail.com>"
LABEL description="scRNA-seq ACP pipeline container"

# System dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

RUN R -e "install.packages(c( \
    'Seurat', 'SeuratObject', 'harmony', 'ggplot2', 'patchwork', \
    'dplyr', 'tidyverse', 'Matrix', 'scales', 'pheatmap', \
    'RColorBrewer', 'tibble', 'future', 'future.apply' \
), repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c( \
    'DropletUtils', 'clusterProfiler', 'ReactomePA', \
    'org.Hs.eg.db', 'enrichplot', 'GO.db', \
    'SingleR', 'celldex', 'glmGamPoi' \
), ask=FALSE)"
