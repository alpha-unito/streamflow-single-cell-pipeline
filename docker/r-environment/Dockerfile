FROM r-base:3.6.0
LABEL maintainer="iacopo.colonnelli@unito.it"

RUN apt-get update -y \
    && apt-get install -y -q \
        libcurl4-openssl-dev \
        libgit2-dev \
        libfontconfig1-dev \
        libmariadb-dev \
        libssh2-1-dev \
        libssl-dev \
        libtbb-dev \
        libudunits2-dev \
        libxml2-dev \
        python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install umap-learn \
    && wget -q \
        -O pandoc-2.7.3-1-amd64.deb https://github.com/jgm/pandoc/releases/download/2.7.3/pandoc-2.7.3-1-amd64.deb \
    && dpkg -i pandoc-2.7.3-1-amd64.deb \
    && rm -f pandoc-2.7.3-1-amd64.deb \
    && R -e "install.packages(c(  \
            'clustree',           \
            'data.table',         \
            'devtools',           \
            'doParallel',         \
            'dplyr',              \
            'e1071',              \
            'ggplot2',            \
            'ggrepel',            \
            'gplots',             \
            'grid',               \
            'gridExtra',          \
            'openxlsx',           \
            'plotly',             \
            'Radviz',             \
            'RColorBrewer',       \
            'reshape2',           \
            'reticulate',         \
            'Seurat',             \
            'stringr'             \
        ),                        \
        dependencies=TRUE,        \
        repos='http://cran.rstudio.com/')" \
    && R -e "devtools::install_github('dviraran/SingleR')" \
    && R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" \
    && R -e "BiocManager::install('clusteerProfiler')" \
    && R -e "BiocManager::install('DOSE')" \
    && R -e "BiocManager::install('MAST')" \
    && R -e "BiocManager::install('org.Hs.eg.db')" \
    && R -e "BiocManager::install('org.Mm.eg.db')"

COPY progs/ /r-scripts
WORKDIR /r-scripts

ENTRYPOINT []
CMD bash
