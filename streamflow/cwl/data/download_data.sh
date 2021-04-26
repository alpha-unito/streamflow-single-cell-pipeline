#!/bin/bash

SCRIPT_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Download genome
curl -fsSL \
    -o "${SCRIPT_DIRECTORY}/refdata-cellranger-GRCh38-3.0.0.tar.gz" \
    https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
tar -xzvf "${SCRIPT_DIRECTORY}/refdata-cellranger-GRCh38-3.0.0.tar.gz" -C ${SCRIPT_DIRECTORY}
rm -f "${SCRIPT_DIRECTORY}/refdata-cellranger-GRCh38-3.0.0.tar.gz"

# Download input
curl -fsSL \
    -o "${SCRIPT_DIRECTORY}/input.zip" \
    https://datacloud.di.unito.it/index.php/s/sXz2brP359sXxFQ/download/input.zip
unzip -d "${SCRIPT_DIRECTORY}" "${SCRIPT_DIRECTORY}/input.zip"
rm -f "${SCRIPT_DIRECTORY}/input.zip"