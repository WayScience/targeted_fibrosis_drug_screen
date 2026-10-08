#!/bin/bash

# initialize and activate env
conda init bash
conda activate fibrosis_cp_env

# convert notebooks to scripts
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# loop through platemap_1 to platemap_4
for i in {1..4}; do
    export PLATEMAP_LAYOUT="platemap_${i}"
    echo ">>> Running analysis for ${PLATEMAP_LAYOUT}"
    python nbconverted/cp_analysis.py
done
