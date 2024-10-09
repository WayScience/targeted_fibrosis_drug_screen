#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the CellProfiler environment
conda activate fibrosis_cp_env

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=scripts/ *.ipynb

# run Python script for IC processing
python scripts/0.illum_correct.py

# Deactivate environment and activate R environment
conda deactivate
conda activate r_fibrosis_env

# run R script to generate QC report
Rscript scripts/1.qc_report.r
