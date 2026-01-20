#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the CellProfiler environment
conda activate fibrosis_cp_env

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# # ----- OPTIONAL: rerun QC metrics calculation notebook -----
# python nbconverted/0.rerun_qc_metrics.py

# run Python script for IC processing
python nbconverted/1.run_image_qc.py

# # deactivate the conda environment
# conda deactivate
# # activate preprocessing environment
# conda activate fibrosis_preprocessing_env

# # run Python script for generating summary plots
# python nbconverted/2.qc_report.py

# # deactivate the conda environment
# conda deactivate
# # activate the R environment
# conda activate fibrosis_r_env

# # run R script for generating platemaps with QC report
# Rscript nbconverted/3.qc_platemaps.r
