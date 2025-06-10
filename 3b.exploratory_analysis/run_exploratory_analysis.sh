#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the preprocessing environment
conda activate fibrosis_preprocessing_env

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

#############################
# RUN UMAP ANALYSIS
#############################
# run Python scripts for UMAP analysis
python nbconverted/0.UMAP.py

# deactivate conda env and activate R based env
conda deactivate
conda activate r_fibrosis_env

# run R script to generate performance plots
Rscript nbconverted/1.visualize_UMAP.r

# deactivate conda env and activate preprocessing env
conda deactivate
conda activate fibrosis_preprocessing_env

#############################
# RUN COEFFICIENT ANALYSIS
#############################
# run Python scripts for coefficient analysis
python nbconverted/2.extract_coefficients.py

# deactivate conda env and activate R based env
conda deactivate
conda activate r_fibrosis_env

# run R script to generate coefficient plots
Rscript nbconverted/3.visualize_coefficients.r
