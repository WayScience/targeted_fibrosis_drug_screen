#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the ML environment
conda activate fibrosis_machine_learning

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=scripts/ *.ipynb

# run Python script for applying models to data
python scripts/0.apply_models.py

# Deactivate environment and activate R environment
conda deactivate
conda activate r_fibrosis_env

# run R script to visualize the model results
Rscript scripts/1.vis_probabilities.r
