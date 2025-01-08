#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the CellProfiler environment
conda activate fibrosis_cp_env

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=scripts/ *.ipynb

# run Python script for CellProfiler segmentation and feature extraction
python scripts/cp_analysis.py
