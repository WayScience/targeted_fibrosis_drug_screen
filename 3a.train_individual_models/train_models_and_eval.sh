#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the machine learning environment
conda activate fibrosis_machine_learning

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# run Python scripts to split data into training and test, train models, and extract performance metrics
python nbconverted/0.split_data.py
python nbconverted/1.train_models.py
python nbconverted/2.extract_model_performance.py
python nbconverted/4.apply_models_to_compounds.py

# deactivate conda env and activate R based env
conda deactivate 
conda activate r_fibrosis_env

# run R script to generate performance plots
Rscript nbconverted/3.vis_model_performance.r
Rscript nbconverted/5.vis_probabilities.r
