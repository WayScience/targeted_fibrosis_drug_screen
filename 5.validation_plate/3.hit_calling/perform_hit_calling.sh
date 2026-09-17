#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the hit calling environment
conda activate fibrosis_machine_learning

# convert notebooks to scripts
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# run the notebooks in order
python nbconverted/0.run_mAP_bulk_profiles.py

echo "Hit calling complete."
