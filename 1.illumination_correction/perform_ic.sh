#!/bin/bash

# initialize shell for conda
conda init bash
conda activate fibrosis_cp_env

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# define the batches you want to process
valid_batches=("batch_1" "batch_2" "batch_3")

# base path where batches live
base_path="/media/18tbdrive/CFReT_screening_data/compound_screen"

# loop over the valid batches
for batch_name in "${valid_batches[@]}"; do
    batch_dir="$base_path/$batch_name"

    # skip if the batch folder doesn't exist
    if [ ! -d "$batch_dir" ]; then
        echo "⚠️  Batch folder $batch_name does not exist, skipping."
        continue
    fi

    echo "Processing batch: $batch_name"

    # run the notebook/script for this batch
    python nbconverted/illum_correct.py "$batch_name"

    echo "Finished batch: $batch_name"
done

# deactivate conda environment
conda deactivate
