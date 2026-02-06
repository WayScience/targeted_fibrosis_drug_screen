#!/bin/bash
shopt -s globstar # enable ** globbing

# initialize shell for conda
conda init bash
conda activate fibrosis_cp_env

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# define the batches you want to process
valid_batches=("batch_1" "batch_2" "batch_3")

# base path where batches live
base_path="/media/18tbdrive/CFReT_screening_data/compound_screen"

# output path for corrected images
output_path="./Corrected_Images"

# check if updated_pipelines already exists and has 11 pipeline files
if [ -d "./pipeline/updated_pipelines" ]; then
    pipeline_count=$(find "./pipeline/updated_pipelines" -type f -name "*.cppipe" | wc -l)
    if [ "$pipeline_count" -eq 11 ]; then
        echo "✅  Found 11 pipelines in ./pipeline/updated_pipelines — skipping 0.generate_pipelines.py"
    else
        echo "⚙️  Running 0.generate_pipelines.py (found $pipeline_count pipelines)"
        python nbconverted/0.generate_pipelines.py
    fi
else
    echo "⚙️  Running 0.generate_pipelines.py (no updated_pipelines directory found)"
    python nbconverted/0.generate_pipelines.py
fi

# loop over the valid batches
for batch_name in "${valid_batches[@]}"; do
    batch_dir="$base_path/$batch_name"
    batch_output_dir="$output_path/$batch_name"

    # skip if the batch folder doesn't exist
    if [ ! -d "$batch_dir" ]; then
        echo "⚠️  Batch folder $batch_name does not exist, skipping."
        continue
    fi

    # skip if the output folder already exists and has any TIFF/TIF files
    if [ -d "$batch_output_dir" ] && find "$batch_output_dir" -type f -iregex '.*\.tif(f)?$' -print -quit | grep -q .; then
        echo "🟡  Output folder for $batch_name already contains TIFF files, skipping."
        continue
    fi

    echo "Processing batch: $batch_name"

    # run the notebook/script for this batch
    python nbconverted/1.illum_correct.py "$batch_name"

    echo "Finished batch: $batch_name"
done

# deactivate conda environment
conda deactivate
echo "All done with illumination correction!"
