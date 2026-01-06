#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the preprocessing environment
conda activate fibrosis_preprocessing_env

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# --- Harmonize CellProfiler outputs to parquet files ---
for i in {1..4}; do
    export PLATEMAP_LAYOUT="platemap_${i}"
    CONVERTED_DIR="./data/${PLATEMAP_LAYOUT}/converted_profiles"

    # Check if folder exists AND has any files/subfolders
    if [ -d "$CONVERTED_DIR" ] && [ "$(find "$CONVERTED_DIR" -mindepth 1 -print -quit)" ]; then
        echo "⚠️ Skipping ${PLATEMAP_LAYOUT}: ${CONVERTED_DIR} already exists and contains files."
        echo "   CytoTable will not be rerun for this platemap."
        continue  # Skip to the next platemap
    fi

    echo ">>> Running analysis for ${PLATEMAP_LAYOUT}"
    python nbconverted/0.convert_cytotable.py
done

echo ">>> Harmonization of CellProfiler outputs complete for platemap_1 to platemap_4."

# --- Single-cell quality control ---
# loop through platemap_1 to platemap_4
for i in {1..4}; do
    export PLATEMAP_LAYOUT="platemap_${i}"
    echo ">>> Running single-cell quality control for ${PLATEMAP_LAYOUT}"
    python nbconverted/1.sc_quality_control.py
done
echo ">>> Single-cell quality control complete for platemap_1 to platemap_4."

# --- Single-cell processing ---
# loop through platemap_1 to platemap_4
for i in {1..4}; do
    export PLATEMAP_LAYOUT="platemap_${i}"
    echo ">>> Running single-cell processing for ${PLATEMAP_LAYOUT}"
    python nbconverted/2.single_cell_processing.py
done
echo ">>> Single-cell processing complete for platemap_1 to platemap_4."

# --- Aggregate single-cell profiles ---
# loop through platemap_1 to platemap_4
for i in {1..4}; do
    export PLATEMAP_LAYOUT="platemap_${i}"
    echo ">>> Aggregating single-cell profiles for ${PLATEMAP_LAYOUT}"
    python nbconverted/3.aggregate_single_cells.py
done
echo ">>> Aggregation of single-cell profiles complete for platemap_1 to platemap_4."
