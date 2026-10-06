#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the CellProfiler environment
conda activate fibrosis_cp_env

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# Check if qc_results already exists and contains 44 Image.csv files
if [ -d "./qc_results" ]; then
    image_csv_count=$(find "./qc_results" -type f -name "Image.csv" | wc -l)
    if [ "$image_csv_count" -eq 44 ]; then
        echo "✅  Found 44 Image.csv files in qc_results — skipping 0.run_image_qc.py"
    else
        echo "⚙️  Running 0.run_image_qc.py (found $image_csv_count Image.csv files)"
        python nbconverted/0.run_image_qc.py
    fi
else
    echo "⚙️  Running 0.run_image_qc.py (no qc_results directory found)"
    python nbconverted/0.run_image_qc.py
fi

# run Python script for IC processing
python nbconverted/1.find_qc_thresholds.py
