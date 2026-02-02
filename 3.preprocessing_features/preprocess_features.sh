#!/bin/bash

# --- SET THE BATCH TO PROCESS HERE ---
BATCH="batch_3"
export BATCH
# -----------------------------------

# initialize the correct shell for your machine to allow conda to work
conda init bash
conda activate fibrosis_preprocessing_env

# convert Jupyter notebook(s) to script (optional)
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# --- Harmonize with CytoTable for all layouts per batch ---

# Set base directory for converted data
CONVERTED_DIR="./data/${BATCH}/"

# Flag to determine if we need to run conversion
RUN_CONVERT=false

# Loop over each layout in the batch
for layout_dir in "$CONVERTED_DIR"/*/; do
    # Check if converted_profiles exists
    converted_profiles_dir="${layout_dir}converted_profiles"
    
    # Count parquet files
    parquet_count=0
    if [ -d "$converted_profiles_dir" ]; then
        parquet_count=$(find "$converted_profiles_dir" -maxdepth 1 -name "*.parquet" | wc -l)
    fi

    # If not 4 parquet files, mark to run conversion
    if [ "$parquet_count" -ne 4 ]; then
        RUN_CONVERT=true
        echo "⚠️ Layout $(basename "$layout_dir") is missing parquet files ($parquet_count found)."
    fi
done

# Run conversion if needed
if [ "$RUN_CONVERT" = true ]; then
    echo ">>> Running harmonization for ${BATCH}"
    python nbconverted/0.convert_cytotable.py
else
    echo "✅ All layouts have 4 parquet files. Skipping CytoTable conversion for ${BATCH}."
fi

# --- Single-cell quality control per layout/plate ---
LAYOUT_PLATE_DIR="../2.cellprofiler_processing/cp_output/${BATCH}"

# Make sure the papermill output folder exists
mkdir -p ./papermill_outputs

# Find all layout/plate combos (two levels down)
LAYOUT_PLATE_LIST=$(find "$LAYOUT_PLATE_DIR" -mindepth 2 -maxdepth 2 -type d)

for layout_plate_dir in $LAYOUT_PLATE_LIST; do
    # Extract layout and plate names from the path
    platemap_layout=$(basename $(dirname "$layout_plate_dir"))
    plate_id=$(basename "$layout_plate_dir")
    
    echo ">>> Running single-cell QC for layout ${platemap_layout}, plate ${plate_id}"

    # Run papermill for QC
    papermill 1.sc_quality_control.ipynb \
              ./papermill_outputs/1.sc_quality_control_${platemap_layout}_${plate_id}.ipynb \
              -p platemap_layout "${platemap_layout}" \
              -p plate_id "${plate_id}"
done

# --- Single-cell processing per batch ---
echo ">>> Running single-cell processing for batch ${BATCH}"

for layout_dir in "$CONVERTED_DIR"/*/; do
    platemap_layout=$(basename "$layout_dir")
    echo ">>> Processing single-cell data for layout ${platemap_layout}"
    python nbconverted/2.single_cell_processing.py --batch "${BATCH}" --layout "${platemap_layout}"
done

# --- Aggregate single-cell profiles per batch ---
echo ">>> Aggregating single-cell profiles for batch ${BATCH}"

for layout_dir in "$CONVERTED_DIR"/*/; do
    platemap_layout=$(basename "$layout_dir")
    echo ">>> Aggregating single-cell profiles for layout ${platemap_layout}"
    python nbconverted/3.aggregate_single_cells.py --batch "${BATCH}" --layout "${platemap_layout}"
done

echo ">>> Processing complete for ${BATCH}"
