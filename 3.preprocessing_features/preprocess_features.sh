#!/bin/bash

# -----------------------------
# Initialize environment
# -----------------------------
conda init bash
conda activate fibrosis_preprocessing_env

# convert notebooks to scripts
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# -----------------------------
# Set the batches explicitly
# -----------------------------
BATCH_LIST=("batch_1" "batch_2" "batch_3")

# Loop over all batches sequentially
for BATCH in "${BATCH_LIST[@]}"; do
    export BATCH
    echo "======================================"
    echo "Processing batch: $BATCH"
    echo "======================================"

    # -----------------------------
    # CytoTable conversion (batch level)
    # -----------------------------
    CONVERTED_DIR="./data/${BATCH}"
    RUN_CONVERT=false

    for layout_dir in "$CONVERTED_DIR"/*/; do
        converted_profiles_dir="${layout_dir}/converted_profiles"
        parquet_count=0
        if [ -d "$converted_profiles_dir" ]; then
            parquet_count=$(find "$converted_profiles_dir" -maxdepth 1 -name "*.parquet" | wc -l)
        fi
        if [ "$parquet_count" -ne 4 ]; then
            RUN_CONVERT=true
            echo "⚠️ Layout $(basename "$layout_dir") missing parquet files ($parquet_count found)"
        fi
    done

    if [ "$RUN_CONVERT" = true ]; then
        echo ">>> Running CytoTable conversion for ${BATCH}"
        python nbconverted/0.convert_cytotable.py
    else
        echo "✅ Conversion already complete"
    fi

    # -----------------------------
    # Plate-level QC (parallel per layout)
    # -----------------------------
    LAYOUT_PLATE_DIR="../2.cellprofiler_processing/cp_output/${BATCH}"
    mkdir -p ./papermill_outputs

    for layout_dir in "$LAYOUT_PLATE_DIR"/*/; do
        platemap_layout=$(basename "$layout_dir")
        echo ">>> Running QC for all plates in layout ${platemap_layout}"

        for plate_dir in "$layout_dir"/*/; do
            plate_id=$(basename "$plate_dir")
            echo "   ⏱ Starting QC for plate ${plate_id}"

            papermill 1.sc_quality_control.ipynb \
                ./papermill_outputs/1.sc_quality_control_${platemap_layout}_${plate_id}.ipynb \
                -p platemap_layout "${platemap_layout}" \
                -p plate_id "${plate_id}" &   # run in background
        done

        # Wait for all plates in this layout to finish
        wait
        echo "✅ QC completed for layout ${platemap_layout}"
    done

    # -----------------------------
    # Single-cell processing (batch)
    # -----------------------------
    echo ">>> Running single-cell processing for ${BATCH}"
    python nbconverted/2.single_cell_processing.py

    # -----------------------------
    # Bulk aggregation (batch)
    # -----------------------------
    echo ">>> Aggregating batch ${BATCH}"
    python nbconverted/3.aggregate_single_cells.py

    echo ">>> Completed batch ${BATCH}"
done

echo ">>> All batches processed ✅"
