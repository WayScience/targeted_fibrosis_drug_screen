import os
import pathlib
import pprint
import sys

import pandas as pd
from pycytominer import feature_select

sys.path.append("../utils")
import position_correction_utils as pcu

# get the batch to process from environment variable
batch_to_process = os.environ.get("BATCH", "batch_1")

# base directory where batches are located
base_dir = pathlib.Path("./data/").resolve(strict=True)

# Decide what to process
if batch_to_process:
    print(f"Processing {batch_to_process}")
    batch_dirs = [base_dir / batch_to_process]
else:
    print("No specific batch set, processing all available batches")
    batch_dirs = [p for p in base_dir.glob("batch_*") if p.is_dir()]

# position correction fit written by 3a.position_correction/fit_position_correction.ipynb
fit_path = pathlib.Path("./3a.position_correction/position_correction_fit.npz").resolve(
    strict=True
)

# operations to perform for feature selection (same as step 2)
# pycytominer applies them in this order, so we first drop features with missing
# values and blocklisted features, and then apply the variance and correlation filters
feature_select_ops = [
    "drop_na_columns",
    "blocklist",
    "variance_threshold",
    "correlation_threshold",
]

fit = pcu.load_fit(fit_path)

print("Platemaps in fit:", sorted(fit["amplitudes"]))
print("Features in fit:", len(fit["features"]))

plate_info_dictionary = {}

for batch_dir in batch_dirs:
    for layout_dir in sorted(p for p in batch_dir.iterdir() if p.is_dir()):
        platemap = int(layout_dir.name.split("_")[1])
        if platemap not in fit["amplitudes"]:
            raise ValueError(
                f"No amplitude for {layout_dir.name} in the position correction fit. "
                "Rerun 3a.position_correction/fit_position_correction.ipynb."
            )

        profile_dir = layout_dir / "single_cell_profiles"
        for normalized_path in sorted(profile_dir.glob("*_sc_normalized.parquet")):
            plate = normalized_path.name.removesuffix("_sc_normalized.parquet")
            plate_info_dictionary[plate] = {
                "normalized_path": normalized_path,
                "platemap": platemap,
                "output_dir": profile_dir,
            }

# View dictionary
print("Number of plates to process:", len(plate_info_dictionary))
pprint.pprint(plate_info_dictionary, indent=4)

for plate, info in plate_info_dictionary.items():
    output_dir = info["output_dir"]
    output_corrected_file = output_dir / f"{plate}_sc_position_corrected.parquet"
    output_feature_select_file = (
        output_dir / f"{plate}_sc_position_corrected_feature_selected.parquet"
    )

    print("Applying position correction to", plate, "(platemap", info["platemap"], ")")
    normalized_df = pd.read_parquet(info["normalized_path"])

    # Step 1: subtract the platemap-scaled tilt from every cell
    corrected_df = pcu.correct_cells(normalized_df, fit, info["platemap"])
    n_corrected = len([f for f in fit["features"] if f in normalized_df.columns])
    print(f"Corrected {n_corrected} features in {len(corrected_df)} cells")
    corrected_df.to_parquet(output_corrected_file, index=False)

    # Step 2: feature selection on the corrected profiles
    print("Performing feature selection for", plate, "...")
    feature_select(
        profiles=corrected_df,
        operation=feature_select_ops,
        na_cutoff=0,
        output_file=str(output_feature_select_file),
        output_type="parquet",
        blocklist_file="./blocklist_features.txt",
    )

    print(f"Position correction and feature selection complete for {plate}")

# Check the last plate: corrected profiles keep the cells and metadata of the normalized profiles
normalized_df = pd.read_parquet(info["normalized_path"])
corrected_df = pd.read_parquet(output_corrected_file)
selected_df = pd.read_parquet(output_feature_select_file)

metadata_cols = [c for c in normalized_df.columns if c.startswith("Metadata_")]
print("Normalized:", normalized_df.shape)
print("Corrected:", corrected_df.shape)
print("Corrected and feature selected:", selected_df.shape)
print(
    "Metadata unchanged:",
    normalized_df[metadata_cols].equals(corrected_df[metadata_cols]),
)

feature_cols = [f for f in fit["features"] if f in normalized_df.columns]
shift = (normalized_df[feature_cols] - corrected_df[feature_cols]).abs()
print(f"Mean absolute shift per feature: {shift.to_numpy().mean():.3f}")
selected_df.head(2)
