#!/usr/bin/env python
# coding: utf-8

# # Generate well-level bulk profiles profiles
# 
# NOTE: We are normalizing the bulk-profile plates to the negative controls as the "standard".

# ## Import libraries

# In[1]:


import os
import pathlib
import pprint

import pandas as pd

from pycytominer import aggregate, annotate, normalize, feature_select
from pycytominer.cyto_utils import infer_cp_features


# ## Set paths and variables

# In[2]:


# get the batch to process from environment variable
batch_to_process = os.environ.get("BATCH", "batch_1")
if batch_to_process is None:
    raise ValueError(
        "Please set the BATCH environment variable before running this script."
    )

# base directory where batches are located
base_dir = pathlib.Path("./data/").resolve(strict=True)

# Decide what to process
if batch_to_process:
    print(f"Processing {batch_to_process}")
    batch_dirs = [base_dir / batch_to_process]
else:
    print("No specific batch set, processing all available batches")
    batch_dirs = [p for p in base_dir.glob("batch_*") if p.is_dir()]

# path for platemap directory
platemap_dir = pathlib.Path("../metadata/updated_platemaps/")

# Load the barcode_platemap file
barcode_platemap_df = pd.read_csv(
    (platemap_dir / "updated_barcode_platemap.csv").resolve()
)

# operations to perform for feature selection
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "blocklist",
    "drop_na_columns",
]


# ## Set dictionary with plates to process

# In[3]:


plate_info_dictionary = {}

# Loop over batches and layouts
for batch_dir in batch_dirs:
    layouts = [p for p in batch_dir.iterdir() if p.is_dir()]  # all layouts
    for layout_dir in layouts:
        qc_labeled_dir = layout_dir / "qc_labeled_profiles"
        output_dir = layout_dir / "bulk_profiles"
        output_spherize_dir = layout_dir / "spherized_bulk_profiles"

        # Create directories once per layout
        for d in [output_dir, output_spherize_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Extract plate names from parquet files
        parquet_files = list(qc_labeled_dir.glob("*.parquet"))
        plate_names = [
            "_".join(f.stem.split("_")[:2]) if len(f.stem.split("_")) >= 2 else f.stem
            for f in parquet_files
        ]

        for name in plate_names:
            # Find the corresponding parquet file
            matching_files = [f for f in parquet_files if name in f.stem]
            if not matching_files:
                continue
            profile_path = matching_files[0].resolve(strict=True)

            # Find corresponding platemap CSV
            platemap_row = barcode_platemap_df.loc[
                barcode_platemap_df["plate_barcode"] == name
            ]
            if platemap_row.empty:
                raise ValueError(f"No platemap found for plate {name}")
            platemap_path = (
                platemap_dir / f"{platemap_row['platemap_file'].values[0]}.csv"
            ).resolve(strict=True)

            # Add to dictionary
            plate_info_dictionary[name] = {
                "profile_path": profile_path,
                "platemap_path": platemap_path,
                "output_dir": output_dir,
                "spherize_output_dir": output_spherize_dir,
            }

# View dictionary
print("Number of plates to process:", len(plate_info_dictionary))
pprint.pprint(plate_info_dictionary, indent=4)


# ## Process data with pycytominer
# 

# In[4]:


for plate, info in plate_info_dictionary.items():
    output_dir = info["output_dir"]
    print("Performing preprocessing on", plate, output_dir)

    # Use the output_dir from the dictionary for this specific plate
    output_dir = info["output_dir"]
    spherized_output_dir = info["spherize_output_dir"]

    # generating all output file paths using the output_dir from the dictionary
    output_annotated_file = str(output_dir / f"{plate}_bulk_annotated.parquet")
    output_normalized_file = str(output_dir / f"{plate}_bulk_normalized.parquet")
    output_feature_select_file = str(
        output_dir / f"{plate}_bulk_feature_selected.parquet"
    )
    output_spherized_file = str(
        spherized_output_dir / f"{plate}_bulk_spherized.parquet"
    )

    profile_df = pd.read_parquet(info["profile_path"])
    platemap_df = pd.read_csv(info["platemap_path"])

    # Drop all rows in the profiles that failed any Metadata_cqc columns
    cqc_columns = [col for col in profile_df.columns if col.startswith("Metadata_cqc")]
    if cqc_columns:
        profile_df = profile_df[~profile_df[cqc_columns].any(axis=1)]

    # Step 1: Aggregate single-cell data to the well-level using the median
    print("Performing aggregation for", plate, "...")
    aggregated_df = aggregate(
        population_df=profile_df,
        operation="median",
        strata=["Image_Metadata_Plate", "Image_Metadata_Well"],
    )

    print("Performing annotation for", plate, "...")
    # Step 2: Annotation
    annotate(
        profiles=aggregated_df,
        platemap=platemap_df,
        join_on=["Metadata_well_position", "Image_Metadata_Well"],
        output_file=output_annotated_file,
        output_type="parquet",
    )

    # Load the annotated parquet file to fix metadata columns names
    annotated_df = pd.read_parquet(output_annotated_file)

    # Rename columns
    annotated_df.rename(columns={"Image_Metadata_Site": "Metadata_Site"}, inplace=True)

    # Save back
    annotated_df.to_parquet(output_annotated_file, index=False)

    # Step 2: Normalization
    print("Performing normalization for", plate, "...")
    neg_control_query = "Metadata_treatment == 'DMSO' and Metadata_cell_type == 'failing'"
    normalized_df = normalize(
        profiles=output_annotated_file,
        method="standardize",
        output_file=output_normalized_file,
        output_type="parquet",
        samples=neg_control_query,
    )

    # Step 3: Feature selection
    print("Performing feature selection for", plate, "...")
    global_fs_df = feature_select(
        profiles=normalized_df,
        operation=feature_select_ops,
        na_cutoff=0,
        output_type="parquet",
        blocklist_file="./blocklist_features.txt",
        corr_threshold=0.95,
        freq_cut=0.05,
    )

    print(
        f"Aggregation, annotation, normalization, and feature selection complete for {plate}"
    )

    # step 4: spherize profiles
    print("Performing spherization for", plate, "...")
    
    # We perform a second feature selection focused specifically on the negative controls.
    # This is required because spherization (whitening) uses the variation observed in 
    # the control group to define the "baseline" for the whole experiment. 
    # If a feature is "static" in the controls (identical feature value across all wells), 
    # there is no baseline variance to measure against, which causes the spherization 
    # calculation to fail.
    feature_select(
        profiles=global_fs_df,
        operation="variance_threshold",
        freq_cut=0.05,
        unique_cut=0.01,
        samples=neg_control_query,
        output_file=output_feature_select_file,
        output_type="parquet",
    )

    # Spherize using the negative controls as the reference population
    normalize(
        profiles=output_feature_select_file,
        method="spherize",
        output_file=output_spherized_file,
        output_type="parquet",
        samples=neg_control_query,
        spherize_center=True,
        spherize_method="ZCA-cor",
        spherize_epsilon=1e-6,
    )


# In[5]:


# Check output file
test_df = pd.read_parquet(output_feature_select_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)


# In[6]:


# Check output file
test_df = pd.read_parquet(output_annotated_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)

