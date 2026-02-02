#!/usr/bin/env python
# coding: utf-8

# # Process single cell profiles
# 
# NOTE: We are normalizing the plates for all samples as we only have three wells associated with the healthy controls, which is insufficient for normalization.

# ## Import libraries

# In[1]:


import os
import pathlib
import pprint

import pandas as pd

from pycytominer import annotate, normalize, feature_select


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
        cleaned_dir = layout_dir / "cleaned_profiles"
        output_dir = layout_dir / "single_cell_profiles"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Extract plate names from parquet files
        parquet_files = list(cleaned_dir.glob("*.parquet"))
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
            }

# View dictionary
print("Number of plates to process:", len(plate_info_dictionary))
pprint.pprint(plate_info_dictionary, indent=4)


# ## Process data with pycytominer

# In[4]:


for plate, info in plate_info_dictionary.items():
    output_dir = info["output_dir"]
    print("Performing preprocessing on", plate, output_dir)

    # Use the output_dir from the dictionary for this specific plate
    output_dir = info["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)

    output_annotated_file = str(output_dir / f"{plate}_sc_annotated.parquet")
    output_normalized_file = str(output_dir / f"{plate}_sc_normalized.parquet")
    output_feature_select_file = str(
        output_dir / f"{plate}_sc_feature_selected.parquet"
    )

    profile_df = pd.read_parquet(info["profile_path"])
    platemap_df = pd.read_csv(info["platemap_path"])

    print("Performing annotation for", plate, "...")
    # Step 1: Annotation
    annotate(
        profiles=profile_df,
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
    normalized_df = normalize(
        profiles=output_annotated_file,
        method="standardize",
        output_file=output_normalized_file,
        output_type="parquet",
        samples="Metadata_treatment == 'DMSO' and Metadata_cell_type == 'failing'",
    )

    # Step 3: Feature selection
    print("Performing feature selection for", plate, "...")
    feature_select(
        normalized_df,
        operation=feature_select_ops,
        na_cutoff=0,
        output_file=output_feature_select_file,
        output_type="parquet",
        blocklist_file="./blocklist_features.txt",
    )

    print(f"Annotation, normalization, and feature selection complete for {plate}")


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

