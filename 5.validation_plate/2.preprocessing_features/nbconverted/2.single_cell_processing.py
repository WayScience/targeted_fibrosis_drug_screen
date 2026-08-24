#!/usr/bin/env python
# coding: utf-8

# # Process single cell profiles
# 
# NOTE: We are normalizing the plates for all samples as we only have three wells associated with the healthy controls, which is insufficient for normalization.

# ## Import libraries

# In[1]:


import pathlib

import pandas as pd

from pycytominer import annotate, normalize, feature_select


# ## Set paths and variables

# In[2]:


# Set the plate to process
plate_id = "CARD-CelIns-CX7_260803130001"

# Directory with QC-labeled profiles
qc_labeled_dir = pathlib.Path("./data/qc_labeled_profiles/").resolve(strict=True)

# Directory to save single-cell profiles
output_dir = pathlib.Path("./data/single_cell_profiles/")
output_dir.mkdir(parents=True, exist_ok=True)

# Path to the platemap for the validation plate
platemap_path = pathlib.Path(
    "../metadata/platemaps/platemap_validation.csv"
).resolve(strict=True)

# Path to the QC-labeled profile for the validation plate
profile_path = (qc_labeled_dir / f"{plate_id}_qc_labeled.parquet").resolve(
    strict=True
)

# operations to perform for feature selection
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "blocklist",
    "drop_na_columns",
]


# ## Process data with pycytominer

# In[3]:


print("Performing preprocessing on", plate_id)

output_annotated_file = str(output_dir / f"{plate_id}_sc_annotated.parquet")
output_normalized_file = str(output_dir / f"{plate_id}_sc_normalized.parquet")
output_feature_select_file = str(
    output_dir / f"{plate_id}_sc_feature_selected.parquet"
)

profile_df = pd.read_parquet(profile_path)
platemap_df = pd.read_csv(platemap_path)

# Drop all rows in the profiles that failed any Metadata_cqc columns
cqc_columns = [col for col in profile_df.columns if col.startswith("Metadata_cqc")]
if cqc_columns:
    profile_df = profile_df[~profile_df[cqc_columns].any(axis=1)]

print("Performing annotation for", plate_id, "...")
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
    samples="all",
)

# Step 3: Feature selection
print("Performing feature selection for", plate_id, "...")
feature_select(
    profiles=normalized_df,
    operation=feature_select_ops,
    na_cutoff=0,
    output_file=output_feature_select_file,
    output_type="parquet",
    blocklist_file="./blocklist_features.txt",
)

print(f"Annotation, normalization, and feature selection complete for {plate_id}")


# In[4]:


# Check output file
test_df = pd.read_parquet(output_feature_select_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)


# In[5]:


# Check output file
test_df = pd.read_parquet(output_annotated_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)

