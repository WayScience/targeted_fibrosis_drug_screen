#!/usr/bin/env python
# coding: utf-8

# # Generate well-level bulk profiles profiles
# 
# NOTE: We are normalizing the bulk-profile plates to the negative controls as the "standard".

# ## Import libraries

# In[1]:


import pathlib

import pandas as pd

from pycytominer import aggregate, annotate, normalize, feature_select


# ## Set paths and variables

# In[2]:


# Set the plate to process
plate_id = "CARD-CelIns-CX7_260803130001"

# Directory with QC-labeled profiles
qc_labeled_dir = pathlib.Path("./data/qc_labeled_profiles/").resolve(strict=True)

# Directory to save bulk profiles
output_dir = pathlib.Path("./data/bulk_profiles/")
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
# 

# In[3]:


print("Performing preprocessing on", plate_id)

# generating all output file paths
output_annotated_file = str(output_dir / f"{plate_id}_bulk_annotated.parquet")
output_normalized_file = str(output_dir / f"{plate_id}_bulk_normalized.parquet")
output_feature_select_file = str(
    output_dir / f"{plate_id}_bulk_feature_selected.parquet"
)

# loading profiles
profile_df = pd.read_parquet(profile_path)
platemap_df = pd.read_csv(platemap_path)

# Drop all rows in the profiles that failed any Metadata_cqc columns
cqc_columns = [col for col in profile_df.columns if col.startswith("Metadata_cqc")]
if cqc_columns:
    profile_df = profile_df[~profile_df[cqc_columns].any(axis=1)]

# Step 1: Aggregate single-cell data to the well-level using the median
print("Performing aggregation for", plate_id, "...")
aggregated_df = aggregate(
    population_df=profile_df,
    operation="median",
    strata=["Image_Metadata_Plate", "Image_Metadata_Well"],
)

# Step 2: Annotation
print("Performing annotation for", plate_id, "...")
annotate(
    profiles=aggregated_df,
    platemap=platemap_df,
    join_on=["Metadata_well_position", "Image_Metadata_Well"],
    output_type="parquet",
    output_file=output_annotated_file,
)

# Load the annotated parquet file to fix metadata columns names
annotated_df = pd.read_parquet(output_annotated_file)

# Rename columns
annotated_df.rename(columns={"Image_Metadata_Site": "Metadata_Site"}, inplace=True)

# Save annotated profiles back to parquet
annotated_df.to_parquet(output_annotated_file, index=False)

# Step 3: Normalization (mad robustize)
# Normalize using the negative controls as the reference population
print("Performing normalization for", plate_id, "...")
neg_control_query = "Metadata_treatment == 'DMSO' and Metadata_cell_type == 'failing'"
normalize(
    profiles=annotated_df,
    method="mad_robustize",
    samples=neg_control_query,
    output_type="parquet",
    output_file=output_normalized_file,
)

# Step 4: Feature selection
print("Performing feature selection for", plate_id, "...")
feature_select(
    profiles=output_normalized_file,
    operation=feature_select_ops,
    na_cutoff=0,
    blocklist_file="./blocklist_features.txt",
    corr_threshold=0.95,
    freq_cut=0.05,
    output_type="parquet",
    output_file=output_feature_select_file,
)

print(
    f"Aggregation, annotation, normalization, and feature selection complete for {plate_id}"
)


# In[4]:


# Check an example output file
test_df = pd.read_parquet(output_feature_select_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)
