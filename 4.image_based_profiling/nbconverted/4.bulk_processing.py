#!/usr/bin/env python
# coding: utf-8

# # Generate well-level bulk profiles
# 
# This is the final step of the module.
# We aggregate the position-corrected single-cell profiles from step 3b to the well level (median per well) and pool the wells of **all plates in all batches**.
# The single-cell profiles are already standardized (z-scored using all cells on each plate) and corrected for plate-position effects, so we do not normalize the bulk profiles again.
# 
# We then perform one feature selection on the pooled profiles.
# 
# **Input** (per plate, from step 3b): `<plate>_sc_position_corrected.parquet`
# 
# **Output** (in `data/bulk_profiles/`):
# 
# - `bulk_position_corrected_feature_selected.parquet`: pooled well-level profiles of all plates after feature selection
# 
# Each well keeps its plate, and we add `Metadata_Batch` and `Metadata_Platemap` (from the folder names) so wells can be grouped by batch and platemap downstream.
# Bulk profiles from earlier runs, which were built from uncorrected cells, are not changed.

# ## Import libraries

# In[ ]:


import pathlib
import pprint

import pandas as pd
from pycytominer import aggregate, feature_select


# ## Set paths and variables

# In[ ]:


# base directory where batches are located
base_dir = pathlib.Path("./data/").resolve(strict=True)

# pooled outputs of all plates
output_dir = base_dir / "bulk_profiles"
output_dir.mkdir(parents=True, exist_ok=True)
output_feature_select_file = (
    output_dir / "bulk_position_corrected_feature_selected.parquet"
)

# operations to perform for feature selection
# pycytominer applies them in this order, so we first drop features with missing
# values and blocklisted features, and then apply the variance and correlation filters
feature_select_ops = [
    "drop_na_columns",
    "blocklist",
    "variance_threshold",
    "correlation_threshold",
]

# columns that identify a well and are constant within it; used as the aggregation
# strata so that the well-level metadata are kept
well_strata = [
    "Metadata_WellRow",
    "Metadata_WellCol",
    "Metadata_heart_number",
    "Metadata_cell_type",
    "Metadata_heart_failure_type",
    "Metadata_treatment",
    "Metadata_Pathway",
    "Metadata_Plate",
    "Metadata_Well",
]

# negative control wells (failing-cell DMSO), counted to check the pooled table
neg_control_query = "Metadata_treatment == 'DMSO' and Metadata_cell_type == 'failing'"


# ## Set list of plates to process
# 
# We process the position-corrected profiles of every plate in every batch.

# In[ ]:


plate_info_list = [
    {
        "profile_path": profile_path,
        "batch": profile_path.parents[2].name,  # e.g. batch_1
        "platemap": profile_path.parents[1].name,  # e.g. platemap_1
    }
    for profile_path in sorted(
        base_dir.glob(
            "batch_*/platemap_*/single_cell_profiles/*_sc_position_corrected.parquet"
        )
    )
]

# View list
print("Number of plates to process:", len(plate_info_list))
pprint.pprint(plate_info_list[:3], indent=4)


# ## Aggregate to the well level and pool all plates
# 
# For each plate, we compute the median profile of each well.
# QC-failed cells were already removed in step 2 and the profiles are already annotated, so the well-level metadata are kept by using them as the aggregation strata.
# We then concatenate the wells of all plates.

# In[ ]:


well_profiles = []

for plate_info in plate_info_list:
    print("Performing aggregation for", plate_info["profile_path"].name, "...")
    profile_df = pd.read_parquet(plate_info["profile_path"])
    well_df = aggregate(
        population_df=profile_df,
        operation="median",
        strata=well_strata,
    )
    well_df.insert(0, "Metadata_Platemap", plate_info["platemap"])
    well_df.insert(0, "Metadata_Batch", plate_info["batch"])
    well_profiles.append(well_df)

pooled_df = pd.concat(well_profiles, ignore_index=True)
print("Pooled well-level profiles:", pooled_df.shape)
print("Plates:", pooled_df["Metadata_Plate"].nunique())
print("Negative control wells:", len(pooled_df.query(neg_control_query)))


# ## Feature selection
# 
# We perform one feature selection on the pooled profiles of all plates, so that all wells share one set of features.
# 

# In[ ]:


print("Feature selecting the pooled profiles...")
feature_select(
    profiles=pooled_df,
    operation=feature_select_ops,
    # drop every feature with a missing value (the pycytominer default allows 5%)
    na_cutoff=0,
    blocklist_file="./blocklist_features.txt",
    # 0.95 is less strict than the pycytominer default of 0.9, so fewer correlated
    # features are removed and more features remain
    corr_threshold=0.95,
    # 0.05 is the pycytominer default for the most common value of a feature
    freq_cut=0.05,
    output_file=output_feature_select_file,
    output_type="parquet",
)

print(f"Saved feature-selected profiles to {output_feature_select_file}")


# In[ ]:


# Check an example output file
test_df = pd.read_parquet(output_feature_select_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)

