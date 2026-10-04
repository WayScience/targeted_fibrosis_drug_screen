#!/usr/bin/env python
# coding: utf-8

# # Generate well-level bulk profiles
# 
# This is the final step of the module.
# We aggregate the position-corrected single-cell profiles from step 3b to the well level (median per well) and pool the wells of **all plates in all batches**.
# The single-cell profiles are already standardized (z-scored using all cells on each plate) and corrected for plate-position effects, so we do not normalize the bulk profiles again.
# 
# We then perform one feature selection and one sphering on the pooled profiles, with the negative controls (failing-cell DMSO) of all plates as the reference population.
# 
# **Input** (per plate, from step 3b): `<plate>_sc_position_corrected.parquet`
# 
# **Outputs** (in `data/bulk_profiles/`):
# 
# - `bulk_position_corrected_feature_selected.parquet`: pooled well-level profiles of all plates after feature selection
# - `bulk_position_corrected_feature_selected_spherized.parquet`: pooled well-level profiles of all plates after sphering
# 
# Each well keeps its plate, and we add `Metadata_Batch` and `Metadata_Platemap` (from the folder names) so wells can be grouped by batch and platemap downstream.
# Bulk profiles from earlier runs, which were built from uncorrected cells, are not changed.

# ## Import libraries

# In[ ]:


import pathlib
import pprint

import pandas as pd
from pycytominer import aggregate, feature_select, normalize

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
output_spherized_file = output_dir / "bulk_position_corrected_feature_selected_spherized.parquet"

# operations to perform for feature selection
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "blocklist",
    "drop_na_columns",
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

# negative control wells used as the reference population
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


# ## Feature selection and sphering
# 
# We perform a two-layer feature selection procedure on the pooled profiles.
# 
# In the first layer, we apply feature selection to the pooled profiles.
# In the second layer, we remove low variance features in the pooled profiles _for the DMSO control wells only_.
# 
# We then apply a sphering transform using this feature selected data.
# We use the negative-control wells of all plates as the reference population to decorrelate the features to place profiles on a shared control-based covariance scale.
# 
# ### Parameters:
# 
# - `neg_control_query`: pandas query string that selects the control wells used for the second feature selection layer and to fit the sphering. Here, the reference population is failing-cell DMSO wells.

# In[ ]:


# step 1: Apply feature selection on the pooled profiles to get a common set of
# features for sphering.
print("Feature selecting the pooled profiles...")
feature_select_df = feature_select(
    profiles=pooled_df,
    operation=feature_select_ops,
    na_cutoff=0,
    blocklist_file="./blocklist_features.txt",
    corr_threshold=0.95,
    freq_cut=0.05,
    output_file=output_feature_select_file,
    output_type="parquet",
)

# step 2: Remove features with too little variation inside the exact control
# population used to fit spherization.
print("Feature selecting with variance threshold within negative controls only...")
zero_negcon_var_fs_df = feature_select(
    profiles=feature_select_df,
    operation="variance_threshold",
    freq_cut=0.05,
    unique_cut=0.01,
    samples=neg_control_query,
)

# step 3: Spherize/whiten all profiles using the pooled negative controls as the
# reference population.
print("Sphering using the pooled negative controls...")
normalize(
    profiles=zero_negcon_var_fs_df,
    method="spherize",
    samples=neg_control_query,
    spherize_center=True,
    spherize_method="ZCA-cor",
    spherize_epsilon=1e-6,
    output_file=output_spherized_file,
    output_type="parquet",
)

print(f"Saved feature-selected profiles to {output_feature_select_file}")
print(f"Saved spherized profiles to {output_spherized_file}")


# In[ ]:


# Check an example output file
test_df = pd.read_parquet(output_spherized_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)

