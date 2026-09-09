#!/usr/bin/env python
# coding: utf-8

# # mAP Analysis bulk profiles
# 
# > Note: This notebook and code are updated from Erik Serrano.
# 
# This notebook computes mean Average Precision (mAP) scores for every compound treatment in the validation plate using bulk morphological profiles.
# 
# What is mAP? mAP measures how reproducibly a compound produces a distinct morphological phenotype. For each treatment, replicate wells are ranked against wells from all other treatments. A score near 1 means the treatment consistently looks different from everything else; a score near 0 means it is indistinguishable from noise.

# In[1]:


import pathlib

import numpy as np
import pandas as pd
from copairs import map
from copairs.matching import assign_reference_index
from pycytominer.cyto_utils import infer_cp_features


# In[ ]:


# Set output directory for mAP scores
output_dir = pathlib.Path("./mAP_scores")
output_dir.mkdir(parents=True, exist_ok=True)

# Set the directory containing the bulk profiles from validation plate
screen_profiles_dir = pathlib.Path("../2.preprocessing_features/data/bulk_profiles").resolve(strict=True)

# Find the one bulk spherized profile file in the directory
profile_files = list(screen_profiles_dir.glob("*bulk_spherized.parquet"))
if len(profile_files) != 1:
    raise ValueError("Expected exactly one profile file in the directory")
profile_file = profile_files[0]


# ## Compute mAP for the Validation Plate
# 
# 1. Split DMSO controls by cell type → `DMSO_failing` and `DMSO_nonfailing`
# 2. Mark DMSO_failing wells as the morphological reference (all other wells get index -1)
# 3. Compute per-well Average Precision (AP) — how well each well's replicates rank above all other treatments
# 4. Aggregate AP into mAP per treatment; correct p-values for multiple testing (threshold = 0.05)
# 
# This returns an mAP score per treatment.

# In[3]:


# Pairing configuration for AP and mAP
reference_col = "Metadata_reference_index"
pos_sameby = ["Metadata_treatment", reference_col]
pos_diffby = []
neg_sameby = []
neg_diffby = ["Metadata_treatment", reference_col]

# load the bulk feature-selected profiles for the validation plate
df = pd.read_parquet(profile_file)

# Split DMSO controls by cell type so reference assignment is biologically consistent.
df.loc[
    (df["Metadata_cell_type"] == "failing")
    & (df["Metadata_treatment"] == "DMSO"),
    "Metadata_treatment",
] = "DMSO_failing"
df.loc[
    (df["Metadata_cell_type"] == "nonfailing")
    & (df["Metadata_treatment"] == "DMSO"),
    "Metadata_treatment",
] = "DMSO_nonfailing"

# assign reference index for DMSO_failing controls
df = assign_reference_index(
    df=df,
    condition="Metadata_treatment == 'DMSO_failing'",
    reference_col=reference_col,
    default_value=-1,
)

# infer cell profiles metadata and morphology features
meta_cols = infer_cp_features(population_df=df, metadata=True)
feature_cols = infer_cp_features(population_df=df, metadata=False)

# calculate average precision scores per replicate and control pairing configuration
activity_ap = map.average_precision(
    meta=df[meta_cols],
    feats=df[feature_cols].to_numpy(),
    pos_sameby=pos_sameby,
    pos_diffby=pos_diffby,
    neg_sameby=neg_sameby,
    neg_diffby=neg_diffby,
).copy()

# Aggregate AP into mAP with copairs built-in function.
activity_map = map.mean_average_precision(
    activity_ap,
    pos_sameby,
    null_size=1000,
    threshold=0.05,
    seed=0,
)

# store all compound mAP scores
activity_map.to_parquet(output_dir / "map_scores.parquet", index=False)

# sort by mAP and reset index
activity_map = activity_map.sort_values(
    "mean_average_precision", ascending=False
).reset_index(drop=True)

# display
print("mAP scores shape: ", activity_map.shape)
activity_map.head(14)

