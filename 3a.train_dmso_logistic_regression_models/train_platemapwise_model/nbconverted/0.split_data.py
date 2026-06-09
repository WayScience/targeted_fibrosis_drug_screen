#!/usr/bin/env python
# coding: utf-8

# # Split the DMSO single cells per each platemap into training and testing splits
# Holding out whole plates across folds

# In[1]:


import pathlib
import random
import json

import pandas as pd
import numpy as np
import polars as pl

from cfret_ml.data_split_utils import (
    string_to_int_seed,
    salt_seed,
    split_summary
)


# ## Set paths and data split config variables

# In[2]:


# Set random state for the whole notebook to ensure reproducibility
random_state = 0
random.seed(random_state)
np.random.seed(random_state)

# Path to directory with feature selected profiles
path_to_feature_selected_data = pathlib.Path().home() / "mnt" / "bandicoot" /\
    "CFReT_screening_data" / "screen_profiles"
if not path_to_feature_selected_data.exists() and\
    not path_to_feature_selected_data.is_dir():
    raise FileNotFoundError(
        f"Directory {path_to_feature_selected_data} does not exist or is not a directory."
    )

# Find all batch folders
batch_folders = list(path_to_feature_selected_data.glob("batch*"))
batch_folders = [folder for folder in batch_folders if folder.is_dir()]
if not batch_folders:
    raise FileNotFoundError(
        f"No batch folders found in {path_to_feature_selected_data}."
    )

# Make directory for split data
datasplit_dir = pathlib.Path(".") / "datasplits"
datasplit_dir.mkdir(exist_ok=True, parents=True)


# ## Data splitting holding out one plate at a time

# ### Collect DMSO profiles for each platemap

# In[3]:


profiles_by_platemap = {}

for batch_folder in batch_folders:

    print(f"Processing batch: {batch_folder.name}")

    platemap_folders = list(
        f 
        for f in batch_folder.glob("*")
        if f.is_dir() 
    )

    for platemap_folder in platemap_folders:

        platemap_repr = platemap_folder.name
        if not platemap_repr.startswith("platemap_"):
            continue
        sc_profile_folder = (platemap_folder / "single_cell_profiles").resolve()
        if not sc_profile_folder.exists() or not sc_profile_folder.is_dir():
            print(f"\tNo single_cell_profiles folder for {platemap_repr}. Skipping.")
            continue
        if not platemap_repr in profiles_by_platemap:
            profiles_by_platemap[platemap_repr] = []

        plate_files = list(
            f
            for f in sc_profile_folder.glob("*_sc_feature_selected.parquet")
            if f.is_file()
        )

        for plate_file in plate_files:

            DMSO_df = (
                pl.scan_parquet(plate_file)
                .filter(pl.col("Metadata_treatment") == "DMSO")
                .collect(engine="cpu")
                .to_pandas()
            )
            if DMSO_df.empty:
                print(f"\t\tNo DMSO rows for {plate_file.stem}. Skipping.")
                continue

            profiles_by_platemap[platemap_repr].append(DMSO_df)


for platemap, profiles in profiles_by_platemap.items():
    print(f"Collected {len(profiles)} profiles for {platemap}")       



# In[11]:


control_cell_counts = []
for platemap, profiles in profiles_by_platemap.items():
    for profile in profiles:
        control_cell_counts.append(
            profile[
                profile["Metadata_cell_type"] == "failing"
            ].groupby(
                ['Metadata_treatment', 'Metadata_Well', 'Metadata_Plate']
            ).size().reset_index(name='row_count')
        )

control_cell_counts = pd.concat(control_cell_counts, ignore_index=True)
control_cell_counts.to_csv(datasplit_dir / 'control_cell_counts.csv', index=False)


# ### Iterate over platemap feature selected data and produce fold splits

# In[4]:


cell_type_classes = []
for platemap, profiles in profiles_by_platemap.items():

    platemap_output_dir = datasplit_dir / platemap
    platemap_output_dir.mkdir(exist_ok=True, parents=True)
    platemap_salt = string_to_int_seed(platemap)
    salted_seed = salt_seed(random_state, platemap_salt)

    # drop all columns that have missing values in any of the profile columns
    # this removes the non-overlapping feature sets across plates of the
    # same platemap
    platemap_DMSO_df = pd.concat(profiles, ignore_index=True)
    cols_to_drop = [
        c for c in platemap_DMSO_df.columns
        if not c.startswith("Metadata_") and platemap_DMSO_df[c].isna().any()
    ]
    platemap_DMSO_df.drop(columns=cols_to_drop, inplace=True)
    print(f"Dropped {len(cols_to_drop)} columns:", cols_to_drop)

    platemap_DMSO_df.to_parquet(platemap_output_dir / f"DMSO.parquet", index=False)
    cell_type_classes.extend(platemap_DMSO_df["Metadata_cell_type"].unique())

    # Holding out one plate at a time for splitting.
    # This is signifcantly simpler compared to the stratified strategy
    # used to create datasplits for single plate models as we only
    # only need to holdout one plate at a time instead of one well per class.
    unique_plates = platemap_DMSO_df["Metadata_Plate"].unique()
    unique_plates = sorted(unique_plates)
    unique_plates_shuffled = unique_plates.copy()
    random.Random(salted_seed).shuffle(unique_plates_shuffled)

    for i, plate in enumerate(unique_plates_shuffled):

        train_index = platemap_DMSO_df[platemap_DMSO_df["Metadata_Plate"] != plate].index.to_numpy()
        test_index = platemap_DMSO_df[platemap_DMSO_df["Metadata_Plate"] == plate].index.to_numpy()

        # human-readable split summary for debugging and documentation
        # processes, counts number of class and group assignments
        fold_record = {
            "fold": i,
            "train_index": train_index.tolist(),
            "test_index": test_index.tolist(),
            "seed": salted_seed,
        }
        fold_record.update(
            split_summary(
                platemap_DMSO_df, 
                train_index, 
                test_index,
                group_col="Metadata_Plate",
                label_col="Metadata_cell_type"
            )
        )

        with open(platemap_output_dir / f"fold_{i}_split.json", "w") as f:
            json.dump(fold_record, f, indent=4)


# ## Save global encoding scheme from collected labels to ensure consistent model fitting and interpretation

# In[5]:


unique_cell_types = sorted(set(cell_type_classes))

# define consistent encodings for downstream modeling
cell_type_encoding = {ct: i for i, ct in enumerate(unique_cell_types)}
print("Cell type encoding:")
for ct, enc in cell_type_encoding.items():
    print(f"\t{ct}: {enc}")

encoding_path = datasplit_dir / "cell_type_encoding.json"
with open(encoding_path, "w") as f:
    json.dump(cell_type_encoding, f)

