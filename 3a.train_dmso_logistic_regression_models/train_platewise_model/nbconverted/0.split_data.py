#!/usr/bin/env python
# coding: utf-8

# # Split the DMSO single cells per each plate into training and testing for plate specific modeling of type label

# In[1]:


import pathlib
import random
import json

import numpy as np
import polars as pl

from cfret_ml.data_split_utils import (
    string_to_int_seed,
    salt_seed,
    stratified_fold_split,
    split_summary
)


# ## Set paths, variables and Helpers

# In[2]:


# Set random state for the whole notebook to ensure reproducibility
random_state = 0
random.seed(random_state)
np.random.seed(random_state)
k_fold_k = 5

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


# ## Iterate over plate feature selected data, write fold datasplit indices along with DMSO subset

# In[3]:


cell_type_classes = []
for batch_folder in batch_folders:

    print(f"Processing batch: {batch_folder.name}")

    feature_selected_files = list(
        batch_folder.glob("**/single_cell_profiles/*_feature_selected.parquet")
    )

    if not feature_selected_files:
        print(f"\tNo profiles found in {batch_folder}. Skipping.")
        continue

    for file in feature_selected_files:

        plate_repr = "_".join(pathlib.Path(file).stem.split("_")[:2])            
        plate_output_dir = datasplit_dir / plate_repr
        plate_output_dir.mkdir(exist_ok=True, parents=True)
        plate_name_salt: int = string_to_int_seed(plate_repr)
        salted_seed: int = salt_seed(random_state, plate_name_salt)

        try:
            DMSO_df = (
                pl.scan_parquet(file)
                .filter(pl.col("Metadata_treatment") == "DMSO")
                .collect(engine="cpu")
                .to_pandas()
            )
        except Exception as e:
            print(f"\tError reading/filtering {file.stem}: {e}. Skipping this file.")
            continue

        print(f"\tProcessing {file.stem} for plate {plate_repr}. ")

        if DMSO_df['Metadata_cell_type'].nunique() < 2:
            print(f"\t>Only one cell type present in {file.stem}. Skipping stratified split.")
            continue

        cell_type_classes.extend(DMSO_df["Metadata_cell_type"].unique())

        split = stratified_fold_split(
            DMSO_df, 
            group_col="Metadata_Well", 
            class_col="Metadata_cell_type", 
            random_state=salted_seed
        )

        for i, (train_index, test_index) in enumerate(split):
            fold_record = {
                "fold": i,
                "train_index": train_index.tolist(),
                "test_index": test_index.tolist(),
                "seed": salted_seed
            }
            fold_record.update(split_summary(
                    DMSO_df, 
                    train_index, 
                    test_index,
                    group_col="Metadata_Well",
                    label_col="Metadata_cell_type"
                )
            )
            fold_documentation_path = plate_output_dir / f"fold_{i}_split.json"
            with open(fold_documentation_path, "w") as f:
                json.dump(fold_record, f, indent=4)

        print("\t>Datasplits written.")

        DMSO_df.to_parquet(plate_output_dir / "DMSO.parquet", index=False)
        print(f"\t>DMSO profiles written for {plate_repr}.")


# ## Save global encoding scheme from collected labels to ensure consistent model fitting and interpretation

# In[4]:


unique_cell_types = sorted(set(cell_type_classes))

# define consistent encodings for downstream modeling
cell_type_encoding = {ct: i for i, ct in enumerate(unique_cell_types)}
print("Cell type encoding:")
for ct, enc in cell_type_encoding.items():
    print(f"\t{ct}: {enc}")

encoding_path = datasplit_dir / "cell_type_encoding.json"
with open(encoding_path, "w") as f:
    json.dump(cell_type_encoding, f)

