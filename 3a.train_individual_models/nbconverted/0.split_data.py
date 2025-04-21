#!/usr/bin/env python
# coding: utf-8

# # Split each plate from the batch data into training, testing, and holdout data

# In[1]:


import pathlib
import random

import pprint

import pandas as pd
from sklearn.model_selection import train_test_split


# ## Set paths and variables

# In[2]:


# Set random state for the whole notebook to ensure reproducibility
random_state = 0
random.seed(random_state)

# Path to directory with feature selected profiles
path_to_feature_selected_data = pathlib.Path(
    "../3.preprocessing_features/data/single_cell_profiles/"
).resolve(strict=True)

# Find all feature selected parquet files
feature_selected_files = list(
    path_to_feature_selected_data.glob("*_feature_selected.parquet")
)

# Make directory for split data
output_dir = pathlib.Path("./data")
output_dir.mkdir(exist_ok=True)


# ## Load in feature selected data

# In[3]:


# Load in all feature selected files as dataframes
feature_selected_dfs_dict = {
    pathlib.Path(file).stem.split("_")[0]: pd.read_parquet(file)
    for file in feature_selected_files
}


# ### Confirm no NaNs in the Pathway column

# In[4]:


for plate, df in feature_selected_dfs_dict.items():
    invalid_rows = df[
        (df["Metadata_treatment"] != "DMSO") & (df["Metadata_Pathway"].isna())
    ]
    if not invalid_rows.empty:
        print(f"NaNs found in Metadata_Pathway for non-DMSO rows in plate: {plate}")
        print(invalid_rows[["Metadata_treatment", "Metadata_Pathway"]])


# ## For each dataframe, take only the DMSO cells and split 70/30 for training and testing

# In[5]:


# Set the ratio of the test data to 30% (training data will be 70%)
test_ratio = 0.30

for plate, df in feature_selected_dfs_dict.items():
    # Filter only the rows with DMSO treatment for model training and testing
    DMSO_df = df[df.Metadata_treatment == "DMSO"]
    print(f"Plate: {plate} contains {DMSO_df.shape[0]} DMSO profiles")

    # Split data into training and test sets
    train_df, test_df = train_test_split(
        DMSO_df,
        test_size=test_ratio,
        stratify=DMSO_df[["Metadata_cell_type"]],
        random_state=random_state,
    )

    # Print the shapes of the training and testing data
    print(f"Training data shape: {train_df.shape}")
    print(f"Testing data shape: {test_df.shape}")

    # Save training and test data
    train_df.to_parquet(output_dir / f"{plate}_train.parquet")
    test_df.to_parquet(output_dir / f"{plate}_test.parquet")


# ## Combine the 4 plates together using the common morphology features

# In[6]:


# Load in the training and testing files from the other plates (do not include combined if it exists or it will be a bug)
train_files = [
    f for f in output_dir.glob("*_train.parquet") if not f.name.startswith("combined")
]
test_files = [
    f for f in output_dir.glob("*_test.parquet") if not f.name.startswith("combined")
]

# Load files
train_dfs = [pd.read_parquet(f) for f in train_files]
test_dfs = [pd.read_parquet(f) for f in test_files]
all_dfs = train_dfs + test_dfs

# Get intersection of feature columns (excluding Metadata_) across the dataframes
common_features = set.intersection(
    *[set(df.columns[~df.columns.str.startswith("Metadata_")]) for df in all_dfs]
)
print(len(common_features), "common features across all dataframes")

# Use metadata columns from first df
metadata_cols = [col for col in all_dfs[0].columns if col.startswith("Metadata_")]
all_cols = metadata_cols + sorted(common_features)

# Reindex with consistent columns
train_dfs = [df.reindex(columns=all_cols) for df in train_dfs]
test_dfs = [df.reindex(columns=all_cols) for df in test_dfs]

# Merge and save
combined_train_df = pd.concat(train_dfs, ignore_index=True)
combined_test_df = pd.concat(test_dfs, ignore_index=True)

combined_train_df.to_parquet(output_dir / "combined_batch1_train.parquet", index=False)
combined_test_df.to_parquet(output_dir / "combined_batch1_test.parquet", index=False)

print("Train shape:", combined_train_df.shape)
print("Test shape:", combined_test_df.shape)

# Print on dataframe to verify
combined_train_df.head()

