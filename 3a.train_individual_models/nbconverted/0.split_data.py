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
random_state=0
random.seed(random_state)

# Path to directory with feature selected profiles
path_to_feature_selected_data = pathlib.Path(
    "../3.preprocessing_features/data/single_cell_profiles/"
).resolve(strict=True)

# Find all feature selected parquet files
feature_selected_files = list(path_to_feature_selected_data.glob("*_feature_selected.parquet"))

# Make directory for split data
output_dir = pathlib.Path("./data")
output_dir.mkdir(exist_ok=True)


# ## Load in feature selected data

# In[3]:


# Load in all feature selected files as dataframes
feature_selected_dfs_dict = {
    pathlib.Path(file).stem.split('_')[0]: pd.read_parquet(file) for file in feature_selected_files
}

pprint.pprint(feature_selected_dfs_dict, indent=4)


# ## For each dataframe, take only the DMSO cells and split 70/30 for training and testing

# In[4]:


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

# In[5]:


# Identify common feature columns (excluding Metadata_ columns)
common_features = set.intersection(*[
    set(df.columns[~df.columns.str.startswith("Metadata_")]) 
    for df in feature_selected_dfs_dict.values()
])

# Concat all filtered dataframes while keeping Metadata columns
combined_df = pd.concat(
    [df[list(common_features) + [col for col in df.columns if col.startswith("Metadata_")]] 
     for df in feature_selected_dfs_dict.values()],
    ignore_index=True
)

print(combined_df.shape)
combined_df.head()


# ## Filter combined dataframe for only DMSO treated cells and split into training and testing dataframes

# In[6]:


# Set the ratio of the test data to 30% (training data will be 70%)
test_ratio = 0.30

# Filter only the DMSO treated cells
DMSO_combined_df = combined_df[combined_df.Metadata_treatment == "DMSO"]

# Split the combined data into training and test
training_data, testing_data = train_test_split(
    DMSO_combined_df,
    test_size=test_ratio,
    stratify=DMSO_combined_df[["Metadata_cell_type"]],
    random_state=random_state,
)

# View shapes and example output
print("The testing data contains", testing_data.shape[0], "single-cells.")
print("The training data contains", training_data.shape[0], "single-cells.")

# Save training and test data
training_data.to_parquet(output_dir / "combined_batch1_train.parquet")
testing_data.to_parquet(output_dir / "combined_batch1_test.parquet")

testing_data.head()

