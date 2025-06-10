#!/usr/bin/env python

# # Extract UMAP embeddings for each plate of CellProfiler features
#
# NOTE: We are using the feature selected profiles per plate.

# In[1]:


import glob
import pathlib

import numpy as np
import pandas as pd
import umap
from pycytominer.cyto_utils import infer_cp_features

# ## Generate Embeddings

# ### Set constant for whole plates

# In[2]:


# Set constants
umap_random_seed = 0
umap_n_components = 2

# Set embeddings directory
output_dir = pathlib.Path("results")
output_dir.mkdir(parents=True, exist_ok=True)


# ### Set paths to all plates

# In[3]:


# Set input path with single-cell profiles
data_dir = pathlib.Path("../3.preprocessing_features/data/single_cell_profiles")

# Select only the feature selected files
file_suffix = "*sc_feature_selected.parquet"

# Obtain file paths for all feature selected plates
fs_files = glob.glob(f"{data_dir}/{file_suffix}")
print(f"There are {len(fs_files)} feature selected files with the following paths:")
fs_files


# ### Generate dictionary with plate and data

# In[4]:


# Load feature data into a dictionary, keyed on plate name without the suffix
cp_dfs = {x.split("/")[-1].split("_sc")[0]: pd.read_parquet(x) for x in fs_files}

# Print out useful information about each dataset
print(cp_dfs.keys())
[cp_dfs[x].shape for x in cp_dfs]


# ### Fit UMAP for whole plates

# In[5]:


# Initialize a dictionary to store UMAP embeddings for each plate
umap_results = {}

# Fit UMAP features per dataset and save
for plate in cp_dfs:
    # Set output file for the UMAP
    output_umap_file = pathlib.Path(output_dir, f"UMAP_{plate}.parquet")

    # Check if the output file already exists
    if output_umap_file.exists():
        print(f"Skipping {output_umap_file.stem} as it already exists.")
        continue

    # Make sure to reinitialize UMAP instance per plate
    umap_fit = umap.UMAP(
        random_state=umap_random_seed, n_components=umap_n_components, n_jobs=1
    )

    # Set dataframe as the current plate
    cp_df = cp_dfs[plate]

    # Process cp_df to separate features and metadata
    cp_features = infer_cp_features(cp_df)
    meta_features = infer_cp_features(cp_df, metadata=True)

    # Subset to only failing + DMSO and healthy + DMSO for fitting (avoid flooding from the compound treated cells)
    fit_subset = cp_df[
        (
            (cp_df["Metadata_cell_type"] == "healthy")
            & (cp_df["Metadata_treatment"] == "DMSO")
        )
        | (
            (cp_df["Metadata_cell_type"] == "failing")
            & (cp_df["Metadata_treatment"] == "DMSO")
        )
    ]

    # Fit on the subset
    umap_fit.fit(fit_subset.loc[:, cp_features])

    # Transform entire dataset and convert to pandas DataFrame
    embeddings = pd.DataFrame(
        umap_fit.transform(cp_df.loc[:, cp_features]),
        columns=[f"UMAP{x}" for x in range(0, umap_n_components)],
    )
    print(f"{embeddings.shape}: {plate}")

    # Combine with metadata
    cp_umap_with_metadata_df = pd.concat(
        [cp_df.loc[:, meta_features], embeddings], axis=1
    )

    # Add treatment type column
    cp_umap_with_metadata_df["Metadata_treatment_type"] = np.select(
        [
            (cp_umap_with_metadata_df["Metadata_cell_type"] == "healthy")
            & (cp_umap_with_metadata_df["Metadata_treatment"] == "DMSO"),
            (cp_umap_with_metadata_df["Metadata_cell_type"] == "failing")
            & (cp_umap_with_metadata_df["Metadata_treatment"] == "DMSO"),
            (cp_umap_with_metadata_df["Metadata_cell_type"] == "failing")
            & (cp_umap_with_metadata_df["Metadata_treatment"] != "DMSO"),
        ],
        ["healthy + DMSO", "failing + DMSO", "failing + compound"],
        default="other",
    )
    # Add new metadata column to the list
    meta_features.append("Metadata_treatment_type")

    # Check and adjust dtypes dynamically
    for col in cp_umap_with_metadata_df.columns:
        if col in meta_features:
            # Try converting to numeric first (if possible), if not, keep as string
            try:
                cp_umap_with_metadata_df[col] = pd.to_numeric(
                    cp_umap_with_metadata_df[col], errors="raise", downcast="integer"
                )
            except ValueError:
                # If can't convert to numeric, keep as string
                cp_umap_with_metadata_df[col] = cp_umap_with_metadata_df[col].astype(
                    str
                )
        else:
            # For UMAP embeddings, ensure they're float
            cp_umap_with_metadata_df[col] = cp_umap_with_metadata_df[col].astype(float)

        # Update the 'Pathway' column for failing + DMSO and healthy + DMSO
        cp_umap_with_metadata_df["Metadata_Pathway"] = cp_umap_with_metadata_df.apply(
            lambda row: (
                "failing + DMSO"
                if row["Metadata_cell_type"] == "failing"
                and row["Metadata_treatment"] == "DMSO"
                else (
                    "healthy + DMSO"
                    if row["Metadata_cell_type"] == "healthy"
                    and row["Metadata_treatment"] == "DMSO"
                    else row["Metadata_Pathway"]
                )
            ),
            axis=1,
        )

    # Store the UMAP result in the dictionary
    umap_results[plate] = cp_umap_with_metadata_df

    # Generate output file, drop unnamed column, and save
    cp_umap_with_metadata_df.to_parquet(output_umap_file, index=False)


# ## Combine all plates together to generate UMAP embeddings

# In[6]:


# Get common features across all plates
common_columns = set.intersection(*[set(df.columns) for df in cp_dfs.values()])

# Use the first plate to identify metadata columns (those starting with 'Metadata_')
first_df = next(iter(cp_dfs.values()))
metadata_columns = [col for col in first_df.columns if col.startswith("Metadata_")]

# Get final columns to keep: intersection of common + metadata
final_columns = list(common_columns.union(metadata_columns))

# Subset each plate's dataframe to only those columns
cp_dfs_subset = {k: df[final_columns] for k, df in cp_dfs.items()}

# Combine all feature-selected plates with common features and metadata
combined_cp_df = pd.concat(cp_dfs_subset.values(), ignore_index=True)

# Verify shape and output
print(f"Combined shape: {combined_cp_df.shape}")
combined_cp_df.head()


# In[7]:


# Process combined_cp_df to separate features and metadata
combined_cp_features = infer_cp_features(combined_cp_df)
combined_meta_features = infer_cp_features(combined_cp_df, metadata=True)

# Subset to only failing + DMSO and healthy + DMSO for fitting
combined_fit_subset = combined_cp_df[
    (
        (combined_cp_df["Metadata_cell_type"] == "healthy")
        & (combined_cp_df["Metadata_treatment"] == "DMSO")
    )
    | (
        (combined_cp_df["Metadata_cell_type"] == "failing")
        & (combined_cp_df["Metadata_treatment"] == "DMSO")
    )
]

# Initialize UMAP instance
combined_umap_fit = umap.UMAP(
    random_state=umap_random_seed, n_components=umap_n_components, n_jobs=1
)

# Fit on the subset
combined_umap_fit.fit(combined_fit_subset.loc[:, combined_cp_features])

# Transform entire dataset and convert to pandas DataFrame
combined_embeddings = pd.DataFrame(
    combined_umap_fit.transform(combined_cp_df.loc[:, combined_cp_features]),
    columns=[f"UMAP{x}" for x in range(0, umap_n_components)],
)

# Combine with metadata
combined_cp_umap_with_metadata_df = pd.concat(
    [combined_cp_df.loc[:, combined_meta_features], combined_embeddings], axis=1
)

# Add treatment type column
combined_cp_umap_with_metadata_df["Metadata_treatment_type"] = np.select(
    [
        (combined_cp_umap_with_metadata_df["Metadata_cell_type"] == "healthy")
        & (combined_cp_umap_with_metadata_df["Metadata_treatment"] == "DMSO"),
        (combined_cp_umap_with_metadata_df["Metadata_cell_type"] == "failing")
        & (combined_cp_umap_with_metadata_df["Metadata_treatment"] == "DMSO"),
        (combined_cp_umap_with_metadata_df["Metadata_cell_type"] == "failing")
        & (combined_cp_umap_with_metadata_df["Metadata_treatment"] != "DMSO"),
    ],
    ["healthy + DMSO", "failing + DMSO", "failing + compound"],
    default="other",
)

# Update the 'Pathway' column for failing + DMSO and healthy + DMSO
combined_cp_umap_with_metadata_df["Metadata_Pathway"] = combined_cp_umap_with_metadata_df.apply(
    lambda row: (
        "failing + DMSO"
        if row["Metadata_cell_type"] == "failing"
        and row["Metadata_treatment"] == "DMSO"
        else (
            "healthy + DMSO"
            if row["Metadata_cell_type"] == "healthy"
            and row["Metadata_treatment"] == "DMSO"
            else row["Metadata_Pathway"]
        )
    ),
    axis=1,
)

# Save the combined UMAP embeddings to a file
combined_output_umap_file = pathlib.Path(output_dir, "UMAP_combined.parquet")
combined_cp_umap_with_metadata_df.to_parquet(combined_output_umap_file, index=False)

print(f"Combined UMAP embeddings saved to {combined_output_umap_file}")


# In[8]:


combined_cp_umap_with_metadata_df.head()
