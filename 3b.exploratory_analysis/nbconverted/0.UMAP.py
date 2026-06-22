#!/usr/bin/env python
# coding: utf-8

# # Extract UMAP embeddings for each plate of CellProfiler features
# 
# NOTE: We are using the feature selected profiles per plate.

# In[1]:


import glob
import pathlib
import pandas as pd
import numpy as np
import umap
import random

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
base_data_dir = pathlib.Path("../3.preprocessing_features/data")

# Discover all batch_#/platemap_# combinations
batch_platemap_dirs = sorted(base_data_dir.glob("batch_*/platemap_*/single_cell_profiles"))

print(f"Found {len(batch_platemap_dirs)} platemap directories:")
for d in batch_platemap_dirs:
    print(f"  {d.relative_to(base_data_dir)}")


# ### Extract platemap information and organize by batch/platemap
# 

# In[4]:


# Create a dictionary to organize data by platemap_#
platemap_data = {}

# Load feature data for each platemap
for data_dir in batch_platemap_dirs:
    # Extract platemap identifier from path (e.g., "batch_1/platemap_1" -> "platemap_1")
    batch_part = data_dir.parent.parent.name  # batch_#
    platemap_part = data_dir.parent.name      # platemap_#
    platemap_key = f"{batch_part}/{platemap_part}"
    
    # Find all feature-selected files in this platemap directory
    file_suffix = "*sc_feature_selected.parquet"
    fs_files = sorted(data_dir.glob(file_suffix))
    
    # Load all files for this platemap
    platemap_dfs = {}
    for f in fs_files:
        plate_name = f.stem.split("_sc")[0]  # Extract plate name from filename
        platemap_dfs[plate_name] = pd.read_parquet(f)
    
    platemap_data[platemap_key] = platemap_dfs
    print(f"{platemap_key}: {len(platemap_dfs)} plates loaded")

print(f"\nTotal platemaps: {len(platemap_data)}")


# ### Fit UMAP for whole plates
# 
# Remove single plate, only color by the 9 top compounds, fix facet labels

# In[5]:


# Process each platemap
for platemap_key, cp_dfs in platemap_data.items():
    print(f"\n{'='*60}")
    print(f"Processing {platemap_key}")
    print(f"{'='*60}")
    
    # Create output directory for this platemap
    platemap_output_dir = pathlib.Path(output_dir, platemap_key)
    platemap_output_dir.mkdir(parents=True, exist_ok=True)
    
    # ===== COMBINED UMAP FOR THIS PLATEMAP =====
    print(f"\nGenerating combined UMAP for {platemap_key}...")
    
    combined_output_umap_file = platemap_output_dir / f"UMAP_combined_{platemap_key.replace('/', '_')}.parquet"
    
    if combined_output_umap_file.exists():
        print(f"  Combined UMAP already exists, skipping.")
    else:
        # Get common features across all plates in this platemap
        common_columns = set.intersection(*[set(df.columns) for df in cp_dfs.values()])
        
        # Use the first plate to identify metadata columns
        first_df = next(iter(cp_dfs.values()))
        metadata_columns = [col for col in first_df.columns if col.startswith("Metadata_")]
        
        # Get final columns to keep
        final_columns = list(common_columns.union(metadata_columns))
        
        # Subset each plate's dataframe to only those columns
        cp_dfs_subset = {k: df[final_columns] for k, df in cp_dfs.items()}
        
        # Combine all plates in this platemap
        combined_cp_df = pd.concat(cp_dfs_subset.values(), ignore_index=True)
        
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
        
        # Initialize and fit UMAP instance
        combined_umap_fit = umap.UMAP(
            random_state=umap_random_seed, n_components=umap_n_components, n_jobs=1
        )
        combined_umap_fit.fit(combined_fit_subset.loc[:, combined_cp_features])
        
        # Transform entire dataset
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
        
        # Update the 'Pathway' column
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
        
        # Save combined UMAP
        combined_cp_umap_with_metadata_df.to_parquet(combined_output_umap_file, index=False)
        print(f"  Saved combined UMAP: {combined_cp_umap_with_metadata_df.shape}")

print(f"\n{'='*60}")
print("All platemaps processed successfully!")
print(f"{'='*60}")

