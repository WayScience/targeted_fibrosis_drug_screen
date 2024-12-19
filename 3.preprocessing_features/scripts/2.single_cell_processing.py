#!/usr/bin/env python
# coding: utf-8

# # Process single cell profiles

# ## Import libraries

# In[1]:


import pathlib
import pprint

import pandas as pd

from pycytominer import annotate, normalize, feature_select


# ## Set paths and variables

# In[2]:


# Path to dir with cleaned data from single-cell QC
converted_dir = pathlib.Path("./data/cleaned_profiles")

# output path for single-cell profiles 
output_dir = pathlib.Path("./data/single_cell_profiles")
output_dir.mkdir(parents=True, exist_ok=True)  

# Extract the plate names from the file name
plate_names = [file.stem.split("_")[0] for file in converted_dir.glob("*.parquet")]

# path for platemap directory
platemap_dir = pathlib.Path("../metadata/original_platemaps/")

# operations to perform for feature selection
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "blocklist",
    "drop_na_columns"
]


# ## Set dictionary with plates to process

# In[3]:


# Load the barcode_platemap file
barcode_platemap_df = pd.read_csv(pathlib.Path(f"{platemap_dir}/barcode_platemap.csv"))

# Create plate info dictionary
plate_info_dictionary = {
    name: {
        "profile_path": (converted_dir / f"{name}_cleaned.parquet").resolve(
            strict=True
        ),
        "platemap_path": (
            platemap_dir
            / f"{barcode_platemap_df.loc[barcode_platemap_df['plate_barcode'] == name, 'platemap_file'].values[0]}.csv"
        ).resolve(strict=True),
    }
    for name in plate_names
}

# View the dictionary to assess that all info is added correctly
pprint.pprint(plate_info_dictionary, indent=4)


# ## Process data with pycytominer

# In[4]:


for plate, info in plate_info_dictionary.items():
    print(f"Performing pycytominer pipeline for {plate}")
    output_annotated_file = str(pathlib.Path(f"{output_dir}/{plate}_sc_annotated.parquet"))
    output_normalized_file = str(pathlib.Path(f"{output_dir}/{plate}_sc_normalized.parquet"))
    output_feature_select_file = str(pathlib.Path(f"{output_dir}/{plate}_sc_feature_selected.parquet"))

    profile_df = pd.read_parquet(info["profile_path"])
    platemap_df = pd.read_csv(info["platemap_path"])

    print("Performing annotation for", plate, "...")
    # Step 1: Annotation
    annotate(
        profiles=profile_df,
        platemap=platemap_df,
        join_on=["Metadata_well_position", "Image_Metadata_Well"],
        output_file=output_annotated_file,
        output_type="parquet",
    )

    # Load the annotated parquet file to fix metadata columns names
    annotated_df = pd.read_parquet(output_annotated_file)

    # Rename columns using the rename() function
    column_name_mapping = {
        "Image_Metadata_Site": "Metadata_Site",
    }

    annotated_df.rename(columns=column_name_mapping, inplace=True)

    # Save the modified DataFrame back to the same location
    annotated_df.to_parquet(output_annotated_file, index=False)

    # Step 2: Normalization
    normalized_df = normalize(
        profiles=output_annotated_file,
        method="standardize",
        output_file=output_normalized_file,
        output_type="parquet",
    )
    
    print("Performing feature selection for", plate, "...")
    # Step 3: Feature selection
    feature_select(
        output_normalized_file,
        operation=feature_select_ops,
        na_cutoff=0,
        output_file=output_feature_select_file,
        output_type="parquet"
    )
    print(f"Annotation, normalization, and feature selection have been performed for {plate}")


# In[5]:


# Check output file
test_df = pd.read_parquet(output_feature_select_file)

print(test_df.shape)
test_df.head(2)

