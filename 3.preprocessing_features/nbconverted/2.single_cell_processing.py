#!/usr/bin/env python
# coding: utf-8

# # Process single cell profiles
# 
# NOTE: We are normalizing the plates for all samples as we only have three wells associated with the healthy controls, which is insufficient for normalization.

# ## Import libraries

# In[1]:


import os
import pathlib
import pprint

import pandas as pd

from pycytominer import annotate, normalize, feature_select


# ## Set paths and variables

# In[ ]:


# Optional: set `PLATEMAP_LAYOUT` env var to process only a single platemap (e.g. 'platemap_1')
platemap_to_process = os.environ.get("PLATEMAP_LAYOUT")
# platemap_to_process = "platemap_1"  # for testing only

# set base directory for where the SQLite files are located (should be local to repo)
base_dir = pathlib.Path("../2.cellprofiler_processing/cp_output/").resolve(strict=True)

# Decide what to process
if platemap_to_process:
    print(f"Processing only {platemap_to_process}")
    layouts = [platemap_to_process]
else:
    print("No specific layout set, processing all available platemaps")
    layouts = [p.name for p in base_dir.glob("platemap_*") if p.is_dir()]

pprint.pprint(layouts)


# In[ ]:


# Path to dir with cleaned data from single-cell QC
converted_dir = pathlib.Path(f"./data/{platemap_to_process}/cleaned_profiles/").resolve(
    strict=True
)

# output path for single-cell profiles
output_dir = pathlib.Path(f"./data/{platemap_to_process}/single_cell_profiles")
output_dir.mkdir(parents=True, exist_ok=True)

# Extract the plate names from the file name
plate_names = [
    "_".join(parts[:2]) if len(parts) >= 2 else parts[0]
    for parts in (file.stem.split("_") for file in converted_dir.glob("*.parquet"))
]


# path for platemap directory
platemap_dir = pathlib.Path("../metadata/updated_platemaps/")

# operations to perform for feature selection
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "blocklist",
    "drop_na_columns",
]


# ## Set dictionary with plates to process

# In[4]:


# Load the barcode_platemap file
barcode_platemap_df = pd.read_csv(
    pathlib.Path(f"{platemap_dir}/updated_barcode_platemap.csv").resolve()
)

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

# In[5]:


for plate, info in plate_info_dictionary.items():
    print(f"Performing pycytominer pipeline for {plate}")
    output_annotated_file = str(
        pathlib.Path(f"{output_dir}/{plate}_sc_annotated.parquet")
    )
    output_normalized_file = str(
        pathlib.Path(f"{output_dir}/{plate}_sc_normalized.parquet")
    )
    output_feature_select_file = str(
        pathlib.Path(f"{output_dir}/{plate}_sc_feature_selected.parquet")
    )

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
        output_type="parquet",
    )
    print(
        f"Annotation, normalization, and feature selection have been performed for {plate}"
    )


# In[6]:


# Check output file
test_df = pd.read_parquet(output_feature_select_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)


# In[7]:


# Check output file
test_df = pd.read_parquet(output_annotated_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)

