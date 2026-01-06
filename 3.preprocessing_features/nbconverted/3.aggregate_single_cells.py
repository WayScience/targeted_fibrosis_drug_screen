#!/usr/bin/env python
# coding: utf-8

# ## Aggregating feature-selected single cells
# 
# In this notebook, we process single-cell feature-selected profiles to generate compound-level aggregated profiles for each plate using the pycytominer. 
# The single-cell profiles are grouped by treatment (Metadata_treatment) and are saved as Parquet files in the aggregated_profiles directory. 
# These aggregated profiles provide concise and interpretable data for downstream analysis at the compound level.

# ## Import libraries

# In[1]:


import pathlib
import pprint
import random
import pandas as pd
import os

from pycytominer import aggregate, annotate


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


# In[3]:


# parameters
sc_fs_tag = "sc_feature_selected"
agg_tag = "aggregated_post_fs"

# setting up paths
data_dir = pathlib.Path("./data").resolve(strict=True)
sc_data_dir = pathlib.Path(
    f"../3.preprocessing_features/data/{platemap_to_process}/single_cell_profiles"
).resolve(strict=True)

# setting metadata paths
metadata_dir = pathlib.Path("../metadata/updated_platemaps").resolve(strict=True)
updated_barcode_path = (metadata_dir / "updated_barcode_platemap.csv").resolve(
    strict=True
)
all_profiles_paths = list(sc_data_dir.glob("*sc_feature_selected.parquet"))

# output files paths
aggregated_dir_path = (
    data_dir / f"{platemap_to_process}/aggregated_profiles"
).resolve()
aggregated_dir_path.mkdir(exist_ok=True)

# Extract the plate names from the file name
plate_names = [
    "_".join(parts[:2]) if len(parts) >= 2 else parts[0]
    for parts in (file.stem.split("_") for file in all_profiles_paths)
]
print(plate_names)


# In[4]:


# Load the barcode_platemap file
barcode_platemap_df = pd.read_csv(updated_barcode_path)

# Create plate info dictionary
plate_info_dictionary = {
    name: {
        "profile_path": (sc_data_dir / f"{name}_sc_feature_selected.parquet").resolve(
            strict=True
        ),
        "platemap_path": (
            metadata_dir
            / f"{barcode_platemap_df.loc[barcode_platemap_df['plate_barcode'] == name, 'platemap_file'].values[0]}.csv"
        ).resolve(strict=True),
    }
    for name in plate_names
}

# View the dictionary to assess that all info is added correctly
pprint.pprint(plate_info_dictionary, indent=4)


# Next, we use the aggregation functionality provided by pycytominer to consolidate single-cell profiles into well-level summaries for each plate. This step groups the data by a specified metadata column and computes aggregate statistics by using the median.

# In[5]:


# Iterate over all profile file paths to process and aggregate data
for plate, info in plate_info_dictionary.items():
    # Load the current plate's feature selected profile data
    plate_path = info["profile_path"]

    # Load the single-cell profile data from the current Parquet file into a DataFrame
    profile_df = pd.read_parquet(plate_path)

    # Move the Well column to the first position
    profile_df = profile_df[
        ["Metadata_Well"]
        + [col for col in profile_df.columns if col != "Metadata_Well"]
    ]

    # Apply the aggregation function using pycytominer to aggregate at the well level
    agg_df = aggregate(
        profile_df,
        strata=["Metadata_Well"],
    )

    # Load the platemap data
    platemap_df = pd.read_csv(info["platemap_path"])

    # Perform annotation to make sure that all metadata is added back
    annotate(
        profiles=agg_df,
        platemap=platemap_df,
        join_on=["Metadata_well_position", "Metadata_Well"],
        output_type="parquet",
        output_file=(aggregated_dir_path / f"{plate}_{agg_tag}.parquet").resolve(),
    )


# In[6]:


# Get a list of Parquet files in the directory
parquet_files = list(aggregated_dir_path.glob("*.parquet"))

# Check if there are any files in the directory
if parquet_files:
    # Randomly select a file
    random_file = random.choice(parquet_files)

    # Load the randomly selected file
    test_df = pd.read_parquet(random_file)
else:
    print(f"No Parquet files found in directory: {aggregated_dir_path}")

# Display information
print(f"Randomly selected file: {random_file.relative_to(pathlib.Path.cwd())}")
print(test_df.shape)
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)

