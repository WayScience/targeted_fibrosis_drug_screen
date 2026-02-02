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

# In[2]:


# get the batch to process from environment variable
batch_to_process = os.environ.get("BATCH", "batch_1")
if batch_to_process is None:
    raise ValueError(
        "Please set the BATCH environment variable before running this script."
    )

# base directory where batches are located
base_dir = pathlib.Path("./data/").resolve(strict=True)

# Decide what to process
if batch_to_process:
    print(f"Processing {batch_to_process}")
    batch_dirs = [base_dir / batch_to_process]
else:
    print("No specific batch set, processing all available batches")
    batch_dirs = [p for p in base_dir.glob("batch_*") if p.is_dir()]


# In[3]:


# parameters
sc_fs_tag = "sc_feature_selected"
agg_tag = "aggregated_post_fs"

# batch to process
batch_dir = pathlib.Path(f"./data/{batch_to_process}").resolve(strict=True)

# find all layouts in this batch
layout_dirs = [p for p in batch_dir.glob("platemap_*") if p.is_dir()]

all_profiles_paths = []
unique_plate_names = set()

print(f"Batch: {batch_to_process} → {len(layout_dirs)} layouts found")

for layout_dir in layout_dirs:
    print(f"  Layout: {layout_dir.name}")

    # the folder containing parquet files
    sc_data_dir = layout_dir / "single_cell_profiles"
    if not sc_data_dir.is_dir():
        print(f"    ⚠️ No 'single_cell_profiles' folder found in {layout_dir.name}")
        continue

    parquet_files = list(sc_data_dir.glob(f"*{sc_fs_tag}.parquet"))
    if not parquet_files:
        print(f"    ⚠️ No feature-selected parquet files found in {layout_dir.name}")
        continue

    # extract plate names from file stems (assumes plate name is first part before _ or CARD prefix)
    plate_names = []
    for f in parquet_files:
        parts = f.stem.split("_")
        plate_name = "_".join(parts[:2])
        plate_names.append(plate_name)
        unique_plate_names.add(plate_name)

    print(
        f"    Found {len(parquet_files)} parquet files → unique plates: {set(plate_names)}"
    )

    # store all parquet files
    all_profiles_paths.extend(parquet_files)


# In[4]:


# Load the barcode_platemap file
barcode_platemap_path = pathlib.Path(
    "../metadata/updated_platemaps/updated_barcode_platemap.csv"
)
barcode_platemap_df = pd.read_csv(barcode_platemap_path)

plate_info_dictionary = {}

# Loop over each layout (platemap_#) in the batch
for layout_dir in batch_dir.glob("platemap_*"):
    sc_data_dir = layout_dir / "single_cell_profiles"

    # Find all feature-selected parquet files for this layout
    parquet_files = list(sc_data_dir.glob("*_sc_feature_selected.parquet"))

    for f in parquet_files:
        plate_name = "_".join(f.stem.split("_")[:2])  # extract unique plate id

        # Find corresponding platemap CSV
        platemap_row = barcode_platemap_df.loc[
            barcode_platemap_df["plate_barcode"] == plate_name
        ]
        if platemap_row.empty:
            raise ValueError(f"No platemap found for plate {plate_name}")

        platemap_path = (
            pathlib.Path("../metadata/updated_platemaps")
            / f"{platemap_row['platemap_file'].values[0]}.csv"
        ).resolve(strict=True)

        # Add to dictionary
        plate_info_dictionary[plate_name] = {
            "profile_path": f.resolve(strict=True),
            "platemap_path": platemap_path,
            "output_dir": batch_dir
            / layout_dir.name
            / "aggregated_profiles",  # keeps each plate in the correct layout folder
        }

# Confirm we have all plates
print(f"Number of plates to process: {len(plate_info_dictionary)}")
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

    # Set up output directory
    aggregated_dir_path = info["output_dir"]
    aggregated_dir_path.mkdir(parents=True, exist_ok=True)

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

