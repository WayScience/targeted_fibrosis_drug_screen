#!/usr/bin/env python
# coding: utf-8

# # Convert SQLite outputs to parquet files with CytoTable

# ## Import libraries

# In[1]:


import pathlib
import pandas as pd
import pprint
import os

# cytotable will merge objects from SQLite file into single cells and save as parquet file
from cytotable import convert, presets

import logging

# Set the logging level to a higher level to avoid outputting unnecessary errors from config file in convert function
logging.getLogger().setLevel(logging.ERROR)


# ## Set paths and variables

# In[ ]:


# get the batch to process from environment variable
batch_to_process = os.environ.get("BATCH", "batch_1")
if batch_to_process is None:
    raise ValueError("Please set the BATCH environment variable before running this script.")

# base directory where batches are located
base_dir = pathlib.Path("../2.cellprofiler_processing/cp_output/").resolve(strict=True)

# Decide what to process
if batch_to_process:
    print(f"Processing {batch_to_process}")
    batch_dirs = [base_dir / batch_to_process]
else:
    print("No specific batch set, processing all available batches")
    batch_dirs = [p for p in base_dir.glob("batch_*") if p.is_dir()]

# Collect platemaps per batch
batch_layouts = {}
for batch_dir in batch_dirs:
    platemaps = [p.name for p in batch_dir.glob("platemap_*") if p.is_dir()]
    batch_layouts[batch_dir.name] = platemaps

pprint.pprint(batch_layouts)


# In[ ]:


# preset configurations based on typical CellProfiler outputs
preset = "cellprofiler_sqlite_pycytominer"

# update preset to include site metadata and cell counts
joins = presets.config[preset]["CONFIG_JOINS"].replace(
    "Image_Metadata_Well,",
    "Image_Metadata_Well, Image_Metadata_Site, Image_Count_Cells,",
)

# Add the PathName columns separately
joins = joins.replace(
    "COLUMNS('Image_FileName_.*'),",
    "COLUMNS('Image_FileName_.*'),\n COLUMNS('Image_PathName_.*'),",
)

# type of file output from cytotable
dest_datatype = "parquet"

# directory for processed data
output_base = pathlib.Path("data")
output_base.mkdir(exist_ok=True)

# confirm batch exists
batch_dir = base_dir / batch_to_process
if not batch_dir.is_dir():
    raise FileNotFoundError(f"Batch directory not found: {batch_dir}")

platemap_dirs = [p for p in batch_dir.glob("platemap_*") if p.is_dir()]

# collect platemap + plate info
batch_info = {}

for platemap_dir in platemap_dirs:
    # find all CARD-prefixed plate folders inside this platemap
    card_dirs = [p for p in platemap_dir.glob("CARD*") if p.is_dir()]

    plate_names = []
    for card_dir in card_dirs:
        # check if the CARD folder contains a SQLite file
        sqlite_files = list(card_dir.glob("*.sqlite"))
        if sqlite_files:
            plate_names.append(card_dir.name)  # use the CARD folder name as plate name

    batch_info[platemap_dir.name] = plate_names

    # create output directory for this platemap
    (output_base / batch_dir.name / platemap_dir.name).mkdir(
        parents=True, exist_ok=True
    )

# print summary
print(f"\nBatch: {batch_to_process} ({len(platemap_dirs)} platemaps)")
for platemap, plates in batch_info.items():
    print(f"  Platemap: {platemap} → {len(plates)} plates: {plates}")


# ## Convert SQLite to parquet files

# In[4]:


# loop through each platemap in the batch
for platemap_name, plate_names in batch_info.items():
    platemap_dir = batch_dir / platemap_name
    output_dir = output_base / batch_dir.name / platemap_name / "converted_profiles"
    output_dir.mkdir(parents=True, exist_ok=True)

    for plate_name in plate_names:
        card_dir = platemap_dir / plate_name
        sqlite_files = list(card_dir.glob("*.sqlite"))
        if not sqlite_files:
            continue  # skip if no sqlite found

        # assume one sqlite per CARD folder
        file_path = sqlite_files[0]
        output_path = output_dir / f"{plate_name}_converted.parquet"

        print(
            "Starting conversion with cytotable for plate:",
            plate_name,
            "from layout:",
            platemap_name,
            "from batch:",
            batch_to_process,
        )
        # Merge single cells and output as parquet file
        convert(
            source_path=str(file_path),
            dest_path=str(output_path),
            dest_datatype=dest_datatype,
            preset=preset,
            joins=joins,
            chunk_size=5000,
        )

print("All plates have been converted with cytotable!")


# # Load in converted profiles to update
# 
# We will rename some of the columns (e.g., location centroids and cell count per FOV) to include Metadata prefix.

# In[5]:


# List of columns to update with the "Metadata_" prefix
metadata_columns_to_update = [
    "Nuclei_Location_Center_X",
    "Nuclei_Location_Center_Y",
    "Cells_Location_Center_X",
    "Cells_Location_Center_Y",
    "Image_Count_Cells",
]

# loop through each platemap in the batch
for platemap_name, plate_names in batch_info.items():
    converted_dir = output_base / batch_dir.name / platemap_name / "converted_profiles"

    for plate_name in plate_names:
        file_path = converted_dir / f"{plate_name}_converted.parquet"
        if not file_path.is_file():
            print(f"Warning: file not found for plate {plate_name} in {platemap_name}")
            continue

        # Load the DataFrame from the Parquet file
        df = pd.read_parquet(file_path)

        # Drop rows where "Metadata_ImageNumber" is NaN
        df = df.dropna(subset=["Metadata_ImageNumber"])

        # Rearrange columns and add "Metadata_" prefix
        df = df[
            metadata_columns_to_update
            + [col for col in df.columns if col not in metadata_columns_to_update]
        ].rename(
            columns=lambda col: (
                "Metadata_" + col if col in metadata_columns_to_update else col
            )
        )

        # Save the processed DataFrame back to the same path
        df.to_parquet(file_path, index=False)
        print(
            f"Processed metadata columns for plate: {plate_name} in platemap: {platemap_name}"
        )

print("All converted profiles have been updated with Metadata columns!")


# ## Check output to confirm process worked
# 
# To confirm the number of single cells is correct, please use any database browser software to see if the number of rows in the "Per_Cells" compartment matches the number of rows in the data frame.

# In[6]:


# pick any platemap and plate to inspect (for example, the first ones)
first_platemap = next(iter(batch_info.keys()))
first_plate = batch_info[first_platemap][0]

converted_path = (
    output_base
    / batch_dir.name
    / first_platemap
    / "converted_profiles"
    / f"{first_plate}_converted.parquet"
)

# Load the selected converted parquet file
converted_df = pd.read_parquet(converted_path)

print(f"Loaded file: {converted_path}")
print(converted_df.shape)
converted_df.head()

