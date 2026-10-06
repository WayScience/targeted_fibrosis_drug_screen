#!/usr/bin/env python
# coding: utf-8

# # Create CellProfiler LoadData CSV
# 
# Create a CellProfiler `LoadData` CSV from raw image paths. Each row is one image set, grouped by plate, well, and site, with one filename/path pair per Cell Painting channel.
# 
# Channel mapping from the main repository README:
# 
# | Image channel | CellProfiler image name |
# | --- | --- |
# | d0 | OrigActin |
# | d1 | OrigMito |
# | d2 | OrigGolgi |
# | d3 | OrigER |
# | d4 | OrigDNA |
# 

# In[1]:


from pathlib import Path
import re

import pandas as pd

import sys

sys.path.append("../../utils")
import loaddata_csv


# ## Set regular expressions for how to find metadata in file and folder names

# In[2]:


# Set the expected image file naming pattern for well site and channel
well_site_channel_pattern = re.compile(
    r"(?P<well>[A-H][0-9]{2})[_-]?f(?P<site>[0-9]{1,2})[_-]?(?P<channel>d[0-4])",
    flags=re.IGNORECASE,
)

# Set the expected pattern to find the plate name
plate_prefix_pattern = re.compile(
    r"^(?P<plate>.+?)[_-](?P<well>[A-H][0-9]{2})[_-]?f[0-9]{1,2}[_-]?d[0-4]",
    flags=re.IGNORECASE,
)

# Set expected folder pattern for plate naming to collect the name
plate_folder_pattern = re.compile(r"CARD-CelIns-CX7_[A-Za-z0-9_-]+")


# ## Configuration
# 
# Update `IMAGE_INPUTS` to point to the hit validation image directory.

# In[3]:


repo_root = Path("..").resolve()

# Path(s) to image directories
IMAGE_INPUTS = [
    Path(
        "/home/jenna/mnt/Way_McKinsey_Cardiac_Fibrosis/Compound_Screen/hit_validation_cell_painting/CARD-CelIns-CX7_260803130001"
    ),
]

# Set path to output LoadData CSVs
OUTPUT_DIR = Path("loaddata_csvs")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Used only if a plate cannot be inferred from the file path or filename.
DEFAULT_PLATE = None

# Set expected image extensions
IMAGE_EXTENSIONS = {".tif", ".tiff"}

# Set channel mapping to assign name in LoadData CSV
CHANNEL_MAP = {
    "d0": "OrigActin",
    "d1": "OrigMito",
    "d2": "OrigGolgi",
    "d3": "OrigER",
    "d4": "OrigDNA",
}

# Define expected columns for LoadData CSVs
LOADDATA_COLUMNS = [
    "Metadata_Plate",
    "Metadata_Well",
    "Metadata_Site",
    *[f"FileName_{image_name}" for image_name in CHANNEL_MAP.values()],
    *[f"PathName_{image_name}" for image_name in CHANNEL_MAP.values()],
]


# ## Build and save the LoadData CSV

# In[4]:


image_paths = loaddata_csv.collect_image_paths(
    IMAGE_INPUTS, image_extensions=IMAGE_EXTENSIONS
)
print(f"Found {len(image_paths):,} image files.")

loaddata_df = loaddata_csv.build_loaddata_csv(
    image_paths,
    channel_map=CHANNEL_MAP,
    loaddata_columns=LOADDATA_COLUMNS,
    default_plate=DEFAULT_PLATE,
    well_site_channel_pattern=well_site_channel_pattern,
    plate_prefix_pattern=plate_prefix_pattern,
    plate_folder_pattern=plate_folder_pattern,
)
print(f"Created {len(loaddata_df):,} image sets.")

missing_channels_df = loaddata_csv.summarize_missing_channels(
    loaddata_df, channel_map=CHANNEL_MAP
)
if not missing_channels_df.empty:
    print(
        f"Warning: {len(missing_channels_df):,} image sets are missing at least one channel."
    )
    display(missing_channels_df.head(20))

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

for plate, plate_loaddata_df in loaddata_df.groupby("Metadata_Plate", sort=True):
    output_csv = OUTPUT_DIR / f"loaddata_{plate}.csv"
    plate_loaddata_df.to_csv(output_csv, index=False)
    print(f"Saved {len(plate_loaddata_df):,} image sets to {output_csv}")

loaddata_df.head()

