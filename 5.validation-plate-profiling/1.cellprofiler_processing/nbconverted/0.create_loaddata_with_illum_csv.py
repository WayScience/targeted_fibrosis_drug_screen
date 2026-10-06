#!/usr/bin/env python
# coding: utf-8

# # Create CellProfiler LoadData CSV with Illumination Functions
# 
# Create a CellProfiler `LoadData` CSV from raw image paths including paths to the illumination correction functions per channel. Each row is one image set, grouped by plate, well, and site, with one filename/path pair per Cell Painting channel and paths to illumination functions (npy).
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

# In[ ]:


from pathlib import Path

import sys

sys.path.append("../../utils")
import loaddata_csv


# ## Set paths and variables

# In[ ]:


# Set path to find original loaddata csv files
orig_loaddata_csvs = Path("../0.illumination_correction/loaddata_csvs").resolve()

# Set output directory for loaddata csvs
output_dir = Path("loaddata_csvs")
output_dir.mkdir(exist_ok=True)

# Set path for illum_directory
illum_directory = Path("../0.illumination_correction/illum_directory")

# Set channel mapping to assign name in LoadData CSV
CHANNEL_MAP = {
    "d0": "OrigActin",
    "d1": "OrigMito",
    "d2": "OrigGolgi",
    "d3": "OrigER",
    "d4": "OrigDNA",
}


# ## Update and save the LoadData CSV with illumination function columns

# In[3]:


updated_dfs = loaddata_csv.add_illumination_columns_to_loaddata_csvs(
    loaddata_csv_paths=list(orig_loaddata_csvs.rglob("loaddata_*.csv")),
    illum_directory=illum_directory,
    channel_map=CHANNEL_MAP,
    output_dir=output_dir,
)

# Print one plate example
display(next(iter(updated_dfs.values())).head())

