#!/usr/bin/env python
# coding: utf-8

# # Run segmentation and feature extraction on data

# ## Import libraries

# In[1]:


import pathlib
import pprint

import sys
import os

sys.path.append("../utils")
import cp_parallel


# ## Set paths and variables

# In[ ]:


# Optional: set `PLATEMAP_LAYOUT` env var to process only a single platemap (e.g. 'platemap_1')
platemap_to_process = os.environ.get("PLATEMAP_LAYOUT")
# platemap_to_process = "platemap_1"  # for testing only

# set base directory for where the correct images are located (should be local to repo)
base_dir = pathlib.Path("../1.illumination_correction/Corrected_Images/").resolve(
    strict=True
)

# Decide what to process
if platemap_to_process:
    print(f"Processing only {platemap_to_process}")
    layouts = [platemap_to_process]
else:
    print("No specific layout set, processing all available platemaps")
    layouts = [p.name for p in base_dir.glob("platemap_*") if p.is_dir()]

pprint.pprint(layouts)


# ### Set the constants

# In[3]:


# set the run type for the parallelization
run_name = "cp_analysis"


# ### Set up paths

# In[4]:


# list for plate names
plate_names = []

# iterate through each platemap folder
for platemap_folder in base_dir.glob("platemap_*"):
    # if PLATEMAP_LAYOUT is set, only process that layout
    if platemap_to_process and platemap_folder.name != platemap_to_process:
        continue
    if platemap_folder.is_dir():
        for plate_folder in platemap_folder.iterdir():
            if plate_folder.is_dir() and plate_folder.name.startswith("CARD"):
                plate_names.append(plate_folder.name)

print("There are a total of", len(plate_names), "plates. The names of the plates are:")
for plate in plate_names:
    print(plate)


# ## Create dictionary with all plate data to run CellProfiler in parallel

# In[5]:


# set path to the analysis pipeline
path_to_pipeline = pathlib.Path("./pipeline/analysis.cppipe").resolve(strict=True)

# set main output dir for all plates if it doesn't exist
output_dir = pathlib.Path("./cp_output").resolve(strict=False)
output_dir.mkdir(exist_ok=True)

# create plate info dictionary
plate_info_dictionary = {}

for platemap_folder in base_dir.glob("platemap_*"):
    # if PLATEMAP_LAYOUT is set, only process that layout
    if platemap_to_process and platemap_folder.name != platemap_to_process:
        continue
    if platemap_folder.is_dir():
        for plate_folder in platemap_folder.iterdir():
            if plate_folder.is_dir() and plate_folder.name.startswith("CARD"):
                # create nested output dir: cp_output/platemap_#/plate
                plate_output_dir = output_dir / platemap_folder.name / plate_folder.name
                # create output dir and set dictionary if plate hasn't been processed
                if not plate_output_dir.exists() or not any(plate_output_dir.iterdir()):
                    plate_output_dir.mkdir(parents=True, exist_ok=True)

                    # add info to dictionary
                    plate_info_dictionary[plate_folder.name] = {
                        "path_to_images": plate_folder.resolve(strict=True),
                        "path_to_output": plate_output_dir.resolve(strict=True),
                        "path_to_pipeline": path_to_pipeline,
                    }
                else:
                    print(
                        f"{plate_output_dir} already exists and contains files, skipping creation and dictionary."
                    )

# view the dictionary to check
pprint.pprint(plate_info_dictionary, indent=4)


# ## Run CellProfiler Parallel
# 
# Note: We do not run this code cell as we will run this process through the script.

# In[ ]:


cp_parallel.run_cellprofiler_parallel(
    plate_info_dictionary=plate_info_dictionary, run_name=run_name, group_level="plate"
)

