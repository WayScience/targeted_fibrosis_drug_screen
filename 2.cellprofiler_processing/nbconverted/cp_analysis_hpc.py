#!/usr/bin/env python
# coding: utf-8

# # Run segmentation and feature extraction on data

# ## Import libraries

# In[1]:


import argparse
import pathlib
import pprint
import random

import sys

sys.path.append("../utils")
import cp_parallel

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False


# ## Set paths and variables

# In[ ]:


# Batch name to process (always contains batch_ prefix then #)
batch_name = "batch_1"

# directory where the corrected images are located within the folder
images_base_dir = pathlib.Path(
    f"../1.illumination_correction/Corrected_Images/{batch_name}"
).resolve(strict=True)

if not in_notebook:
    print("Running as script")
    # set up arg parser
    parser = argparse.ArgumentParser(
        description="CellProfiler segmentation and feature extraction"
    )

    parser.add_argument(
        "--image_dir",
        type=str,
        help="Path to the image directory to process corrected images",
    )

    args = parser.parse_args()
    images_dir = pathlib.Path(args.image_dir).resolve(strict=True)
else:
    print("Running in a notebook")

    # list all platemap layout folders in the batch
    platemap_folders = [p for p in images_base_dir.iterdir() if p.is_dir()]

    # gather all plates (folders starting with "CARD") inside all platemaps
    plate_folders = []
    for platemap in platemap_folders:
        plate_folders.extend(
            [p for p in platemap.iterdir() if p.is_dir() and p.name.startswith("CARD")]
        )

    if not plate_folders:
        raise ValueError(f"No plate folders starting with 'CARD' found in {batch_name}")

    # randomly select one plate
    images_dir = random.choice(plate_folders).resolve(strict=True)

    # extract plate name and platemap layout from folder structure
    plate_name = images_dir.name
    platemap_layout = images_dir.parent.name

    print(f"Processing plate: {plate_name}")
    print(f"Platemap layout: {platemap_layout}")

# set the run type for the parallelization
run_name = "cp_analysis"

# set path for CellProfiler pipeline
path_to_pipeline = pathlib.Path("./pipeline/analysis.cppipe").resolve(strict=True)

# set main output dir for all plates if it doesn't exist
output_dir = pathlib.Path("./sqlite_outputs")
output_dir.mkdir(exist_ok=True)


# ## Create dictionary with all plate data to run CellProfiler in parallel

# In[ ]:


# set path to the analysis pipeline
path_to_pipeline = pathlib.Path("./pipeline/analysis.cppipe").resolve(strict=True)

# set main output dir for all plates if it doesn't exist
output_dir = pathlib.Path("./cp_output").resolve(strict=False)
output_dir.mkdir(exist_ok=True)

# plate_info_dictionary for this run (single plate)
plate_info_dictionary = {}

# images_dir should be set from --image_dir in script mode or notebook selection
# create nested output dir: cp_output/batch/platemap_#/plate
plate_output_dir = output_dir / batch_name / images_dir.parent.name / images_dir.name
plate_output_dir.mkdir(parents=True, exist_ok=True)

# add info to dictionary
plate_info_dictionary[images_dir.name] = {
    "path_to_images": images_dir.resolve(strict=True),
    "path_to_output": plate_output_dir.resolve(strict=True),
    "path_to_pipeline": path_to_pipeline,
}

# print info for verification
print(f"Processing single plate: {images_dir.name}")
print(f"Platemap layout: {images_dir.parent.name}")
print(f"Output folder: {plate_output_dir.resolve(strict=True)}")

# view the dictionary
pprint.pprint(plate_info_dictionary, indent=4)


# ## Run CellProfiler Parallel
# 
# Note: We do not run this code cell as we will run this process through the script.

# In[ ]:


# if dictionary is not empty, run CellProfiler in parallel
if plate_info_dictionary:
    cp_parallel.run_cellprofiler_parallel(
        plate_info_dictionary=plate_info_dictionary,
        run_name=run_name,
        group_level="plate",
    )
else:
    print("No new plates to process. Exiting script.")

