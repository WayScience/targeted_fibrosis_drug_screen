#!/usr/bin/env python
# coding: utf-8

# # CellProfiler segmentation and feature extraction
# 
# Note: We name this notebook as `analysis` to follow similar conventions set by CellProfiler.

# ## Import libraries

# In[1]:


import pathlib
import pprint
import requests

import sys

sys.path.append("../utils")
import cp_parallel


# ## Set paths and variables

# In[2]:


# set the run type for the parallelization
run_name = "analysis"

# path to output for SQLite database files per plate folder (create if does not already exist)
output_dir = pathlib.Path("./cp_output/")
output_dir.mkdir(exist_ok=True)

# Directory where all images are separated by folder per plate
images_dir = pathlib.Path("../1.illumination_correction/Corrected_Images").resolve(strict=True)

# list for plate names based on folders to use to create dictionary
plate_names = []

# iterate and append plate names from folders
for file_path in images_dir.iterdir():
    plate_names.append(str(file_path.stem))

print("There are a total of", len(plate_names), "plates. The names of the plates are:")
for plate in plate_names:
    print(plate)


# ## Load in CellProfiler analysis file to process data

# In[3]:


# Define the GitHub raw file link (this link will get whatever file is in main)
github_url = "https://raw.githubusercontent.com/WayScience/cellpainting_predicts_cardiac_fibrosis/refs/heads/main/2.cellprofiler_processing/pipeline/CFReT_project_CL.cppipe"

# Create the pipeline directory if it doesn't exist
pipeline_dir = pathlib.Path("pipeline")
pipeline_dir.mkdir(exist_ok=True)

# Download the file
response = requests.get(github_url)
response.raise_for_status()  # Raise an error for bad responses (4xx, 5xx)

# Save the file contents
file_path = pipeline_dir / github_url.split("/")[-1]
file_path.write_bytes(response.content)

print(f"File downloaded successfully to {file_path}")

# Create a variable to store the resolved path
path_to_pipeline = file_path.resolve(strict=True)


# ## Update the `Threshold correction factor` when segmenting nuclei and cells from the original file
# 
# When manually evaluating how the parameters are working with this dataset, I noticed some issues with the segmentation:
# 
# 1. Under-segmentation of whole cells
# 2. Segmentation of nuclei from empty images
# 
# Both of these segmentation parameters are updated from their original value to 0.5 to be more "lenient" per the documentation, but I have found it improves both of these issues. 
# Segmentation is never perfect but this makes a improvement from eye.
# 
# Note: The code below is hardcoded to change the parameters from a specific integer. This currently is what works, we will look to make this more generalizable in the future.

# In[4]:


# Read and modify the file
path_to_pipeline = file_path.resolve(strict=True)
with open(path_to_pipeline, 'r') as file:
    lines = file.readlines()

# Variables to keep track of where we are in the file
in_identify_secondary_objects = False
in_identify_primary_objects = False
threshold_correction_factor_found_secondary = False
threshold_correction_factor_found_primary = False

# Modify the content
with open(path_to_pipeline, 'w') as file:
    for line in lines:
        # Check if we're in the IdentifySecondaryObjects section
        if "IdentifySecondaryObjects" in line:
            in_identify_secondary_objects = True
            in_identify_primary_objects = False  # Ensure we're only in one section at a time
        
        # Check if we're in the IdentifyPrimaryObjects section
        if "IdentifyPrimaryObjects" in line:
            in_identify_primary_objects = True
            in_identify_secondary_objects = False  # Ensure we're only in one section at a time
        
        # If in the IdentifySecondaryObjects section and find the Threshold correction factor
        if in_identify_secondary_objects and "Threshold correction factor" in line:
            # Replace the value of the threshold correction factor (optimize segmentation for cells)
            line = line.replace("0.8", "0.5")
            threshold_correction_factor_found_secondary = True
        
        # If in the IdentifyPrimaryObjects section and find the Threshold correction factor
        if in_identify_primary_objects and "Threshold correction factor" in line:
            # Replace the value of the threshold correction factor (prevent non-existent nuclei from being segmented)
            line = line.replace("0.9", "0.3")
            threshold_correction_factor_found_primary = True
        
        # Write the line back to the file
        file.write(line)
        
        # Exit sections after processing the Threshold correction factor
        if in_identify_secondary_objects and threshold_correction_factor_found_secondary:
            in_identify_secondary_objects = False

        if in_identify_primary_objects and threshold_correction_factor_found_primary:
            in_identify_primary_objects = False

print(f"File updated successfully at {path_to_pipeline}")


# ## Create dictionary with all of the necessary paths to run CellProfiler analysis

# In[5]:


# create plate info dictionary with all parts of the CellProfiler CLI command to run in parallel
plate_info_dictionary = {
    name: {
        "path_to_images": pathlib.Path(list(images_dir.rglob(name))[0]).resolve(
            strict=True
        ),
        "path_to_output": pathlib.Path(f"{output_dir}/{name}"),
        "path_to_pipeline": path_to_pipeline,

    }
    for name in plate_names
}

# view the dictionary to assess that all info is added correctly
pprint.pprint(plate_info_dictionary, indent=4)


# ## Run CellProfiler analysis on all plates
# 
# **Note:** This code cell will not be run in this notebook due to the instability of jupyter notebooks compared to running as a python script. All CellProfiler SQLite outputs will have the same name but outputted into their respective plate folder (due to parallelization).

# In[ ]:


cp_parallel.run_cellprofiler_parallel(
    plate_info_dictionary=plate_info_dictionary, run_name=run_name
)

