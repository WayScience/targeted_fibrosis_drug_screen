#!/usr/bin/env python

# # Run illumination correction on data
#
# Note: We load in the CellProfiler IC pipeline to use for this process.

# ## Import libraries

# In[1]:


import pathlib
import pprint
import sys

import requests

sys.path.append("../utils")
import cp_parallel

# ## Set paths and variables

# In[2]:


# set the run type for the parallelization
run_name = "illum_correction"

# set main output dir for all plates if it doesn't exist
output_dir = pathlib.Path("./Corrected_Images")
output_dir.mkdir(exist_ok=True)

# directory where images are located within folders (held in DropBox specific directory)
images_dir = pathlib.Path("../../Way Science Lab Dropbox/Jenna  Tomkinson/McKinseyLab_WayLab_CardiacFibroblasts/Compound Screen")

# list for plate names based on folders to use to create dictionary
plate_names = []
# iterate through 0.download_data and append plate names from folder names that contain image data from that plate
for file_path in images_dir.iterdir():
    plate_names.append(str(file_path.stem))

print("There are a total of", len(plate_names), "plates. The names of the plates are:")
for plate in plate_names:
    print(plate)


# ## Load in `illum.cppipe` file to process data

# In[3]:


# Define the GitHub raw file link (this link will get whatever file is in main)
github_url = "https://raw.githubusercontent.com/WayScience/cellpainting_predicts_cardiac_fibrosis/refs/heads/main/1.preprocessing_data/pipelines/illum.cppipe"

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


# ## Create dictionary with all plate data to run CellProfiler in parallel

# In[4]:


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


# ## Run CellProfiler Parallel
#
# Note: We do not run this code cell as we will run this process through the script.

# In[ ]:


cp_parallel.run_cellprofiler_parallel(
    plate_info_dictionary=plate_info_dictionary, run_name=run_name
)
