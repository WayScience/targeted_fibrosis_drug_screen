#!/usr/bin/env python
# coding: utf-8

# # Run illumination correction on data
# 
# Note: We load in the CellProfiler IC pipeline to use for this process.

# ## Import libraries

# In[1]:


import pathlib
import pprint

import sys

sys.path.append("../utils")
import cp_parallel


# ## Set paths and variables

# ### Set the constants

# In[2]:


# set the run type for the parallelization
run_name = "illum_correction"

# set up the batch name for the plate(s) being processed
batch_name = "batch_2"


# ### Set up paths

# In[3]:


# set main output dir for all plates if it doesn't exist
output_dir = pathlib.Path("./illum_directory")
output_dir.mkdir(exist_ok=True)

# ensure the directory for the batch_name exists within the illum_directory
batch_output_dir = output_dir / batch_name
batch_output_dir.mkdir(exist_ok=True)

# set base directory for where the images are located (WILL NEED TO CHANGE ON YOUR LOCAL MACHINE)
base_dir = pathlib.Path(
    "../../Way Science Lab Dropbox/Jenna  Tomkinson/McKinseyLab_WayLab_CardiacFibroblasts/Compound Screen/"
).resolve(strict=True)

# folder where images are located within folders
images_dir = pathlib.Path(f"{base_dir}/Plate 2").resolve(strict=True)

# list for plate names based on folders to use to create dictionary
plate_names = []
# iterate through 0.download_data and append plate names from folder names that contain image data from that plate
for file_path in images_dir.iterdir():
    plate_names.append(str(file_path.stem))

print("There are a total of", len(plate_names), "plates. The names of the plates are:")
for plate in plate_names:
    print(plate)


# ## Create dictionary with all plate data to run CellProfiler in parallel

# In[4]:


# set path to the illum pipeline
path_to_pipeline = pathlib.Path("./pipeline/illum.cppipe").resolve(strict=True)

# create plate info dictionary with all parts of the CellProfiler CLI command to run in parallel
plate_info_dictionary = {
    name: {
        "path_to_images": pathlib.Path(list(images_dir.rglob(name))[0]).resolve(
            strict=True
        ),
        "path_to_output": pathlib.Path(f"{output_dir}/{batch_name}/{name}"),
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
    plate_info_dictionary=plate_info_dictionary, run_name=run_name, group_level="plate"
)

