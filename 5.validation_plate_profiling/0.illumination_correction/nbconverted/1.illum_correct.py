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

sys.path.append("../../utils")
import cp_parallel


# ## Set paths and variables

# ### Set the constants

# In[2]:


# set the run type for the parallelization
run_name = "illum_correction"


# ### Set up paths

# In[3]:


# set main output dir for all plates if it doesn't exist
output_dir = pathlib.Path("./illum_directory")
output_dir.mkdir(exist_ok=True)

# set base directory for where the images are located (WILL NEED TO CHANGE ON YOUR LOCAL MACHINE)
base_dir = pathlib.Path(
    "/home/jenna/mnt/Way_McKinsey_Cardiac_Fibrosis/Compound_Screen/hit_validation_cell_painting"
).resolve(strict=True)

# list for plate names based on folders found in the base directory
plate_names = sorted(folder.name for folder in base_dir.iterdir() if folder.is_dir())

print("Found", len(plate_names), "plates:")
for plate in plate_names:
    print(plate)


# ## Create dictionary with all plate data to run CellProfiler in parallel

# In[4]:


# set path to the illum pipeline
path_to_pipeline = pathlib.Path("./pipeline/illum.cppipe").resolve(strict=True)

# set path to loaddata csv files
loaddata_dir = pathlib.Path("./loaddata_csvs").resolve(strict=True)

# create plate info dictionary with all parts of the CellProfiler CLI command to run in parallel
plate_info_dictionary = {
    name: {
        "path_to_loaddata": list(loaddata_dir.rglob(f"loaddata_{name}.csv"))[0].resolve(
            strict=True
        ),
        "path_to_output": pathlib.Path(f"{output_dir}/{name}/"),
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

