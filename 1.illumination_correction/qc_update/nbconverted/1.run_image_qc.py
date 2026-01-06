#!/usr/bin/env python
# coding: utf-8

# # Run image QC on data
# 
# Note: We load in the CellProfiler QC pipeline to use for this process.

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
run_name = "image_qc"


# ### Set up paths

# In[3]:


# set base directory for where the images are located
base_dir = pathlib.Path(
    "/media/18tbdrive/CFReT_screening_data/compound_screen/"
).resolve(strict=True)

# list for plate names
plate_names = []

# iterate through each platemap folder
for platemap_folder in base_dir.glob("platemap_*"):
    if platemap_folder.is_dir():
        # iterate through the layout folders inside platemap
        for layout_folder in platemap_folder.iterdir():
            if layout_folder.is_dir():
                # now iterate through CARD plates inside the layout folder
                for plate_folder in layout_folder.iterdir():
                    if plate_folder.is_dir() and plate_folder.name.startswith("CARD"):
                        plate_names.append(plate_folder.name)

print("There are a total of", len(plate_names), "plates. The names of the plates are:")
for plate in plate_names:
    print(plate)


# ## Create dictionary with all plate data to run CellProfiler in parallel

# In[4]:


# set path to the illum pipeline
path_to_pipeline = pathlib.Path("./pipeline/image_qc.cppipe").resolve(strict=True)

# set main output dir for all plates if it doesn't exist
output_dir = pathlib.Path("./qc_results")
output_dir.mkdir(exist_ok=True)

# create plate info dictionary
plate_info_dictionary = {}

for platemap_folder in base_dir.glob("platemap_*"):
    if platemap_folder.is_dir():
        for alias_folder in platemap_folder.iterdir():
            if alias_folder.is_dir():
                for plate_folder in alias_folder.iterdir():
                    if plate_folder.is_dir() and plate_folder.name.startswith("CARD"):
                        # create nested output dir: qc_results/platemap_#/plate
                        plate_output_dir = (
                            output_dir / platemap_folder.name / plate_folder.name
                        )

                        # create output dir and set dictionary if plate hasn't been processed
                        if not plate_output_dir.exists() or not any(
                            plate_output_dir.iterdir()
                        ):
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


# if dictionary is not empty, run CellProfiler in parallel
if plate_info_dictionary:
    cp_parallel.run_cellprofiler_parallel(
        plate_info_dictionary=plate_info_dictionary, run_name=run_name, group_level="plate"
    )
else:
    print("No new plates to process. Exiting script.")

