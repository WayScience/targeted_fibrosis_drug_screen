#!/usr/bin/env python
# coding: utf-8

# # Convert SQLite outputs to parquet files with CytoTable

# ## Import libraries

# In[1]:


import pathlib
import pandas as pd
import pprint

# cytotable will merge objects from SQLite file into single cells and save as parquet file
from cytotable import convert, presets

import logging

# Set the logging level to a higher level to avoid outputting unnecessary errors from config file in convert function
logging.getLogger().setLevel(logging.ERROR)


# ## Set paths and variables

# In[2]:


# base directory where the validation plate's CellProfiler outputs are located
base_dir = pathlib.Path("../1.cellprofiler_processing/cp_output/").resolve(strict=True)

# find all plate folders that contain a CellProfiler SQLite output
plate_names = [
    p.name for p in base_dir.glob("CARD*") if p.is_dir() and list(p.glob("*.sqlite"))
]

print(f"There are a total of {len(plate_names)} plate(s) to convert:")
for plate in plate_names:
    print(plate)


# In[3]:


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
converted_dir = output_base / "converted_profiles"
converted_dir.mkdir(parents=True, exist_ok=True)

# create plate info dictionary with paths needed for conversion
plate_info_dictionary = {
    name: {
        "path_to_sqlite": pathlib.Path(
            next(iter((base_dir / name).glob("*.sqlite")))
        ).resolve(strict=True),
        "path_to_output": converted_dir / f"{name}_converted.parquet",
    }
    for name in plate_names
}

# view the dictionary to assess that all info is added correctly
pprint.pprint(plate_info_dictionary, indent=4)

converted_profiles_exist = bool(plate_info_dictionary) and all(
    info["path_to_output"].is_file() for info in plate_info_dictionary.values()
)

if converted_profiles_exist:
    print(
        "Converted profiles already exist for all plates. Conversion cells will be skipped."
    )


# ## Convert SQLite to parquet files

# In[4]:


if converted_profiles_exist:
    print(
        "Skipping SQLite to parquet conversion because converted profiles already exist."
    )
else:
    for plate_name, info in plate_info_dictionary.items():
        print("Starting conversion with cytotable for plate:", plate_name)
        # Merge single cells and output as parquet file
        convert(
            source_path=str(info["path_to_sqlite"]),
            dest_path=str(info["path_to_output"]),
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


if converted_profiles_exist:
    print("Skipping metadata updates because converted profiles already exist.")
else:
    # List of columns to update with the "Metadata_" prefix
    metadata_columns_to_update = [
        "Nuclei_Location_Center_X",
        "Nuclei_Location_Center_Y",
        "Cells_Location_Center_X",
        "Cells_Location_Center_Y",
        "Image_Count_Cells",
    ]

    for plate_name, info in plate_info_dictionary.items():
        file_path = info["path_to_output"]
        if not file_path.is_file():
            print(f"Warning: file not found for plate {plate_name}")
            continue

        # Load the DataFrame from the Parquet file
        df = pd.read_parquet(file_path)

        # Ensure Metadata_Plate contains only one unique value (occurs due to failure acquiring plate during run)
        if "Metadata_Plate" in df.columns and df["Metadata_Plate"].nunique() != 1:
            df["Metadata_Plate"] = plate_name

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
        print(f"Processed metadata columns for plate: {plate_name}")

    print("All converted profiles have been updated with Metadata columns!")


# ## Check output to confirm process worked
# 
# To confirm the number of single cells is correct, please use any database browser software to see if the number of rows in the "Per_Cells" compartment matches the number of rows in the data frame.

# In[6]:


# pick a plate to inspect (there is currently only one validation plate)
first_plate = next(iter(plate_info_dictionary.keys()))
converted_path = plate_info_dictionary[first_plate]["path_to_output"]

# Load the selected converted parquet file
converted_df = pd.read_parquet(converted_path)

print(f"Loaded file: {converted_path}")
print(converted_df.shape)
converted_df.head()

