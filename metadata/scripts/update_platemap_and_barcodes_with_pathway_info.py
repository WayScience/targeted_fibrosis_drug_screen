#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import pathlib
import pandas as pd

sys.path.append("../")
from utils import io_utils


# Setting up paths

# In[2]:


# setting tag for files that have been augmented with pathway data
tag = "with_pathways"
# setting up working directory and where the pathway_platemap file paths
original_platemaps_path = pathlib.Path("./original_platemaps").resolve(strict=True)
pathways_path = (original_platemaps_path / "pathways_platemap.csv").resolve(strict=True)

# settting original barcode path
barcode_path = (original_platemaps_path / "barcode_platemap.csv").resolve(strict=True)

# setting all the platemap paths
all_platemap_paths = list(
    original_platemaps_path.glob("Target_Selective_Library_Screen_*.csv")
)

# creating output directory
updated_platemaps_dir = pathlib.Path("./updated_platemaps")
updated_platemaps_dir.mkdir(exist_ok=True)

# update barcode output path
updated_barcode_path = (
    updated_platemaps_dir / "updated_barcode_platemap.csv"
).resolve()


# This process adds pathway metadata to experimental platemaps to provide more biological context for single-cell image-based profiles. Here's how it works:
# 
# 1. Pathway metadata, which includes treatment identifiers, well positions, and pathway information, is filtered for each specific plate.
# 2. This filtered metadata is merged with the platemap, linking treatments and well positions to their corresponding pathways.
# 3. The final augmented platemaps offer a clear view that connects experimental treatments to their associated pathways.

# In[3]:


# loading in the pathway_platemap metadata
pathway_meta_df = pd.read_csv(pathways_path)

# Iterate over all platemap files to augment them with pathway MOA information
for platemap_path in all_platemap_paths:
    # Generate an output filename for the augmented platemap
    # Example: Extract prefix and plate ID, then append a "with_moa" tag
    prefix = platemap_path.stem.rsplit("_", 2)[0]
    plate_id = platemap_path.stem.split("_", 4)[-1]
    outname = f"{prefix}_{plate_id}_{tag}.csv"

    # Load the current platemap into a DataFrame
    plate_meta_df = pd.read_csv(platemap_path)

    # Filter the pathway metadata to include only rows corresponding to the current plate
    # Select relevant columns: 'UCD ID' (treatment), 'Well' (well position), and 'Pathway' (moa)
    plate_pathway_meta_df = pathway_meta_df.loc[pathway_meta_df["Plate"] == plate_id][
        ["UCD ID", "Well", "Pathway"]
    ]

    # Merge the pathway metadata into the platemap using 'treatment' and 'well_position'
    # as keys. Perform a left join to retain all rows from the platemap
    merged_df = plate_meta_df.merge(
        plate_pathway_meta_df,
        left_on=["treatment", "well_position"],
        right_on=["UCD ID", "Well"],
        how="left",
    )

    # Drop the UCD ID and Well columns after the merge
    # contains identical information within the platemap
    merged_df = merged_df.drop(columns=["UCD ID", "Well"])

    # Ensure no rows or columns in the original platemap were modified during the merge
    # This checks that only the added MOA information changes the dataframe,
    # while the original content remains identical.
    assert merged_df[
        plate_meta_df.columns
    ].equals(
        plate_meta_df
    ), "The merged DataFrame's original rows or columns have been unexpectedly modified."

    # Save the augmented platemap with MOA information to the output directory
    merged_df.to_csv(updated_platemaps_dir / outname, index=False)


# In[4]:


# Load the original barcode metadata
barcodes = io_utils.load_barcodes(barcode_path=barcode_path)

# Initialize a list to store the updated barcode entries
updated_barcodes = []

# Iterate through each batch and its associated plates in the barcode metadata
for batch_id, batch_profiles in barcodes.items():
    # Generate the updated platemap name with the new tag
    updated_platemap_name = f"{list(batch_profiles.keys())[0]}_{tag}"

    # Iterate through all profiles associated with the current batch of plates
    for profile_names in batch_profiles.values():
        for profile_name in profile_names:
            # Append the updated profile-barcode mapping
            updated_barcodes.append([profile_name, updated_platemap_name])

# Convert the updated barcode list into a DataFrame
updated_barcodes_df = pd.DataFrame(
    updated_barcodes, columns=["plate_barcode", "platemap_file"]
)

# Save dataframe into the module specific result folder
updated_barcodes_df.to_csv(updated_barcode_path, index=False)

# Display new barcode file
updated_barcodes_df.head()

