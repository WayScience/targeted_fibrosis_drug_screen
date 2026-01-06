#!/usr/bin/env python
# coding: utf-8

# # Perform single-cell quality control
# 
# In this notebook, we perform single-cell quality control using coSMicQC. We filter the single cells by identifying outliers with z-scores, and use either combinations of features or one feature for each condition. We use features from the AreaShape and Intensity modules to assess the quality of the segmented single-cells:
# 
# ### Assessing poor nuclei segmentation
# 
# Due to high confluence, sometimes nuclei overlap on top of each other, creating highly intense clusters within the Hoechst channel. To identify these nuclei, we use:
# 
# - **Nuclei Area:** This metric quantifies the number of pixels in a nucleus segmentation. 
# We detect nuclei that are abnormally large, which likely indicates poor nucleus segmentation where overlapping nuclei are merged into one segmentation. 
# - **Nuclei Intensity:** This metric quantifies the total intensity of all pixels in a nucleus segmentation. 
#   
# In combination with abnormally large nuclei, we detect nuclei that are also highly intense, likely indicating that this a group of overlapped nuclei.
# 
# `We utilize the same thresholds for this section as set in the [cellpainting_predicts_cardiac_fibroblasts](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis) repository.`
# 
# We decided in this notebook to also include finding nuclei with very **LOW** intensity which more than likely correlates with mis-segmented nuclei from the background. 
# These occur due to very cytotoxic compounds that kill all the cells, leaving empty FOVs where the segmentation parameters sometimes decides to segment nuclei from nothing.
# 
# `We determine our own threshold using the first batch of data (4 plates from layout one) that we will use for the rest of the plates.`
# 
# ### Assessing poor cell segmentation
# 
# Also due to high confluence, images with large, intense clusters of cells leads to errors in the segmentation algorithm that causes cells around the cluster to segmented incorrectly. 
# When this happens, a cell is segmented around the same segmentation as the nucleus, giving it the same area which is very small for a normal cardiac fibroblast cell. To detect poorly segmented cells, we use:
# 
# - **Cells area:** The cells Area metric quantifies the number of pixels in a cell segmentation.
# 
# `We update the threshold for this dataset from the original thresholds since we determined that there were too many good cells being removed. We update the threshold to be more loose.`

# In[ ]:


import os
import pathlib
import pprint

import matplotlib
matplotlib.use("Agg")  # Ensures no GUI windows are opened
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from cosmicqc import find_outliers


# ## Set paths and variables

# In[ ]:


# Optional: set `PLATEMAP_LAYOUT` env var to process only a single platemap (e.g. 'platemap_1')
platemap_to_process = os.environ.get("PLATEMAP_LAYOUT")
# platemap_to_process = "platemap_1"  # for testing only

# set base directory for where the SQLite files are located (should be local to repo)
base_dir = pathlib.Path("../2.cellprofiler_processing/cp_output/").resolve(strict=True)

# Decide what to process
if platemap_to_process:
    print(f"Processing only {platemap_to_process}")
    layouts = [platemap_to_process]
else:
    print("No specific layout set, processing all available platemaps")
    layouts = [p.name for p in base_dir.glob("platemap_*") if p.is_dir()]

pprint.pprint(layouts)


# In[3]:


# Directory with data
data_dir = pathlib.Path(f"./data/{platemap_to_process}/converted_profiles/").resolve(
    strict=True
)

# Directory to save cleaned data
cleaned_dir = pathlib.Path(f"./data/{platemap_to_process}/cleaned_profiles/")
cleaned_dir.mkdir(parents=True, exist_ok=True)

# Directory to save qc figures
qc_fig_dir = pathlib.Path(f"./qc_figures/{platemap_to_process}")
qc_fig_dir.mkdir(parents=True, exist_ok=True)

# Create an empty dictionary to store data frames for each plate
all_qc_data_frames = {}

# Metadata columns to include in output data frame
metadata_columns = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Metadata_Nuclei_Location_Center_X",
    "Metadata_Nuclei_Location_Center_Y",
]


# ## Load in plates to perform QC on

# In[4]:


# Create a dictionary to store DataFrames with plate names as keys
plate_data_dict = {}

# Iterate through all plates in the folder (assuming they are in Parquet format)
for plate_file in data_dir.glob("*_converted.parquet"):
    # Extract the plate name (stem) from the file path
    plate = plate_file.stem.replace("_converted", "")

    # Load the converted plate data
    plate_df = pd.read_parquet(plate_file)

    # Store the DataFrame in the dictionary with the plate name as the key
    plate_data_dict[plate] = {
        "converted_df": plate_df  # Store the loaded DataFrame under 'converted_df'
    }

    # Print the shape of the DataFrame for each plate
    print(f"Loaded plate: {plate}, Shape: {plate_df.shape}")


# ## Perform QC and removed failed single-cells across al plates

# In[5]:


# Process each plate in the dictionary and store QC results
for plate_name, plate_data in plate_data_dict.items():
    print(f"Processing QC for plate: {plate_name}")

    # Access the original DataFrame
    plate_df = plate_data["converted_df"]

    # OVER-SEGMENTED AND OVER-SATURATED NUCLEI
    ###########################################
    # Set outlier threshold for large nuclei and high intensity
    feature_thresholds_large_nuclei_high_int = {
        "Nuclei_AreaShape_Area": 2,
        "Nuclei_Intensity_IntegratedIntensity_DNA": 2,
    }

    # Find large nuclei and high intensity outliers
    large_nuclei_high_int_outliers = find_outliers(
        df=plate_df,
        metadata_columns=metadata_columns,
        feature_thresholds=feature_thresholds_large_nuclei_high_int,
    )

    # MIS-SEGMENTED BACKGROUND AS NUCLEI
    ###########################################
    # Set feature thresholds for low total intensity
    feature_thresholds_low_intensity = {
        "Nuclei_Intensity_IntegratedIntensity_DNA": -2,
    }

    # Find low intensity nuclei (most likely background)
    low_intensity_outliers = find_outliers(
        df=plate_df,
        metadata_columns=metadata_columns,
        feature_thresholds=feature_thresholds_low_intensity,
    )

    # UNDER-SEGMENTED CELLS
    ###########################################
    # Set feature thresholds for small cells
    feature_thresholds_small_cells = {
        "Cells_AreaShape_Area": -1.5,
    }

    # Find small cells outliers
    small_cells_outliers = find_outliers(
        df=plate_df,
        metadata_columns=metadata_columns,
        feature_thresholds=feature_thresholds_small_cells,
    )

    # Append the outliers data to the existing dictionary entry for this plate
    plate_data_dict[plate_name].update(
        {
            "low_intensity_outliers": low_intensity_outliers,
            "large_nuclei_high_int_outliers": large_nuclei_high_int_outliers,
            "small_cells_outliers": small_cells_outliers,
        }
    )

    # Find the outliers indices to determine failed single-cells
    outlier_indices = pd.concat(
        [large_nuclei_high_int_outliers, small_cells_outliers, low_intensity_outliers]
    ).index

    # Remove rows with outlier indices from the plate DataFrame
    plate_df_cleaned = plate_df.drop(outlier_indices)

    # Calculate the total percentage of nuclei that failed QC
    total_nuclei = plate_df.shape[0]
    total_failed = len(outlier_indices)
    percent_failed = (total_failed / total_nuclei) * 100 if total_nuclei > 0 else 0

    # Save cleaned data for this plate
    plate_cleaned_name = plate_df["Image_Metadata_Plate"].iloc[0]
    plate_df_cleaned.to_parquet(f"{cleaned_dir}/{plate_cleaned_name}_cleaned.parquet")

    # Verify the result and include the percentage of failed QC
    print(
        f"Cleaned data saved for plate {plate_cleaned_name}. Shape: {plate_df_cleaned.shape}. Total percent nuclei failed QC: {percent_failed:.2f}%"
    )


# In[6]:


# Create scatterplots per plate for nuclei area versus intensity labelled by QC condition failed or passed
for plate_name, plate_data in plate_data_dict.items():
    print(f"Creating QC plot for plate: {plate_name}")

    # Access the original dataframe
    plate_df = plate_data["converted_df"]  # Correctly assign the converted_df

    # Access the large nuclei high intensity outliers
    large_nuclei_high_int_outliers = plate_data["large_nuclei_high_int_outliers"]

    # Access the low intensity nuclei outliers
    low_intensity_outliers = plate_data["low_intensity_outliers"]

    # Set the default value to 'Single-cell passed QC'
    plate_df["Outlier_Status"] = "Single-cell passed QC"

    # Update the 'Outlier_Status' column based on the over-segmented nuclei outliers DataFrame using index
    plate_df.loc[
        plate_df.index.isin(large_nuclei_high_int_outliers.index), "Outlier_Status"
    ] = "Over-segmented nuclei"

    # Update the 'Outlier_Status' column based on the low intensity nuclei outliers DataFrame using index
    plate_df.loc[
        plate_df.index.isin(low_intensity_outliers.index), "Outlier_Status"
    ] = "Mis-segmented nuclei"

    # Create scatter plot
    plt.figure(figsize=(10, 6))
    plot = sns.scatterplot(
        data=plate_df,
        x="Nuclei_AreaShape_Area",
        y="Nuclei_Intensity_IntegratedIntensity_DNA",
        hue="Outlier_Status",
        palette={
            "Single-cell passed QC": "#006400",
            "Over-segmented nuclei": "#990090",
            "Mis-segmented nuclei": "#D55E00",
        },  # Specify color-blind friendly colors
        alpha=0.6,
    )

    # Add threshold lines
    plt.axvline(
        x=large_nuclei_high_int_outliers["Nuclei_AreaShape_Area"].min(),
        color="r",
        linestyle="--",
        label="Min. threshold for Nuclei Area",
    )
    plt.axhline(
        y=large_nuclei_high_int_outliers[
            "Nuclei_Intensity_IntegratedIntensity_DNA"
        ].min(),
        color="b",
        linestyle="--",
        label="Min. threshold for Nuclei Intensity",
    )

    # Add title and labels
    plt.title(f"Nuclei Area vs. Nuclei Integrated Intensity for {plate_name}")
    plt.xlabel("Nuclei Area")
    plt.ylabel("Nuclei Integrated Intensity (DNA)")
    plt.tight_layout()

    # Show the legend
    plt.legend(loc="upper left", bbox_to_anchor=(0, 1.0), prop={"size": 10})

    # Save the figure
    plt.savefig(qc_fig_dir / f"{plate_name}_nuclei_outliers.png", dpi=500)

    # Close the plot to prevent it from displaying
    plt.close()


# In[7]:


# Create density plots for each plate using outlier data for cells
for plate_name, plate_data in plate_data_dict.items():
    print(f"Creating density plot for plate: {plate_name}")

    # Access the original dataframe
    plate_df = plate_data["converted_df"]

    # Access the small cells outliers
    small_cells_outliers = plate_data["small_cells_outliers"]

    # Filter for rows that pass QC
    filtered_plate_df = plate_df[plate_df["Outlier_Status"] == "Single-cell passed QC"]

    # Create a density plot
    plt.figure(figsize=(10, 6))
    sns.kdeplot(x="Cells_AreaShape_Area", data=filtered_plate_df, fill=True)

    # Add threshold line
    plt.axvline(
        x=small_cells_outliers["Cells_AreaShape_Area"].max(),
        color="r",
        linestyle="--",
        label=f'Threshold for Outliers: < {small_cells_outliers["Cells_AreaShape_Area"].max()}',
    )

    # Set labels and title
    plt.ylabel("Count")
    plt.xlabel("Cells Area")
    plt.title(f"Distribution of Cells Area for {plate_name}")
    plt.legend()
    plt.tight_layout()

    # Save figure
    plt.savefig(qc_fig_dir / f"{plate_name}_cells_outliers.png", dpi=500)

    # Close the plot to prevent it from displaying
    plt.close()

