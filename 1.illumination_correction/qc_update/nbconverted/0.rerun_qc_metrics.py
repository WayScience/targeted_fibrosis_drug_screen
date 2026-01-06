#!/usr/bin/env python
# coding: utf-8

# ## Update QC metrics for blur to better catch poor quality images
# 
# Using the data from the `cellpainting_predicts_cardiac_fibrosis` repository, the QC thresholds will be updated to better detect poor quality images.
# 
# In the previous experiment, the thresholds were working well in those conditions.
# To apply to this experiment, we will rerun the QC to find optimal thresholds for blur for all channels for only one side of the distribution.
# More negative represents blur, so we only need one threshold to catch these conditions.
# More positive/close to 0 looks to represent empty images but we don't need to catch that condition.
# 
# As well, the QC for saturation has already been updated to have a universal threshold of 0.10 or 10% of pixels can be at the maximum value. 
# This is a stricter threshold that better accounts for the FOVs where cells are growing on top of each other.

# In[1]:


import pandas as pd


# In[2]:


# URL of the CSV file on GitHub
github_url = "https://raw.githubusercontent.com/WayScience/cellpainting_predicts_cardiac_fibrosis/main/1.preprocessing_data/qc_results/localhost231120090001/Image.csv"

# Load the CSV file into a pandas DataFrame
qc_df = pd.read_csv(github_url)


# In[3]:


# Calculate thresholds for each channel
channels = ["OrigActin", "OrigDNA", "OrigER", "OrigMito", "OrigPM"]
blur_thresholds = {}

for channel in channels:
    col = f"ImageQuality_PowerLogLogSlope_{channel}"
    # Check if the column exists in the DataFrame
    if col in qc_df.columns:
        # Calculate the 25th and 75th percentiles and the IQR
        Q1 = qc_df[col].quantile(0.25)
        Q3 = qc_df[col].quantile(0.75)
        IQR = Q3 - Q1
        # Calculate the blur threshold using IQR method (any value very negative)
        blur_thresholds[channel] = Q1 - 1.5 * IQR

# Display the calculated thresholds
print("Calculated blur thresholds:")
for channel, threshold in blur_thresholds.items():
    print(f"{channel}: {threshold}")

