#!/usr/bin/env python
# coding: utf-8

# ## Extract coefficient values from all plates

# ## Import libraries

# In[1]:


from joblib import load
import pathlib
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt


# In[2]:


# set model directory
model_dir = pathlib.Path("../3a.train_individual_models/models").resolve(strict=True)

# load final models in to extract coefficients from
model_files = list(model_dir.glob("*final_downsample.joblib"))

# load in the original model (Circulation paper)
orig_circ_model = pd.read_csv(
    pathlib.Path(
        "/media/18tbdrive/1.Github_Repositories/cellpainting_predicts_cardiac_fibrosis/5.machine_learning/2.interpret_coefficients/coeff_data/all_coeffs.csv"
    ).resolve(strict=True)
)

# add a column for model name to match the newly trained models
orig_circ_model["Model"] = "original_circ_model"

# output directory for the coefficients
output_dir = pathlib.Path("coefficients")
output_dir.mkdir(parents=True, exist_ok=True)

# set directory for figure output
fig_dir = pathlib.Path("figures")
fig_dir.mkdir(parents=True, exist_ok=True)


# In[3]:


# loop through each model and extract the coefficients
for model_file in model_files:
    # load the model
    model = load(model_file)

    # get the model name
    parts = model_file.stem.split("_")
    model_name = "_".join(parts[:2]) if parts[0] == "combined" else parts[0]

    # get the coefficients
    coefs = pd.DataFrame(
        {"Feature": model.feature_names_in_, "Coefficient": model.coef_.flatten()}
    )

    # add model name to the coefficients
    coefs["Model"] = model_name

    # save the coefficients to a csv file
    coefs.to_csv(output_dir / f"{model_name}_coefficients.csv", index=False)


# ## Print top coefficients for combined model (batch1)

# In[4]:


# print off combined model coefficients in order from most important for healthy prediction
combined_coefs = pd.read_csv(output_dir / "combined_batch1_coefficients.csv")
print(f"Top 10 coefficients for healthy prediction for model {model_name}:")
combined_coefs.sort_values(by="Coefficient", ascending=False).head(10)


# In[5]:


# print off combined model coefficients in order from most important for failing prediction
print(f"Top 10 coefficients for failing prediction for model {model_name}:")
combined_coefs.sort_values(by="Coefficient", ascending=True).head(10)


# ## Generate distribution plots of the top actin feature for the batch1 data and the original model data (media only)

# ### Load in all plates from the batch

# In[6]:


# directory with normalized data
data_dir = pathlib.Path(
    "../3.preprocessing_features/data/single_cell_profiles"
).resolve(strict=True)

# get all of the files with normalized data and concat
data_files = list(data_dir.glob("*_sc_annotated.parquet"))
norm_combined_df = pd.concat(
    [pd.read_parquet(f) for f in data_files], ignore_index=True
)

# show dataframe
print(norm_combined_df.shape)
norm_combined_df.head()


# ### Generate the plot for the plates in the batch

# In[7]:


# Increase font sizes globally
plt.rcParams.update(
    {
        "font.size": 16,
        "axes.titlesize": 18,
        "axes.labelsize": 16,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
    }
)

palette_dict = {"healthy": "#004400", "failing": "#a0004b"}

# Filter the dataframe for DMSO treatment
filtered_df = norm_combined_df[norm_combined_df["Metadata_treatment"] == "DMSO"]

# Select a feature to plot
feature = "Cells_Intensity_IntegratedIntensity_Actin"

# Create the KDE plot
plt.figure(figsize=(15, 5))
sns.kdeplot(
    data=filtered_df[filtered_df["Metadata_cell_type"] == "healthy"],
    x=feature,
    label="Healthy",
    fill=True,
    alpha=0.5,
    color=palette_dict["healthy"],
)
sns.kdeplot(
    data=filtered_df[filtered_df["Metadata_cell_type"] == "failing"],
    x=feature,
    label="Failing",
    fill=True,
    alpha=0.5,
    color=palette_dict["failing"],
)

# Add labels and legend
plt.title("Actin intensity distribution across batch 1 DMSO cells")
plt.xlabel(feature)
plt.ylabel("Density")
plt.legend()

# Save the plot
plt.savefig(
    f"{fig_dir}/batch1_actin_intensity_distribution.png", dpi=500, bbox_inches="tight"
)
plt.show()


# ### Load in the plate data from the original model

# In[8]:


# directory with normalized data
orig_data_dir = pathlib.Path(
    "/media/18tbdrive/1.Github_Repositories/cellpainting_predicts_cardiac_fibrosis/3.process_cfret_features/data/single_cell_profiles/"
)

# get all of the files with normalized data and concat
orig_annot_df = pd.read_parquet(
    pathlib.Path(orig_data_dir / "localhost231120090001_sc_annotated.parquet")
)

# show dataframe
print(orig_annot_df.shape)
orig_annot_df.head()


# ### Generate plot with the original data but filter to only the same hearts as the screening data

# In[9]:


# Increase font sizes globally
plt.rcParams.update(
    {
        "font.size": 16,
        "axes.titlesize": 18,
        "axes.labelsize": 16,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
    }
)

palette_dict = {"healthy": "#004400", "failing": "#a0004b"}

# Filter the dataframe only media treatment and hearts 7 and 19
filtered_orig_df = orig_annot_df[orig_annot_df["Metadata_treatment"] != "DMSO"]
filtered_orig_df = filtered_orig_df[
    filtered_orig_df["Metadata_heart_number"].isin([7, 19])
]

# Select a feature to plot
feature = "Cells_Intensity_IntegratedIntensity_Actin"

# Create the KDE plot
plt.figure(figsize=(15, 5))
sns.kdeplot(
    data=orig_annot_df[orig_annot_df["Metadata_cell_type"] == "Healthy"],
    x=feature,
    label="Healthy",
    fill=True,
    alpha=0.5,
    color=palette_dict["healthy"],
)
sns.kdeplot(
    data=orig_annot_df[orig_annot_df["Metadata_cell_type"] == "Failing"],
    x=feature,
    label="Failing",
    fill=True,
    alpha=0.5,
    color=palette_dict["failing"],
)

# Add labels and legend
plt.title("Actin intensity distribution for the original model data (media only)")
plt.xlabel(feature)
plt.ylabel("Density")
plt.legend()
plt.savefig(
    f"{fig_dir}/orig_model_actin_intensity_distribution.png",
    dpi=500,
    bbox_inches="tight",
)
plt.show()


# ## Perform a outer merge of the circulation model and combined model coefficients
# 
# To avoid losing important features from each model, we are using outer merge so all unique features from each model are included. Given the each model did not use all the same features, NaNs will be added where there is not a match. We change these NaNs to 0's so we can still compare all the features across models.

# In[10]:


# Perform an outer merge of the circulation model and the combined model
merged_coefs = pd.merge(
    orig_circ_model,
    combined_coefs,
    on="Feature",
    how="outer",
    suffixes=("_orig_circ_model", "_combined_batch1_model"),
)

# Drop model column as it doesn't add any information
merged_coefs = merged_coefs.drop(
    columns=[col for col in merged_coefs.columns if "Model" in col]
)

# Fill NaN values with 0
merged_coefs.fillna(0, inplace=True)

# Save the merged coefficients to a CSV file
merged_coefs.to_csv(
    output_dir / "merged_coefficients_circ_model_combined_batch1.csv", index=False
)

# Display the merged dataframe
print(merged_coefs.shape)
merged_coefs.head()

