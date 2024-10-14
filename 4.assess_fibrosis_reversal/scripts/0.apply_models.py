#!/usr/bin/env python
# coding: utf-8

# # Apply final and shuffled models to the data
# 
# We will be extracting the healthy and failing probabilities for each single-cell into a parquet file to use for figure generation in the next notebook.

# ## Import libraries

# In[1]:


import pathlib
import pprint
import requests

from io import BytesIO
import pandas as pd
import pyarrow.parquet as pq
from joblib import load

import sys

sys.path.append("../utils")
from training_utils import get_X_y_data


# ## Set function to download files from GitHub

# In[2]:


# Function to download and save a file from a GitHub URL
def download_file(github_url: str, save_dir: pathlib.Path) -> None:
    """Download a file from a GitHub raw URL and load it to a specified directory

    Args:
        github_url (str): string of a GitHub raw URL of the file to download
        save_dir (pathlib.Path): path to directory to save files to
    """
    response = requests.get(github_url)
    response.raise_for_status()  # Raise an error for bad responses (4xx, 5xx)

    file_name = github_url.split("/")[-1]
    file_path = save_dir / file_name
    file_path.write_bytes(response.content)

    print(f"File downloaded successfully to {file_path}")
    return file_path.resolve(strict=True)


# ## Load in the joblib files for the final and shuffled baseline models

# In[3]:


# Define the GitHub raw file links for final and shuffled models
final_model_url = "https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis/raw/refs/heads/main/5.machine_learning/models/log_reg_fs_plate_4_final_downsample.joblib"
shuffled_model_url = "https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis/raw/refs/heads/main/5.machine_learning/models/log_reg_fs_plate_4_shuffled_downsample.joblib"
label_encoder_url = "https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis/raw/refs/heads/main/5.machine_learning/encoder_results/label_encoder_log_reg_fs_plate_4.joblib"

# Create the models directory if it doesn't exist
models_dir = pathlib.Path("models")
models_dir.mkdir(exist_ok=True)

# Create the models directory if it doesn't exist
label_encoder_dir = pathlib.Path("label_encoder")
label_encoder_dir.mkdir(exist_ok=True)

# Download the final and shuffled models
final_model_path = download_file(final_model_url, models_dir)
shuffled_model_path = download_file(shuffled_model_url, models_dir)
encoder_path = download_file(label_encoder_url, label_encoder_dir)


# ## Get list of features for models to filter data

# In[4]:


# Define the GitHub raw link for the feature selected parquet file for the plate used with the models
github_url = "https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis/raw/refs/heads/main/3.process_cfret_features/data/single_cell_profiles/localhost231120090001_sc_feature_selected.parquet"

# Load the parquet file into memory (not locally)
response = requests.get(github_url)
response.raise_for_status()  # Ensure the request was successful

# Convert the response content into a BytesIO object for pandas to read
parquet_file = BytesIO(response.content)

# Use pyarrow to read only the schema (column names) from the Parquet file
schema = pq.read_schema(parquet_file)

# Extract columns that do not start with "Metadata"
model_features = [col for col in schema.names if not col.startswith("Metadata")]

# Output the number of columns to check if it matches correctly
print(len(model_features))


# ## Set paths and variables

# In[5]:


# Directory for probability data to be saved
prob_dir = pathlib.Path("./prob_data")
prob_dir.mkdir(exist_ok=True)

# Directory with normalized plate datasets
data_dir = pathlib.Path("../3.preprocessing_features/data/single_cell_profiles")

# Use rglob to search for files with the suffix *_sc_normalized.parquet
parquet_files = list(data_dir.rglob("*_sc_normalized.parquet"))

# Load all matching parquet files into a list of DataFrames
dfs = [pd.read_parquet(file) for file in parquet_files]

# Print the file paths and check the content if needed
for file in parquet_files:
    print(f"Loaded file: {file}")


# ## Filter the data for the model features and add to dictionary to process

# In[6]:


# Dictionary to store plate data
plate_data_dict = {}

# Loop through and process each parquet file
for parquet_file in parquet_files:
    # Extract the plate name from the file name
    plate_name = parquet_file.name.split("_")[0]
    print(f"Processing plate: {plate_name}")

    # Load the Parquet file
    df = pd.read_parquet(parquet_file)

    # Drop rows with NaN values in feature columns that the model uses
    df = df.dropna(subset=model_features)

    # Capitalize the cell type values to match the model
    df["Metadata_cell_type"] = df["Metadata_cell_type"].str.capitalize()

    # Extract metadata columns
    metadata_columns = [col for col in df.columns if col.startswith("Metadata_")]

    # Extract feature columns that don't start with "Metadata_"
    feature_columns = [col for col in df.columns if not col.startswith("Metadata_")]

    # Filter columns in the data frame to only include those in the model
    filtered_feature_columns = [col for col in feature_columns if col in model_features]

    # Filter the DataFrame to keep only the desired columns
    model_df = df[metadata_columns + filtered_feature_columns]

    # Store the processed DataFrame in the dictionary under the key "model_df"
    plate_data_dict[plate_name] = {"model_df": model_df}

    # Print info about the processed DataFrame
    print(
        f"Number of unique treatments in {plate_name}: {df['Metadata_treatment'].nunique()}"
    )
    print(f"Shape of the model DataFrame: {model_df.shape}")


# ## Extract final model predicted probabilities for each treatment

# In[7]:


# Create a list to store probability DataFrames from each loop iteration
prob_dfs = []

# Loop through each model in the models directory
for model_path in models_dir.iterdir():
    model_type = model_path.stem.split("_")[5]  # Get the model type

    # Process each plate's data from the plate_data_dict
    for plate_name, info in plate_data_dict.items():
        print(f"Extracting {model_type} probabilities from {plate_name} data...")

        # Load in model to apply to datasets
        model = load(model_path)

        # Load in label encoder
        le = load(encoder_path)

        # Get unique cell types and their corresponding encoded values
        unique_labels = le.classes_
        encoded_values = le.transform(unique_labels)

        # Create a dictionary mapping encoded values to original labels
        label_dict = dict(zip(encoded_values, unique_labels))

        # Load in the DataFrame associated with the plate
        data_df = info["model_df"].reset_index(drop=True)

        # Load in X data to get predicted probabilities
        X, _ = get_X_y_data(df=data_df, label="Metadata_cell_type")

        # Predict class probabilities for morphology feature data
        predicted_probs = model.predict_proba(X)

        # Storing probabilities in a pandas DataFrame
        prob_df = pd.DataFrame(predicted_probs, columns=model.classes_)

        # Update column names in prob_df using the dictionary and add suffix "_probas"
        prob_df.columns = [label_dict[col] + "_probas" for col in prob_df.columns]

        # Add a new column called predicted_label for each row
        prob_df["predicted_label"] = prob_df.apply(
            lambda row: row.idxmax()[:-7], axis=1
        )

        # Select metadata columns from the data
        metadata_columns = data_df.filter(like="Metadata_")

        # Combine metadata columns with predicted probabilities DataFrame based on index
        prob_df = prob_df.join(metadata_columns)

        # Add a new column for model_type
        prob_df["model_type"] = model_type

        # Append the DataFrame to the list instead of combining immediately
        prob_dfs.append(prob_df)

# Combine all DataFrames from the list into a single DataFrame after the loop
combined_prob_df = pd.concat(prob_dfs, ignore_index=True)

# Save combined probability data
combined_prob_df.to_csv(f"{prob_dir}/combined_batch_1_predicted_proba.csv", index=False)


# ## Display counts for correctly predicted cells across the cell types

# In[8]:


# Filter rows where Metadata_treatment is 'DMSO'
dmso_rows = combined_prob_df[combined_prob_df["Metadata_treatment"] == "DMSO"]

# Calculate counts and percentage of correct predictions for each Metadata_cell_type
result_counts = (
    dmso_rows.groupby("Metadata_cell_type")
    .apply(
        lambda x: pd.Series(
            {
                "correct_count": (
                    x["predicted_label"] == x["Metadata_cell_type"]
                ).sum(),
                "fail_count": (x["predicted_label"] != x["Metadata_cell_type"]).sum(),
                "total_count": x.shape[0],
            }
        )
    )
    .reset_index()
)

# Calculate the percentage of correct predictions
result_counts["percentage_correct"] = (
    result_counts["correct_count"] / result_counts["total_count"]
) * 100

# Print the number of correctly predicted cells with percentage for each Metadata_cell_type
for idx, row in result_counts.iterrows():
    print(
        f"{row['Metadata_cell_type']}: {row['correct_count']} correct predictions out of "
        f"{row['total_count']} cells ({row['percentage_correct']:.2f}%)"
    )

# Display the results
result_counts

