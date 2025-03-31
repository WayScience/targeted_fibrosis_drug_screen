#!/usr/bin/env python
# coding: utf-8

# # Extract model performance metrics
# 
# In this notebook, we extract metrics to evaluate performance such as:
# 
# 1. Precision-recall
# 2. Predicted probabilities

# ## Import libraries

# In[1]:


import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from joblib import load
from sklearn.metrics import precision_recall_curve

sys.path.append("../utils")
from training_utils import get_X_y_data


# # Set paths

# In[2]:


# Directory with the training and testing datasets per plate (or combined per batch)
data_dir = pathlib.Path("data")

# Directory with the trained models
model_dir = pathlib.Path("models")

# Directory with encoder
encoder_dir = pathlib.Path("encoder_results")

# Directory with the training indices
train_indices_dir = pathlib.Path("training_indices")

# Directory with the normalized datasets
normalized_data_path = pathlib.Path("../3.preprocessing_features/data/single_cell_profiles")

# Output directory the performance metrics
performance_metrics_dir = pathlib.Path("performance_metrics")
performance_metrics_dir.mkdir(exist_ok=True)


# ## Create dictionary with all relevant paths per plate to extract metrics

# In[3]:


# Get the list of encoder files
encoder_files = list(encoder_dir.glob("label_encoder*"))

# Extract plate names by removing the 'label_encoder' prefix
plate_names = [file.stem.replace("label_encoder_", "") for file in encoder_files]

# Create a nested dictionary with info per plate
plates_dict = {}
for plate in plate_names:
    plates_dict[plate] = {
        "training_data": data_dir / f"{plate}_train.parquet",
        "testing_data": data_dir / f"{plate}_test.parquet",
        "final_model": model_dir / f"{plate}_final_downsample.joblib",
        "shuffled_model": model_dir / f"{plate}_shuffled_downsample.joblib",
        "encoder_result": encoder_dir / f"label_encoder_{plate}.joblib",
        "training_indices": train_indices_dir / f"{plate}_training_data_indices.csv"
    }


# ## Extract metrics from the train and testing datasets applied to their respective plates

# In[4]:


# Label being predicted
label = "Metadata_cell_type"

# Initialize results list
pr_results = []

# Iterate through each plate
for plate, paths in plates_dict.items():
    # Load models
    final_model = load(paths["final_model"])
    shuffled_model = load(paths["shuffled_model"])
    
    # Load encoder
    label_encoder = load(paths["encoder_result"])
    
    # Load data
    train_df = pd.read_parquet(paths["training_data"])
    test_df = pd.read_parquet(paths["testing_data"])
    
    # Load training indices for filtering training data
    training_indices = pd.read_csv(paths["training_indices"])["Index"]
    train_df_filtered = train_df.loc[training_indices]
    
    # Define datasets
    datasets = {"train": train_df_filtered, "test": test_df}
    
    for dataset_name, dataset in datasets.items():
        # Retain Metadata_treatment before collecting X and y data
        metadata_treatment = dataset["Metadata_treatment"].values 
        # Extract features and labels
        X, y = get_X_y_data(df=dataset, label=label, shuffle=False)
        print(f"{plate}: {dataset_name} shape:", X.shape, )

        # Encode the labels
        y_encoded = label_encoder.transform(y)

        for model_name, model in [("final", final_model), ("shuffled", shuffled_model)]:
            # Predict probabilities for each cell
            y_scores = model.predict_proba(X)[:, 1]  # Assuming binary classification

            # Compute PR curve
            precision, recall, _ = precision_recall_curve(y_encoded, y_scores)
            
            # Append results, ensuring each row represents a single cell
            for actual, pred_prob, p, r, treatment in zip(y, y_scores, precision, recall, metadata_treatment):
                pr_results.append({
                    "model_type": model_name,
                    "dataset": dataset_name,
                    "plate_trained": plate,
                    "actual_label": actual,
                    "predicted_probability": pred_prob,
                    "precision": p,
                    "recall": r,
                    "Metadata_treatment": treatment
                })

# Convert results to a single dataframe
combined_train_test_df = pd.DataFrame(pr_results)

# Check output
print(combined_train_test_df.shape)
combined_train_test_df.head(2)


# In[5]:


# Iterate through each plate
cross_plate_results = {}

for plate, paths in plates_dict.items():
    if plate == "combined_batch1":
        continue  # Skip the combined_batch1 plate

    # Print current plate models being applied
    print(f"Plate models (final & shuffled) being applied is: {plate}")

    cross_plate_results[plate] = []  # Initialize an empty list for each plate

    # Load models
    final_model = load(paths["final_model"])
    shuffled_model = load(paths["shuffled_model"])

    # Load encoder
    label_encoder = load(paths["encoder_result"])

    # Get the feature names from the model
    model_features = final_model.feature_names_in_

    # Iterate over other plates using normalized data (holdout sets)
    for other_plate in plates_dict:
        if other_plate == plate or other_plate == "combined_batch1":
            continue  # Skip the same plate and combined_batch1

        other_normalized_path = f"{normalized_data_path}/{other_plate}_sc_normalized.parquet"
        other_combined_df = pd.read_parquet(other_normalized_path)

        # Filter out only the cells with Metadata_treatment as DMSO
        other_combined_df = other_combined_df[other_combined_df["Metadata_treatment"] == "DMSO"]

        # Drop rows with NaNs based on model features
        other_combined_df_drop_nans = other_combined_df.dropna(subset=model_features)

        # Get the metadata columns (those starting with 'Metadata_')
        metadata_columns = [col for col in other_combined_df_drop_nans.columns if col.startswith('Metadata_')]

        # Get model features that exist in the dataframe
        model_features_in_df = [col for col in model_features if col in other_combined_df_drop_nans.columns]

        # Filter the dataframe to keep only model features and metadata columns
        other_combined_df_filtered = other_combined_df_drop_nans[metadata_columns + model_features_in_df]

        # Retain Metadata_treatment before collecting X and y data
        metadata_treatment = other_combined_df_filtered["Metadata_treatment"].values 

        # Extract features and labels
        X, y = get_X_y_data(df=other_combined_df_filtered, label=label)

        # Assert that the columns in X match the features in the model
        assert all(col in model_features for col in X
                   ), "Features in the model do not match the columns in the dataset"

        # Encode the labels
        y_encoded = label_encoder.transform(y)

        for model_name, model in [("final", final_model), ("shuffled", shuffled_model)]:
            # Predict probabilities
            y_scores = model.predict_proba(X)[:, 1]

            # Compute PR curve
            precision, recall, thresholds = precision_recall_curve(y_encoded, y_scores)

            # Append the results to the list associated with the plate
            for actual, pred_prob, p, r, treatment in zip(y, y_scores, precision, recall, metadata_treatment):
                cross_plate_results[plate].append({
                    "model_type": model_name,
                    "dataset": f"holdout_{other_plate}",  # mark these other plates as holdout
                    "actual_label": actual,
                    "predicted_probability": pred_prob,
                    "precision": p,
                    "recall": r,
                    "plate_trained": plate,
                    "Metadata_treatment": treatment
                })
                
# Convert the results to DataFrame
cross_plate_results_dfs = {
    plate: pd.DataFrame(data) for plate, data in cross_plate_results.items()
}

# Print the shape of each DataFrame
for plate, df in cross_plate_results_dfs.items():
    print(f"Plate: {plate}, Shape: {df.shape}")

# Combine all the dataframes into one dataframe and add a 'plate' column to track origin
combined_holdout_df = pd.concat(cross_plate_results_dfs.values(), ignore_index=True)

# Check the output
print(combined_holdout_df.shape)
combined_holdout_df.head(2)


# In[6]:


# Label being predicted
label = "Metadata_cell_type"

# Initialize results dictionary
pr_results = {}

# Only process the combined_batch1 data
paths = plates_dict["combined_batch1"]

# Load models
final_model = load(paths["final_model"])
shuffled_model = load(paths["shuffled_model"])

# Load encoder
label_encoder = load(paths["encoder_result"])

# Load testing data
test_df = pd.read_parquet(paths["testing_data"])

# Split testing data by Metadata_Plate
test_groups = test_df.groupby("Metadata_Plate")

for plate, dataset in test_groups:
    pr_results[plate] = []  # Initialize an empty list for each testing plate

    # Retain Metadata_treatment before collecting X and y data
    metadata_treatment = dataset["Metadata_treatment"].values 

    # Extract features and labels
    X, y = get_X_y_data(df=dataset, label=label, shuffle=False)

    # Encode the labels
    y_encoded = label_encoder.transform(y)

    for model_name, model in [("final", final_model), ("shuffled", shuffled_model)]:
        # Predict probabilities
        y_scores = model.predict_proba(X)[:, 1]  # Assuming binary classification

        # Compute PR curve
        precision, recall, _ = precision_recall_curve(y_encoded, y_scores)

        # Store results as individual rows (one per precision-recall pair)
        for actual, pred_prob, p, r, treatment in zip(y, y_scores, precision, recall, metadata_treatment):
            pr_results[plate].append({
                "model_type": model_name,
                "dataset": f"test_{plate}",
                "actual_label": actual,
                "predicted_probability": pred_prob,
                "precision": p,
                "recall": r,
                "plate_trained": "combined_batch1",
                "Metadata_treatment": treatment
            })
            
# Convert results to dataframes
pr_results_dfs = {
    plate: pd.DataFrame(data) for plate, data in pr_results.items()
}

# Combine all results into one dataframe
combined_test_df = pd.concat(pr_results_dfs.values(), ignore_index=True)

# Check the output
print(combined_test_df.shape)
combined_test_df.head(2)


# In[7]:


# Combine all relevant DataFrames into one main DataFrame
performance_metrics_df = pd.concat(
    [combined_train_test_df, combined_holdout_df, combined_test_df],
    ignore_index=True
)

# Save the combined DataFrame as a parquet file
performance_metrics_df.to_parquet(f"{performance_metrics_dir}/performance_metrics.parquet", index=False)

# Check the shape of the final DataFrame
print(performance_metrics_df.shape)
performance_metrics_df.head()

