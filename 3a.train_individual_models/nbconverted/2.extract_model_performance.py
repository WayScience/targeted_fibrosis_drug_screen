#!/usr/bin/env python

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

import pandas as pd
from joblib import load
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve
from sklearn.preprocessing import LabelEncoder

sys.path.append("../utils")
from training_utils import get_X_y_data

# ## Helper function to collect precision-recall results and predicted probabilities

# In[2]:


def get_pr_curve_results(
    model: LogisticRegression, df: pd.DataFrame, label: str, label_encoder: LabelEncoder
) -> pd.DataFrame:
    """Collect the precision-recall curve results from a model and dataset.

    Args:
        model (LogisticRegression): loaded in logistic regression model to collect results from
        df (pd.DataFrame): dataframe containing the data to apply model to
        label (str): label with the class being predicted
        label_encoder (LabelEncoder): encoder to transform the labels to integers

    Returns:
        pd.DataFrame: dataframe with the PR curve results for that data and model
    """
    # Get X and y data for the model
    X, y = get_X_y_data(df=df, label=label, shuffle=False)

    assert all(col in model.feature_names_in_ for col in X), \
        "Features in the model do not match the columns in the dataset"

    y_encoded = label_encoder.transform(y)
    y_scores = model.predict_proba(X)[:, 1]

    precision, recall, _ = precision_recall_curve(y_encoded, y_scores)

    return pd.DataFrame({
        "precision": precision[:-1],  # remove last to align with thresholds
        "recall": recall[:-1],
    })


# In[3]:


def get_predicted_probabilities(
    model: LogisticRegression, df: pd.DataFrame, label: str, label_encoder: LabelEncoder
) -> pd.DataFrame:
    """Collect predicted probabilities per single-cell from the model and dataset.

    Args:
        model (LogisticRegression): loaded in logistic regression model to collect results from
        df (pd.DataFrame): dataframe containing the data to apply model to
        label (str): label with the class being predicted
        label_encoder (LabelEncoder): encoder to transform the labels to integers

    Returns:
        pd.DataFrame: dataframe with the predicted probabilities per single-cell
    """
    metadata_treatment = df["Metadata_treatment"].values
    X, y = get_X_y_data(df=df, label=label, shuffle=False)

    assert all(col in model.feature_names_in_ for col in X), \
        "Features in the model do not match the columns in the dataset"

    # y_encoded = label_encoder.transform(y)
    y_scores = model.predict_proba(X)[:, 1]

    return pd.DataFrame({
        "actual_label": y,
        "predicted_probability": y_scores,
        "Metadata_treatment": metadata_treatment
    })


# # Set paths

# In[4]:


# Directory with the training and testing datasets per plate (or combined per batch)
data_dir = pathlib.Path("data")

# Directory with the trained models
model_dir = pathlib.Path("models")

# Directory with encoder
encoder_dir = pathlib.Path("encoder_results")

# Directory with the training indices
train_indices_dir = pathlib.Path("training_indices")

# Directory with the normalized datasets
normalized_data_path = pathlib.Path(
    "../3.preprocessing_features/data/single_cell_profiles"
)

# Output directory the performance metrics
performance_metrics_dir = pathlib.Path("performance_metrics")
performance_metrics_dir.mkdir(exist_ok=True)

# Label being predicted
label = "Metadata_cell_type"


# ## Create dictionary with all relevant paths per plate to extract metrics

# In[5]:


# Get the list of encoder files
encoder_dir = pathlib.Path("./encoder_results")

# Extract plate names from model filenames
plate_names = {
    f.stem.replace("_final_downsample", "")
    for f in model_dir.glob("*_final_downsample.joblib")
}

# Create a nested dictionary with info per plate
plates_dict = {}
for plate in plate_names:
    plates_dict[plate] = {
        "training_data": data_dir / f"{plate}_train.parquet",
        "testing_data": data_dir / f"{plate}_test.parquet",
        "final_model": model_dir / f"{plate}_final_downsample.joblib",
        "shuffled_model": model_dir / f"{plate}_shuffled_downsample.joblib",
        "encoder_result": encoder_dir / "label_encoder_global.joblib",
        "training_indices": train_indices_dir / f"{plate}_training_data_indices.csv",
    }


# ## Extract metrics from the train and testing datasets applied to their respective plates

# In[6]:


# Initialize results list
test_train_pr_results = []
test_train_probability_results = []

# Run through each plate and get the PR results for training and testing data only
for plate, paths in plates_dict.items():
    # Load the models and data
    final_model = load(paths["final_model"])
    shuffled_model = load(paths["shuffled_model"])
    label_encoder = load(paths["encoder_result"])
    train_df = pd.read_parquet(paths["training_data"])
    test_df = pd.read_parquet(paths["testing_data"])

    # Filter the training data to only include the indices used in training the model
    training_indices = pd.read_csv(paths["training_indices"])["Index"]
    train_df_filtered = train_df.loc[training_indices]

    # Set dictionary with the training and testing data
    datasets = {"train": train_df_filtered, "test": test_df}

    # Loop through both datasets and models
    for dataset_name, dataset in datasets.items():
        for model_name, model in [("final", final_model), ("shuffled", shuffled_model)]:
            # Get per-sample predicted probabilities
            prob_df = get_predicted_probabilities(
                model=model, df=dataset, label=label, label_encoder=label_encoder
            )
            prob_df["model_type"] = model_name
            prob_df["dataset"] = dataset_name
            prob_df["plate_trained"] = plate
            test_train_probability_results.append(prob_df)

            # Get PR curve results (global)
            pr_df = get_pr_curve_results(
                model=model, df=dataset, label=label, label_encoder=label_encoder
            )
            pr_df["model_type"] = model_name
            pr_df["dataset"] = dataset_name
            pr_df["plate_trained"] = plate
            test_train_pr_results.append(pr_df)

            print(f"{model_name.upper()} | {plate} | {dataset_name} → Done")

# Combine all results into one dataframe
train_test_all_models_pr_results_df = pd.concat(test_train_pr_results, ignore_index=True)
train_test_all_models_probabilities_df = pd.concat(test_train_probability_results, ignore_index=True)

# Check output
print(train_test_all_models_probabilities_df.shape)
train_test_all_models_probabilities_df.head(2)


# ## Apply the individual plate models to the three other plates in the batch (considered holdout)

# In[7]:


# Initialize lists to store PR curve results and predicted probabilities
holdout_plate_individual_models_results = []
holdout_plate_individual_models_probability_results = []

# Iterate through each plate
for plate, paths in plates_dict.items():
    if plate == "combined_batch1":
        continue  # Skip the combined_batch1 plate

    # Print current plate models being applied
    print(f"Plate models being applied is: {plate}")

    # Load models
    final_model = load(paths["final_model"])
    shuffled_model = load(paths["shuffled_model"])

    # Load encoder
    label_encoder = load(paths["encoder_result"])

    # Get the feature names from the model
    model_features = final_model.feature_names_in_

    # Iterate over other plates using normalized data (holdout sets)
    for holdout_plate in plates_dict:
        # Skip processing for the model plate or combined_batch1
        if holdout_plate == plate or holdout_plate == "combined_batch1":
            continue

        # Load in the normalized data for the holdout plate
        holdout_normalized_path = (
            f"{normalized_data_path}/{holdout_plate}_sc_normalized.parquet"
        )
        holdout_norm_df = pd.read_parquet(holdout_normalized_path)

        # Filter out only the cells with Metadata_treatment as DMSO
        holdout_norm_df = holdout_norm_df[
            holdout_norm_df["Metadata_treatment"] == "DMSO"
        ]

        # Drop rows with NaNs based on model features
        holdout_norm_df = holdout_norm_df.dropna(subset=model_features)

        # Get the metadata columns (those starting with 'Metadata_')
        metadata_columns = [
            col
            for col in holdout_norm_df.columns
            if col.startswith("Metadata_")
        ]

        # Get model features from the dataframe
        model_features_in_df = [
            col for col in model_features if col in holdout_norm_df.columns
        ]

        # Filter the dataframe to keep only model features and metadata columns
        holdout_norm_df = holdout_norm_df[
            metadata_columns + model_features_in_df
        ]

        # Loop through models (final and shuffled)
        for model_name, model in [("final", final_model), ("shuffled", shuffled_model)]:
            # Get predicted probabilities (per-sample) using the function
            prob_df = get_predicted_probabilities(
                model=model, df=holdout_norm_df, label=label, label_encoder=label_encoder
            )
            prob_df["model_type"] = model_name
            prob_df["dataset"] = f"holdout_{holdout_plate}"  # Mark as holdout plate
            prob_df["plate_trained"] = plate
            holdout_plate_individual_models_probability_results.append(prob_df)

            # Get PR curve results using the function
            pr_df = get_pr_curve_results(
                model=model, df=holdout_norm_df, label=label, label_encoder=label_encoder
            )

            # Add context columns
            pr_df["model_type"] = model_name
            pr_df["dataset"] = f"holdout_{holdout_plate}"  # Mark as holdout plate
            pr_df["plate_trained"] = plate

            # Append to results list for the plate
            holdout_plate_individual_models_results.append(pr_df)

            # Print the shape of pr_df per holdout plate
            print(f"{model_name} applied to {holdout_plate} shape: {pr_df.shape}")

# Combine all results into one dataframe for both PR results and probabilities
holdout_plate_individual_models_pr_df = pd.concat(holdout_plate_individual_models_results, ignore_index=True)
holdout_plate_individual_models_probability_df = pd.concat(holdout_plate_individual_models_probability_results, ignore_index=True)

print(f"Final combined Probability DataFrame shape: {holdout_plate_individual_models_probability_df.shape}")
holdout_plate_individual_models_probability_df.head(2)


# ## Apply the combined model to the testing dataset split by plate to evaluate performance

# In[8]:


# Initialize lists to store PR curve results and predicted probabilities
combined_model_split_test_pr_results = []
combined_model_split_test_probability_results = []

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
    for model_name, model in [("final", final_model), ("shuffled", shuffled_model)]:
        # Get PR curve results
        pr_df = get_pr_curve_results(
            model=model, df=dataset, label=label, label_encoder=label_encoder
        )
        # Add model and dataset context
        pr_df["model_type"] = model_name
        pr_df["dataset"] = f"test_{plate}"
        pr_df["plate_trained"] = "combined_batch1"

        # Append PR curve results
        combined_model_split_test_pr_results.append(pr_df)

        # Get predicted probabilities (per-sample) using the function
        prob_df = get_predicted_probabilities(
            model=model, df=dataset, label=label, label_encoder=label_encoder
        )
        # Add model and dataset context
        prob_df["model_type"] = model_name
        prob_df["dataset"] = f"test_{plate}"
        prob_df["plate_trained"] = "combined_batch1"

        # Append predicted probabilities
        combined_model_split_test_probability_results.append(prob_df)

        # Print the shape of pr_df per test plate split
        print(f"{model_name} applied to {plate} shape: {pr_df.shape}")

# Combine all results into one dataframe for PR curve results and probabilities
combined_model_split_test_pr_results_df = pd.concat(combined_model_split_test_pr_results, ignore_index=True)
combined_model_split_test_probability_results_df = pd.concat(combined_model_split_test_probability_results, ignore_index=True)

# Check the output
print(f"Final combined Probability DataFrame shape: {combined_model_split_test_probability_results_df.shape}")
combined_model_split_test_probability_results_df.head(2)


# In[9]:


# Combine all relevant DataFrames into one main DataFrame
pr_results_df = pd.concat(
    [train_test_all_models_pr_results_df, holdout_plate_individual_models_pr_df, combined_model_split_test_pr_results_df], ignore_index=True
)
probabilities_df = pd.concat(
    [train_test_all_models_probabilities_df, holdout_plate_individual_models_probability_df, combined_model_split_test_probability_results_df], ignore_index=True
)


# Save the combined DataFrame as a parquet file
pr_results_df.to_parquet(
    f"{performance_metrics_dir}/batch1_pr_curve_results.parquet", index=False
)
probabilities_df.to_parquet(
    f"{performance_metrics_dir}/batch1_probabilities_DMSO_results.parquet", index=False
)

# Check the shape of the final DataFrame
print(probabilities_df.shape)
probabilities_df.head()
