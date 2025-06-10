#!/usr/bin/env python

# # Apply models to the cells from wells treated with a small-molecule compound and extract probabilities

# ## Import libraries

# In[1]:


import pathlib
import pprint
import sys

import pandas as pd
from joblib import load

sys.path.append("../utils")
from training_utils import get_X_y_data

# ## Set paths and variables

# In[2]:


# path to models
models_dir = pathlib.Path("./models")

# directory with normalized datasets
norm_profiles_dir = pathlib.Path(
    "../3.preprocessing_features/data/single_cell_profiles"
)

# create a list of the plate names for gathering paths
plate_names = {
    (
        "_".join(f.stem.split("_")[:2])
        if f.stem.startswith("combined")
        else f.stem.split("_")[0]
    )
    for f in models_dir.glob("*.joblib")
}

for name in sorted(plate_names):
    print(name)


# ## Create dictionary with paths per model

# In[3]:


# Create a dictionary to store paths for each plate
model_paths_dir = {}

# Iterate through each plate name
for plate in plate_names:
    if plate.startswith("combined"):
        norm_matches = list(
            norm_profiles_dir.glob("*sc_normalized.parquet")
        )  # Collect paths to each plate
    else:
        norm_matches = list(
            norm_profiles_dir.glob(f"{plate}*sc_normalized.parquet")
        )  # Only collect path for that plate

    # Create a dictionary for each plate with paths to final and shuffled models
    model_paths_dir[plate] = {
        "final_model": list(models_dir.glob(f"{plate}*final_downsample.joblib")),
        "shuffled_model": list(models_dir.glob(f"{plate}*shuffled_downsample.joblib")),
        "norm_path": norm_matches,
    }

pprint.pprint(model_paths_dir)


# ## Extract probabilities for the compound treated cells

# In[4]:


# Initialize a list to store all probability DataFrames
all_predictions = []

for plate, items in model_paths_dir.items():
    print(f"Processing model: {plate}")
    # Check if the model paths exist
    if not items["final_model"] or not items["shuffled_model"]:
        print(f"Model paths for {plate} do not exist.")
        continue

    # Load the models
    final_model = load(items["final_model"][0])
    shuffled_model = load(items["shuffled_model"][0])

    # Load the normalized data
    if plate == "combined_batch1":
        # Concatenate all plate paths for the combined case
        norm_data_frames = [pd.read_parquet(path) for path in items["norm_path"]]
        norm_data = pd.concat(norm_data_frames, ignore_index=True)

        # Save the cell counts per well per plate w/ metadata
        cell_counts = (
            norm_data.groupby(["Metadata_Plate", "Metadata_Well"])
            .size()
            .reset_index(name="cell_count")
        )
        cell_counts = cell_counts.merge(
            norm_data.drop_duplicates(subset=["Metadata_Plate", "Metadata_Well"])[
                [
                    col
                    for col in norm_data.columns
                    if col.startswith("Metadata_")
                    and not any(
                        x in col for x in ["Image", "Location", "Parent", "Object"]
                    )
                ]
            ],
            on=["Metadata_Plate", "Metadata_Well"],
            how="left",
        )
        cell_counts.to_csv("./performance_metrics/batch1_cell_counts.csv", index=False)
    else:
        # Load the single parquet file for non-combined plates (same plate as model)
        norm_path = items["norm_path"][0]
        norm_data = pd.read_parquet(norm_path)

    # Filter norm_data for only the rows with compound treatment (remove DMSO)
    norm_data = norm_data[norm_data["Metadata_treatment"] != "DMSO"]

    # Get the feature columns used in the model to filter the data
    model_columns = (
        final_model.feature_names_in_
    )  # Feature names are same in shuffled model

    # Filter the data to only include the model columns including metadata columns
    meta_cols = norm_data.columns[norm_data.columns.str.startswith("Metadata_")]
    norm_data = norm_data[list(meta_cols) + list(model_columns)]

    # Drop any rows with NaN in the model columns after filtering
    norm_data = norm_data.dropna(subset=model_columns)

    # Get X and y data
    X, y = get_X_y_data(norm_data, label="Metadata_cell_type")

    # Make predictions (probabilities)
    y_pred_final_probs = final_model.predict_proba(X)
    y_pred_shuffled_probs = shuffled_model.predict_proba(X)

    # Add predicted probabilities back to the metadata columns
    predictions_df = norm_data[list(meta_cols)].copy()  # Copy metadata columns

    # Create DataFrames for final and shuffled model predictions
    final_predictions_df = predictions_df.copy()
    final_predictions_df["predicted_probas"] = y_pred_final_probs[
        :, 1
    ]  # Probability of the positive class (healthy)
    final_predictions_df["model_type"] = "final"

    shuffled_predictions_df = predictions_df.copy()
    shuffled_predictions_df["predicted_probas"] = y_pred_shuffled_probs[
        :, 1
    ]  # Probability of the positive class (healthy)
    shuffled_predictions_df["model_type"] = "shuffled"

    # Concatenate the two DataFrames
    predictions_df = pd.concat(
        [final_predictions_df, shuffled_predictions_df], ignore_index=True
    )

    # Add the model_name column
    predictions_df["model_name"] = plate

    # Append the DataFrame to the list
    all_predictions.append(predictions_df)

# Concatenate all predictions into a single DataFrame
all_predictions_df = pd.concat(all_predictions, ignore_index=True)

# Save all predictions to parquet file
all_predictions_df.to_parquet(
    "./performance_metrics/batch1_compound_predictions.parquet", index=False
)

# Print the resulting DataFrame for verification
print(all_predictions_df.shape)
all_predictions_df.head()
