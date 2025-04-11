#!/usr/bin/env python
# coding: utf-8

# ## Train machine learning models to predict failing or healthy cell status
# 
# Each model will be trained on an individual plate or all plates from a batch combined.

# ## Import libraries

# In[1]:


import pathlib
import pprint
import sys
import warnings

import numpy as np
import pandas as pd
from joblib import dump
from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import parallel_backend

sys.path.append("../utils")
from training_utils import downsample_data, get_X_y_data


# ## Set paths and variables

# In[2]:


# set numpy seed to make sure any random operations performs are reproducible
np.random.seed(0)

# path to training/testing datasets
training_data_path = pathlib.Path("./data")

# Find all training datasets
training_files = list(training_data_path.glob("*_train.parquet"))

# Metadata column used for prediction class
label = "Metadata_cell_type"

# Directory for models to be outputted
model_dir = pathlib.Path("./models")
model_dir.mkdir(exist_ok=True, parents=True)

# Directory for label encoder
encoder_dir = pathlib.Path("./encoder_results")
encoder_dir.mkdir(exist_ok=True, parents=True)

# Directory for training indices
training_indices_dir = pathlib.Path("./training_indices")
training_indices_dir.mkdir(exist_ok=True, parents=True)


# ## Load in training data

# In[3]:


# Load in all training files as dataframes and store them under 'orig_train_df'
training_data_dfs_dict = {
    "_".join(parts[:2]) if parts[0] == "combined" else parts[0]: {
        "orig_train_df": pd.read_parquet(file)
    }
    for file in training_files
    if (parts := pathlib.Path(file).stem.split("_"))
}

# Pretty print the dictionary
pprint.pprint(training_data_dfs_dict, indent=4)


# ## Perform downsampling on training data and output as data frame

# In[4]:


# for loop to process each dataframe in the dict
for plate, info in training_data_dfs_dict.items():
    # load in training plate 4 data as downsampled to lowest class
    downsample_df = downsample_data(data=info["orig_train_df"], label=label)

    # Store the downsampled dataframe under 'downsample_train_df'
    training_data_dfs_dict[plate]["downsample_train_df"] = downsample_df

    # Export sample indices used in training the model to a new one-column CSV file
    output_file = f"{training_indices_dir}/{plate}_training_data_indices.csv"
    pd.DataFrame(downsample_df.index, columns=["Index"]).to_csv(
        output_file, index=False
    )

    print(f"CSV file created at {output_file} with {len(downsample_df.index)} entries.")

    print(downsample_df.shape)
    print(downsample_df["Metadata_cell_type"].value_counts())


# ## Get X and y data and label encoder for final and shuffled models for all plates

# In[5]:


# Encode classes
le = LabelEncoder()

for plate, info in training_data_dfs_dict.items():
    # Get downsampled dataframe
    downsample_df = info["downsample_train_df"]

    # Get not shuffled training data from downsampled df (e.g., "final")
    X_train, y_train = get_X_y_data(df=downsample_df, label=label, shuffle=False)

    # Print out the number of features the model will train on per plate
    print(f"Number of features for plate {plate}: {X_train.shape[1]}")

    # Fit the LabelEncoder on the non-shuffled labels
    le.fit(y_train)

    # Encode the labels for both non-shuffled data
    y_train_encoded = le.transform(y_train)

    # Get shuffled training data from downsampled df(e.g., "shuffled_baseline")
    X_shuffled_train, y_shuffled_train = get_X_y_data(
        df=downsample_df, label=label, shuffle=True
    )

    # Encode the labels for the shuffled labels
    y_shuffled_train_encoded = le.transform(y_shuffled_train)

    # Store the X and y data under respective keys
    training_data_dfs_dict[plate]["X_train"] = X_train
    training_data_dfs_dict[plate]["y_train"] = y_train
    training_data_dfs_dict[plate]["X_shuffled_train"] = X_shuffled_train
    training_data_dfs_dict[plate]["y_shuffled_train"] = y_shuffled_train_encoded

    # Save label encoder
    dump(le, f"{encoder_dir}/label_encoder_{plate}.joblib")

# Print the class mapping to see the encoding
class_mapping = dict(zip(le.classes_, le.transform(le.classes_)))
print("Class Mapping:")
print(class_mapping)


# ## Train the models
# 
# These hyperparameters are set based on the model training from the [`cellpainting_predicts_cardiac_fibrosis` repository](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis). 
# The following model training code is derived from the [model training notebook](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis/blob/main/5.machine_learning/0.train_logistic_regression/1.train_models.ipynb).
# We will be using RandomizedSearchCV to hyperparameterize the model since that is how the original model was trained and we want to remain consistent.

# ### Set up the model and hyper parameter method

# In[6]:


# Set folds for k-fold cross validation (default is 5)
straified_k_folds = StratifiedKFold(n_splits=10, shuffle=False)

# Set Logistic Regression model parameters (use default for max_iter)
logreg_params = {
    "penalty": "elasticnet",
    "solver": "saga",
    "max_iter": 1000,
    "n_jobs": -1,
    "random_state": 0,
    "class_weight": "balanced",
}

# Define the hyperparameter search space for RandomizedSearchCV
param_dist = {
    "C": np.logspace(-3, 3, 7),
    "l1_ratio": np.linspace(0, 1, 11),
}

# Set the random search hyperparameterization method parameters (used default for "cv" and "n_iter" parameter)
random_search_params = {
    "param_distributions": param_dist,
    "scoring": "f1_weighted",
    "random_state": 0,
    "n_jobs": -1,
    "cv": straified_k_folds,
}


# ### Train final and shuffled models per plate and combined batch

# In[7]:


# Initialize Logistic Regression and RandomizedSearchCV
logreg = LogisticRegression(**logreg_params)
random_search = RandomizedSearchCV(logreg, **random_search_params)

# Loop through the training data dictionary for both non-shuffled and shuffled data
for plate, info in training_data_dfs_dict.items():
    # Get the non-shuffled and shuffled data for the current feature type
    X_train = info["X_train"]
    y_train = info["y_train"]
    X_shuffled_train = info["X_shuffled_train"]
    y_shuffled_train = info["y_shuffled_train"]

    # Prevent the convergence warning in sklearn, it does not impact the result
    with parallel_backend("multiprocessing"):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=ConvergenceWarning, module="sklearn"
            )
            ########################################################
            # Train the model for non-shuffled (final) training data
            ########################################################
            print(f"Training model for {plate} features (final)...")
            final_random_search = clone(random_search)
            final_random_search.fit(X_train, y_train)
            print(
                f"Optimal parameters for {plate} features (final):",
                final_random_search.best_params_,
            )

            # Save the model for non-shuffled/final data using joblib
            final_model_filename = model_dir / f"{plate}_final_downsample.joblib"
            dump(final_random_search.best_estimator_, final_model_filename)
            print(f"Model saved as: {final_model_filename}")

            ########################################################
            # Train the model for shuffled training data
            ########################################################
            print(f"Training model for {plate} features (shuffled)...")
            shuffled_random_search = clone(random_search)
            shuffled_random_search.fit(X_shuffled_train, y_shuffled_train)
            print(
                f"Optimal parameters for {plate} features (shuffled):",
                shuffled_random_search.best_params_,
            )

            # Save the final model for shuffled data using joblib
            shuffled_final_model_filename = (
                model_dir / f"{plate}_shuffled_downsample.joblib"
            )
            dump(shuffled_random_search.best_estimator_, shuffled_final_model_filename)
            print(f"Model saved as: {shuffled_final_model_filename}")

