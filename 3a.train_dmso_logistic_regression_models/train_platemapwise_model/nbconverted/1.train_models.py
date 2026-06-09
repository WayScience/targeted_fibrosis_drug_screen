#!/usr/bin/env python
# coding: utf-8

# ## Train platemap level logistic regression models to predict failing or healthy cell status
# 
# NOTE: Healthy refers to healthy patient heart and failing refers to patient with dilated cardiomyopathy (heart disease).
# 
# Each training unit consumes a fold split of all pooled plates under the same platemap
# 1. Perform complete quasi-separation check and remove features perfectly predictive of label.
# 2. Perform iterative Recursive Feature Elimination (RFE) to cut down the number of features until the 1 in 10 or 20 rule is satisified per split and plate based on the minority class size.
# 3. Fit the final logistic regression model with post RFE features. 
# 4. Evaluate on testing split and produce ROC and PRC visualizations.
# 
# Per every model fitted as described, a shuffled model is also trained following the exact same procedure and same split data except with labels shuffled.   

# ## Import libraries

# In[1]:


import pathlib
import json
from joblib import Parallel, delayed

import pandas as pd
import numpy as np

from cfret_ml.data_utils import split_and_prep_data
from cfret_ml.orchestrator import process_model_fitting


# ## Pathing and global parameters

# In[2]:


random_state = 0
np.random.seed(random_state)
metadata_prefix = "Metadata_"
label_col = "Metadata_cell_type"

datasplit_dir = pathlib.Path(".") / "datasplits"
if not datasplit_dir.exists():
    raise FileNotFoundError(f"Datasplit directory not found: {datasplit_dir}")

fitted_model_dir = pathlib.Path(".") / "models"
fitted_model_dir.mkdir(exist_ok=True)

eval_plot_dir = pathlib.Path(".") / "eval_plots"
eval_plot_dir.mkdir(exist_ok=True)

encoding_path = datasplit_dir / "cell_type_encoding.json"
if not encoding_path.exists():
    raise FileNotFoundError(f"Cell type encoding file not found: {encoding_path}")


# In[3]:


encoding_dict = json.loads(encoding_path.read_text())
print(f"Loaded cell type encoding for {len(encoding_dict)} cell types.")

platemap_level_splits = [
    p for p in datasplit_dir.iterdir() if p.is_dir()
]
print(f"Found {len(platemap_level_splits)} platemap level splits")
if not platemap_level_splits:
    raise ValueError(f"No platemap level splits found in {datasplit_dir}")


# ## Parallelized model fitting due to RFE being really slow

# In[4]:


tasks = []
split_rows = []

for platemap_dir in platemap_level_splits:

    # Load split files
    split_json_files = list(platemap_dir.glob("*.json"))
    if not split_json_files:
        continue
    dmso_parquet = platemap_dir / "DMSO.parquet"
    if not dmso_parquet.exists():
        continue
    dmso_df = pd.read_parquet(dmso_parquet)
    dmso_df['Metadata_cell_type'] = dmso_df['Metadata_cell_type'].map(encoding_dict)

    platemap_repr = platemap_dir.name
    print(f"Queueing tasks for plate {platemap_repr} with {len(split_json_files)} splits")

    # Make directories for fitted models and eval plots for this platemap
    platemap_fitted_model_dir = fitted_model_dir / platemap_repr
    platemap_fitted_model_dir.mkdir(exist_ok=True)

    platemap_eval_plot_dir = eval_plot_dir / platemap_repr 
    platemap_eval_plot_dir.mkdir(exist_ok=True)

    # Initialize train/test splits with None to ensure correct ordering by fold index
    train_splits = [None] * len(split_json_files)
    test_splits = [None] * len(split_json_files)

    # Process each split JSON file to extract train/test indices and summary statistics
    for split_json in split_json_files:

        with open(split_json, "r") as f:
            split_info = json.load(f)

        train_splits[split_info['fold']] = split_info["train_index"]
        test_splits[split_info['fold']] = split_info["test_index"]

        split_row = {
            "platemap": platemap_repr,
            "fold": split_info["fold"],
            "train_n": len(split_info["train_index"]),
            "test_n": len(split_info["test_index"]),
            **{
                f"{split}_{label}": split_info[split][label]
                for split in ["train", "test"]
                for label in split_info[split]
            }
        }
        split_rows.append(split_row)

    # Iterate through the splits in order of fold index to create model fitting tasks
    for fold_idx, (train_idx, test_idx) in enumerate(zip(train_splits, test_splits)):

        if train_idx is None or test_idx is None:
            continue

        # Helper function to create the train/test profiles, labels and shuffled labels for a given split
        (
            train_profiles,
            test_profiles,
            train_labels,
            test_labels,
            train_labels_shuffled,
        ) = split_and_prep_data(
            dmso_df, 
            train_idx, 
            test_idx, 
            shuffle_random_state=random_state,
            metadata_prefix=metadata_prefix,
            label_col=label_col
        )

        # ensure both train and test sets have at least some representation of both classes 
        n_pos = train_labels.sum()
        n_neg = len(train_labels) - n_pos
        min_class = min(n_pos, n_neg)
        if min_class <= 50 or (n_pos + n_neg <= 100):
            print(f"\tNot enough train samples in platemap {platemap_repr} fold {fold_idx}")
            continue

        n_pos_test = test_labels.sum()
        n_neg_test = len(test_labels) - n_pos_test
        if n_pos_test == 0 or n_neg_test == 0:
            print(f"\tTest set missing a class in platemap {platemap_repr} fold {fold_idx}")
            continue

        # Create two tasks for this split - one with original labels and one with shuffled labels
        for shuffle_status, labels in zip(
            ["original", "shuffled"],
            [train_labels, train_labels_shuffled]
        ):
            tasks.append({
                "train_profiles": train_profiles,
                "test_profiles": test_profiles,
                "labels": labels,
                "test_labels": test_labels,
                "shuffle_status": shuffle_status,
                "plate_repr": platemap_repr,
                "fold": fold_idx,
                "plate_fitted_model_dir": platemap_fitted_model_dir,
                "plate_eval_plot_dir": platemap_eval_plot_dir,
                "random_state": random_state
            })

# Run all the model fitting tasks in parallel and collect results into a dataframe, then merge with the split summary statistics and save to CSV
print(f"Executing {len(tasks)} model fitting tasks in parallel (n_jobs=8)...")
results = Parallel(n_jobs=8)(delayed(process_model_fitting)(**kwargs) for kwargs in tasks)

results_df = pd.DataFrame(results)
split_df = pd.DataFrame(split_rows)
enriched_df = pd.merge(
    results_df.rename(columns={"plate": "platemap"}),
    split_df,
    on=["platemap", "fold"],
    how="left",
)
enriched_df.to_csv(eval_plot_dir / "model_fit_summary.csv", index=False)

