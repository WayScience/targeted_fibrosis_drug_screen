#!/usr/bin/env python
# coding: utf-8

# ## Train plate specific logistic regression models to predict failing or healthy cell status
# 
# Each training unit consumes a fold split of a plate and does the following:
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

# In[ ]:


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


# ## Parallelized model fitting due to RFE being really slow

# In[3]:


encoding_path = datasplit_dir / "cell_type_encoding.json"
if not encoding_path.exists():
    raise FileNotFoundError(f"Cell type encoding file not found: {encoding_path}")
encoding_dict = json.loads(encoding_path.read_text())
print(f"Loaded cell type encoding for {len(encoding_dict)} cell types.")

plate_level_splits = [
    p for p in datasplit_dir.iterdir() if p.is_dir()
]
print(f"Found {len(plate_level_splits)} plate level splits")
if not plate_level_splits:
    raise ValueError(f"No plate level splits found in {datasplit_dir}")

tasks = []
split_rows = []
for plate_dir in plate_level_splits:
    split_json_files = list(plate_dir.glob("*.json"))
    if not split_json_files:
        continue
    dmso_parquet = plate_dir / "DMSO.parquet"
    if not dmso_parquet.exists():
        continue
    dmso_df = pd.read_parquet(dmso_parquet)
    dmso_df['Metadata_cell_type'] = dmso_df['Metadata_cell_type'].map(encoding_dict)

    plate_repr = plate_dir.name
    print(f"Queueing tasks for plate {plate_repr} with {len(split_json_files)} splits")

    plate_fitted_model_dir = fitted_model_dir / plate_repr
    plate_fitted_model_dir.mkdir(exist_ok=True)

    plate_eval_plot_dir = eval_plot_dir / plate_repr
    plate_eval_plot_dir.mkdir(exist_ok=True)

    train_splits = [None] * len(split_json_files)
    test_splits = [None] * len(split_json_files)
    for split_json in split_json_files:
        with open(split_json, "r") as f:
            split_info = json.load(f)

        train_splits[split_info['fold']] = split_info["train_index"]
        test_splits[split_info['fold']] = split_info["test_index"]

        split_row = {
            "plate": plate_repr,
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

    for fold_idx, (train_idx, test_idx) in enumerate(zip(train_splits, test_splits)):

        if train_idx is None or test_idx is None:
            continue

        # basic data prep
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

        n_pos = train_labels.sum()
        n_neg = len(train_labels) - n_pos
        min_class = min(n_pos, n_neg)
        if min_class <= 50 or (n_pos + n_neg <= 100):
            print(f"\tNot enough train samples in plate {plate_repr} fold {fold_idx}")
            continue

        n_pos_test = test_labels.sum()
        n_neg_test = len(test_labels) - n_pos_test
        if n_pos_test == 0 or n_neg_test == 0:
            print(f"\tTest set missing a class in plate {plate_repr} fold {fold_idx}")
            continue

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
                "plate_repr": plate_repr,
                "fold": fold_idx,
                "plate_fitted_model_dir": plate_fitted_model_dir,
                "plate_eval_plot_dir": plate_eval_plot_dir,
                "random_state": random_state
            })

print(f"Executing {len(tasks)} model fitting tasks in parallel (n_jobs=8)...")
results = Parallel(n_jobs=8)(delayed(process_model_fitting)(**kwargs) for kwargs in tasks)

results_df = pd.DataFrame(results)
split_df = pd.DataFrame(split_rows)
enriched_df = pd.merge(
    results_df,
    split_df,
    on=["plate", "fold"],
    how="left",
)
enriched_df.to_csv(eval_plot_dir / "model_fit_summary.csv", index=False)

