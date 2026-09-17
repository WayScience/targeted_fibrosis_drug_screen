#!/usr/bin/env python
# coding: utf-8

# # Inference for treatments in validation plate
# 
# > Note: Code is adapted from Weishan Li.

# In[1]:


import json
import pathlib

import numpy as np
import pandas as pd
import polars as pl
import statsmodels.api as sm
import yaml
from joblib import load


# In[3]:


random_state = 0
metadata_prefix = "Metadata_"
label_col = "Metadata_cell_type"

# Path to directory with feature selected profiles for the validation plate
path_to_feature_selected_data = pathlib.Path().home() / "mnt" / "bandicoot" /\
    "CFReT_screening_data" / "validation_profiles"
if not path_to_feature_selected_data.exists() and\
    not path_to_feature_selected_data.is_dir():
    raise FileNotFoundError(
        f"Directory {path_to_feature_selected_data} does not exist or is not a directory."
    )

# Find the single-cell feature selected profile for the validation plate
sc_profile_folder = (path_to_feature_selected_data / "single_cell_profiles").resolve()
if not sc_profile_folder.exists() or not sc_profile_folder.is_dir():
    raise FileNotFoundError(
        f"Directory {sc_profile_folder} does not exist or is not a directory."
    )

plate_files = [
    f
    for f in sc_profile_folder.glob("*_sc_feature_selected.parquet")
    if f.is_file()
]
if len(plate_files) != 1:
    raise ValueError(
        f"Expected exactly one feature selected profile for the validation plate, found {len(plate_files)}."
    )
plate_file = plate_files[0]

## Other needed paths
fitted_model_dir = pathlib.Path(".") / "models"
if not fitted_model_dir.exists():
    raise FileNotFoundError(f"Fitted model directory not found: {fitted_model_dir}")

datasplit_dir = pathlib.Path(".") / "datasplits"
if not datasplit_dir.exists():
    raise FileNotFoundError(f"Datasplit directory not found: {datasplit_dir}")

encoding_path = datasplit_dir / "cell_type_encoding.json"
if not encoding_path.exists():
    raise FileNotFoundError(f"Cell type encoding file not found: {encoding_path}")

platemap_csv_dir = pathlib.Path.cwd().parent / "metadata" / "platemaps"
if not platemap_csv_dir.exists():
    raise FileNotFoundError(
        f"Directory {platemap_csv_dir} does not exist."
    )

platemap_files = list(platemap_csv_dir.glob("*.csv"))
if len(platemap_files) != 1:
    raise ValueError(
        f"Expected exactly one platemap CSV file for the validation plate, found {len(platemap_files)}."
    )
platemap_file = platemap_files[0]

# output path
output_dir = pathlib.Path(".") / "inference_results"
output_dir.mkdir(exist_ok=True)


# In[4]:


encoding_dict = json.loads(encoding_path.read_text())
print(f"Loaded cell type encoding for {len(encoding_dict)} cell types.")


# In[ ]:


platemap_df = pd.read_csv(platemap_file)

treatment_well_df = platemap_df.loc[
    :,
    ["well_position", "cell_type", "treatment"]
]
treatment_well_df = treatment_well_df[
    (treatment_well_df["treatment"] != "DMSO") &
    (treatment_well_df["cell_type"] == "failing")
].rename(columns={"well_position": "Metadata_Well"})

# Validation plate has a single plate barcode, derived from the feature selected profile filename
treatment_well_df["Metadata_Plate"] = plate_file.stem.removesuffix("_sc_feature_selected")

treatment_well_df.to_csv(output_dir / "treatment_well.csv")
treatment_well_df.head()


# In[ ]:


score_dfs = []
cell_counts = []
missing_profile_wells = {}

plate_repr = plate_file.stem

# Resolve fitted model checkpoints for the validation plate
saved_model_files = sorted(
    (fitted_model_dir / plate_repr).resolve().glob("*original_statsmodels_logit.joblib")
)
if not saved_model_files:
    raise FileNotFoundError(
        f"No fitted models found for {plate_repr} in {fitted_model_dir}."
    )

# Only inference on non-DMSO treated failing cells
# Note that training is done on cells that are DMSO treated and both failing and healthy
df = (
    pl.scan_parquet(plate_file)
    .filter(
        (pl.col("Metadata_treatment") != "DMSO") & 
        (pl.col("Metadata_cell_type") == "failing")
    )
    .collect(engine="cpu")
    .to_pandas()
)
if df.empty:
    raise ValueError(f"No non-DMSO failing cell rows found for {plate_repr}.")

cell_counts.append(
    df.groupby(['Metadata_treatment', 'Metadata_Well', 'Metadata_Plate']).size().reset_index(name='row_count')
)

# Compute missing treatment wells by cross-referencing with the platemap
plate_name = df['Metadata_Plate'].unique()
if len(plate_name) != 1:
    raise ValueError(f"Multiple plate names found in {plate_repr}: {plate_name}.")
plate_name = plate_name[0]

treatment_wells = treatment_well_df["Metadata_Well"].unique()
missing_treatment_wells = set(treatment_wells) - set(df["Metadata_Well"].unique())
missing_treatments = treatment_well_df[
    treatment_well_df['Metadata_Well'].isin(missing_treatment_wells)
]['treatment'].unique()
# collect missing profile wells info for this plate
missing_profile_wells[plate_name] = list(missing_treatments)

# Iterate over fitted models for this plate and compute scores
scores = []
treatment = df.loc[:, 'Metadata_treatment'].copy()
for saved_model_file in saved_model_files:

    model = load(saved_model_file)

    feats = list(model.params.index)
    feats = [feat for feat in feats if feat != 'const']
    X = df.loc[:, feats].copy()
    X = sm.add_constant(X, has_constant="add")

    score = model.predict(X)
    if encoding_dict['nonfailing'] == 0:
        # model predicts probability of the class encoded as 1 (nonfailing); invert so
        # a larger score means more nonfailing
        score = 1 - score
    elif encoding_dict['nonfailing'] != 1:
        raise ValueError(
            f"Unexpected encoding for nonfailing cell type: {encoding_dict['nonfailing']}"
        )
    scores.append(score)

score_df = pd.concat(
    [treatment] + scores,
    axis=1
)
score_df.columns = [score_df.columns[0]] + [
    f"fold_{i}_score" for i in range(len(scores))
]
score_df['mean_score'] = score_df.iloc[:, 1:].mean(axis=1)
score_df['plate'] = plate_repr
score_df['platemap'] = platemap_file.stem
score_dfs.append(score_df)


# In[ ]:


with (datasplit_dir / "missing_profiles.yaml").open("w") as f:
    yaml.safe_dump(
        missing_profile_wells,
        f,
        default_flow_style=False,
        sort_keys=False,
    )

pd.concat(cell_counts, ignore_index=True).to_csv(
    datasplit_dir / "cell_counts.csv", index=False)

all_score_df = pd.concat(
    score_dfs,
    axis=0,
    ignore_index=True
)
all_score_df.to_parquet(output_dir / 'all_logit_scores.parquet')
all_score_df.head()

