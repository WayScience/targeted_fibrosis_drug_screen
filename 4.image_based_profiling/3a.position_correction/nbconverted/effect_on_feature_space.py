#!/usr/bin/env python
# coding: utf-8

# # Effect of the plate-position correction on the feature space
# 
# `fit_position_correction.ipynb` estimates a plate-position effect (the "tilt") in the single-cell features, and `../3b.apply_position_correction.ipynb` subtracts it from every cell.
# The tilt is an additive row plus column effect of the well position, estimated for each feature.
# The correction subtracts the tilt, scaled by one amplitude for each platemap, from every cell.
# 
# This notebook computes what we need to show the effect of the correction on the features and on the well profiles, and saves it as tables.
# `plot_effect_on_feature_space.ipynb` draws the figures:
# 
# 1. The size of the tilt for each feature family (compartment, measurement type, and channel).
# 2. The wells in the two principal components that depend most on the plate row, before and after the correction, and the displacement of each well.
# 3. The principal component that depends most on the plate row, drawn on the plate grid before and after the correction.
# 
# **Inputs** (in this folder, from `fit_position_correction.ipynb`)
# 
# - `position_correction_fit.npz`: the tilt map and the amplitudes of the platemaps.
# - `well_medians.parquet`: the median profile of every well of every plate before the correction. The fit notebook caches it, and git ignores it.
# 
# The correction shifts all cells of a well by the same vector, so the median profile of a well after the correction is its median before the correction minus that vector.
# We therefore compute the corrected wells from the well medians, without reading single cells.
# 
# **Outputs** (in this folder)
# 
# - `tilt_size_by_feature.csv`: the size of the tilt and the family of every feature.
# - `well_pca_scores.csv`: the principal component scores of every well before and after the correction.
# - `pc_variance_explained_by_row.csv`: the fraction of each component that plate row explains, before and after.

# ## Import libraries

# In[1]:


import pathlib
import sys

import pandas as pd

sys.path.append("../../utils")
import position_correction_effect_utils as pce
import position_correction_utils as pcu


# ## Set paths and parameters

# In[2]:


fit_path = pathlib.Path("position_correction_fit.npz")
well_medians_path = pathlib.Path("well_medians.parquet")
for path in (fit_path, well_medians_path):
    if not path.exists():
        raise FileNotFoundError(f"Input does not exist: {path}. Run fit_position_correction.ipynb first.")

# folder for the tables
results_dir = pathlib.Path(".")

# number of principal components of the uncorrected wells
n_components = 10
# seed of the PCA, so the components are reproducible
random_state = 0


# ## Load the fit and the well medians
# 
# Each row of the well medians is one well of one plate, with its platemap, treatment, and cell type.

# In[3]:


fit = pcu.load_fit(fit_path)
medians = pd.read_parquet(well_medians_path)

missing = set(fit["features"]) - set(medians.columns)
if missing:
    raise KeyError(f"{len(missing)} features of the fit are not in the well medians")

print("Features in the fit:", len(fit["features"]))
print("Wells:", len(medians), "on", medians["plate"].nunique(), "plates")
print(
    "Amplitude of each platemap:",
    {k: round(v, 2) for k, v in fit["amplitudes"].items()},
)


# ## Size of the tilt by feature family
# 
# The size of a feature is the root mean square of its tilt over the 60 positions of the plate (6 rows by 10 columns), with an amplitude of 1.
# The features are z-scored on each plate, so the size is in units of the standard deviation of a single cell.

# In[4]:


sizes = pce.tilt_size_by_feature(fit)
sizes.sort_values("tilt_rms", ascending=False).to_csv(
    results_dir / "tilt_size_by_feature.csv", index=False
)
sizes.sort_values("tilt_rms", ascending=False).head(10)


# ## Wells in principal components, before and after the correction
# 
# We fit the principal components on the uncorrected wells, and project both the uncorrected and the corrected wells into them.
# This puts the two versions in the same space, so we can compare the wells before and after the correction.
# The wells include the compound wells and the DMSO controls.

# In[5]:


features = fit["features"]
corrected = pce.correct_well_medians(medians, fit)

pca = pce.fit_pca(medians[features].to_numpy(dtype=float), n_components, random_state)
before = pca.transform(medians[features].to_numpy(dtype=float))
after = pca.transform(corrected[features].to_numpy(dtype=float))

rows, cols = pcu.well_index(medians["well"])
variance_by_row = pd.DataFrame(
    {
        "component": [f"PC{i + 1}" for i in range(n_components)],
        "variance_explained": pca.explained_variance_ratio_,
        "row_explains_before": [
            pce.variance_explained_by(before[:, i], rows) for i in range(n_components)
        ],
        "row_explains_after": [
            pce.variance_explained_by(after[:, i], rows) for i in range(n_components)
        ],
    }
)
variance_by_row.to_csv(results_dir / "pc_variance_explained_by_row.csv", index=False)

pce.well_score_table(medians, before, after).to_csv(
    results_dir / "well_pca_scores.csv", index=False
)

top_components = variance_by_row.nlargest(2, "row_explains_before")["component"]
print("Components that depend most on plate row:", sorted(top_components))
variance_by_row.round(3)

