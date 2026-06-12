#!/usr/bin/env python
# coding: utf-8

# # Generate well-level bulk profiles profiles
# 
# NOTE: We are normalizing the bulk-profile plates to the negative controls as the "standard".

# ## Import libraries

# In[1]:


import os
import pathlib
import pprint

import pandas as pd

from pycytominer import aggregate, annotate, normalize, feature_select


# ## Set paths and variables

# In[ ]:


# get the batch to process from environment variable
batch_to_process = os.environ.get("BATCH", "batch_1")
if batch_to_process is None:
    raise ValueError(
        "Please set the BATCH environment variable before running this script."
    )

# base directory where batches are located
base_dir = pathlib.Path("./data/").resolve(strict=True)

# Decide what to process
if batch_to_process:
    print(f"Processing {batch_to_process}")
    batch_dirs = [base_dir / batch_to_process]
else:
    print("No specific batch set, processing all available batches")
    batch_dirs = [p for p in base_dir.glob("batch_*") if p.is_dir()]

# path for platemap directory
platemap_dir = pathlib.Path("../metadata/updated_platemaps/")

# Load the barcode_platemap file
barcode_platemap_df = pd.read_csv(
    (platemap_dir / "updated_barcode_platemap.csv").resolve()
)

# operations to perform for feature selection
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "blocklist",
    "drop_na_columns",
]


# ## Set dictionary with plates to process

# In[3]:


plate_info_dictionary = {}

# Loop over batches and layouts
for batch_dir in batch_dirs:
    layouts = [p for p in batch_dir.iterdir() if p.is_dir()]  # all layouts
    for layout_dir in layouts:
        qc_labeled_dir = layout_dir / "qc_labeled_profiles"
        output_dir = layout_dir / "bulk_profiles"
        output_spherize_dir = layout_dir / "spherized_bulk_profiles"

        # Create directories once per layout
        for d in [output_dir, output_spherize_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Extract plate names from parquet files
        parquet_files = list(qc_labeled_dir.glob("*.parquet"))
        plate_names = [
            "_".join(f.stem.split("_")[:2]) if len(f.stem.split("_")) >= 2 else f.stem
            for f in parquet_files
        ]

        for name in plate_names:
            # Find the corresponding parquet file
            matching_files = [f for f in parquet_files if name in f.stem]
            if not matching_files:
                continue
            profile_path = matching_files[0].resolve(strict=True)

            # Find corresponding platemap CSV
            platemap_row = barcode_platemap_df.loc[
                barcode_platemap_df["plate_barcode"] == name
            ]
            if platemap_row.empty:
                raise ValueError(f"No platemap found for plate {name}")
            platemap_path = (
                platemap_dir / f"{platemap_row['platemap_file'].values[0]}.csv"
            ).resolve(strict=True)

            # Add to dictionary
            plate_info_dictionary[name] = {
                "profile_path": profile_path,
                "platemap_path": platemap_path,
                "output_dir": output_dir,
                "spherize_output_dir": output_spherize_dir,
            }

# View dictionary
print("Number of plates to process:", len(plate_info_dictionary))
pprint.pprint(plate_info_dictionary, indent=4)


# ## Process data with pycytominer
# 

# In[4]:


for plate, info in plate_info_dictionary.items():
    output_dir = info["output_dir"]
    print("Performing preprocessing on", plate, output_dir)

    # Use the output_dir from the dictionary for this specific plate
    output_dir = info["output_dir"]

    # generating all output file paths using the output_dir from the dictionary
    output_annotated_file = str(output_dir / f"{plate}_bulk_annotated.parquet")
    output_normalized_file = str(output_dir / f"{plate}_bulk_normalized.parquet")
    output_feature_select_file = str(
        output_dir / f"{plate}_bulk_feature_selected.parquet"
    )

    # loading profiles
    profile_df = pd.read_parquet(info["profile_path"])
    platemap_df = pd.read_csv(info["platemap_path"])

    # Drop all rows in the profiles that failed any Metadata_cqc columns
    cqc_columns = [col for col in profile_df.columns if col.startswith("Metadata_cqc")]
    if cqc_columns:
        profile_df = profile_df[~profile_df[cqc_columns].any(axis=1)]

    # Step 1: Aggregate single-cell data to the well-level using the median
    print("Performing aggregation for", plate, "...")
    aggregated_df = aggregate(
        population_df=profile_df,
        operation="median",
        strata=["Image_Metadata_Plate", "Image_Metadata_Well"],
    )

    # Step 2: Annotation
    print("Performing annotation for", plate, "...")
    annotate(
        profiles=aggregated_df,
        platemap=platemap_df,
        join_on=["Metadata_well_position", "Image_Metadata_Well"],
        output_type="parquet",
        output_file=output_annotated_file,
    )

    # Load the annotated parquet file to fix metadata columns names
    annotated_df = pd.read_parquet(output_annotated_file)

    # Rename columns
    annotated_df.rename(columns={"Image_Metadata_Site": "Metadata_Site"}, inplace=True)

    # Save annotated profiles back to parquet
    annotated_df.to_parquet(output_annotated_file, index=False)

    # Step 3: Normalization (mad robustize)
    # Per plate normlization using the negative controls as the reference population
    print("Performing normalization for", plate, "...")
    neg_control_query = "Metadata_treatment == 'DMSO' and Metadata_cell_type == 'failing'"
    normalize(
        profiles=annotated_df,
        method="mad_robustize",
        samples=neg_control_query,
        output_type="parquet",
        output_file=output_normalized_file,
    )

    # Step 4: Per plate feature selection
    print("Performing feature selection for", plate, "...")
    feature_select(
        profiles=output_normalized_file,
        operation=feature_select_ops,
        na_cutoff=0,
        blocklist_file="./blocklist_features.txt",
        corr_threshold=0.95,
        freq_cut=0.05,
        output_type="parquet",
        output_file=output_feature_select_file,
    )

    print(
        f"Aggregation, annotation, normalization, and feature selection complete for {plate}"
    )


# ## spherizing profiles
# This pools all normalized replicate plates of the same platmap, performs features, and then spherizes/whitens the profiles. Spherization uses the negative-control wells as the reference population to decorrelate the features and put profiles on a shared control-based covariance scale.
# 
# Parameters:
# - `neg_control_query`: pandas query string that selects the control wells used to fit both standardization and spherization. Here, the reference population is failing-cell DMSO wells.
# 
# - `negcon_fs`: if True, run an extra feature-selection pass using only the negative-control wells. This removes features that are static or nearly static inside the exact control population used for spherization. Keep False unless spherization fails because of zero-variance control features.

# In[5]:


# parameters
neg_control_query = "Metadata_treatment == 'DMSO' and Metadata_cell_type == 'failing'"
negcon_fs = True


# In[6]:


# Group annotated replicate profiles by plate map for pooled normalization and 
# spherization outputs. 
# {"platemap_1": {"replicate_paths": [path1, path2], "spherized_output_dir": path}, ...}
normalized_replicate_plates = {}

for plate_barcode, info in plate_info_dictionary.items():
    # Example key: platemap_1, platemap_2, etc.
    plate_key = info["output_dir"].parent.name

    # Path to the well-level normalized profile generated in the previous step.
    annotated_plate_path = (
        info["output_dir"] / f"{plate_barcode}_bulk_normalized.parquet"
    )

    # Start a new group the first time we see this plate map.
    if plate_key not in normalized_replicate_plates:
        normalized_replicate_plates[plate_key] = {
            "replicate_paths": [],
            "spherized_output_dir": info["spherize_output_dir"],
        }

    # Add this physical replicate plate to its plate-map group.
    normalized_replicate_plates[plate_key]["replicate_paths"].append(
        annotated_plate_path
    )

# sort based on plate key
normalized_replicate_plates = dict(sorted(normalized_replicate_plates.items()))


# In[7]:


# Spherizing step:
# Process each plate map after pooling its replicate plates.
for plate_key, plate_info in normalized_replicate_plates.items():
    replicate_paths = sorted(plate_info["replicate_paths"])
    spherized_output_dir = plate_info["spherized_output_dir"]

    # All replicate paths in this group share the same bulk output directory.
    output_dir = replicate_paths[0].parent
    spherized_output_dir.mkdir(parents=True, exist_ok=True)

    # Platemap-level outputs created from pooled replicate plates.
    output_normalized_file = (
        output_dir / f"{plate_key}_replicate_bulk_normalized.parquet"
    )
    output_feature_select_file = (
        output_dir / f"{plate_key}_replicate_bulk_feature_selected.parquet"
    )
    output_spherized_file = (
        spherized_output_dir / f"{plate_key}_replicate_bulk_spherized.parquet"
    )

    # step 1: concat mad-normalized replicate plates before normalization to avoid 
    # centering
    # each physical plate's small DMSO control set independently.
    print(f"Concatenating replicate plates for {plate_key}...")
    concat_replicate_df = pd.concat(
        [pd.read_parquet(path) for path in replicate_paths],
        ignore_index=True,
    ).reset_index(drop=True)

    # step 2a: Apply feature selection on the pooled replicate plates to get a common 
    # set of features for spherization.
    print(f"Feature selecting {plate_key}...")
    feature_select(
        profiles=concat_replicate_df,
        operation=feature_select_ops,
        na_cutoff=0,
        blocklist_file="./blocklist_features.txt",
        corr_threshold=0.95,
        freq_cut=0.05,
        output_type="parquet",
        output_file=output_feature_select_file,
    )

    if negcon_fs:
        # step 2b: Remove features with too little variation inside the exact control
        # population used to fit spherization.
        print(f"Feature selecting {plate_key} with variance threshold...")
        # only use this if spherization fails due to features with zero variance in the 
        # control population.
        feature_select(
            profiles=output_feature_select_file,
            operation="variance_threshold",
            freq_cut=0.05,
            unique_cut=0.01,
            samples=neg_control_query,
            output_file=output_feature_select_file,
            output_type="parquet",
        )

    # step 3: Spherize/whiten all profiles using the pooled negative controls as the
    # reference population.
    print(f"Spherizing {plate_key} using pooled negative controls...")
    normalize(
        profiles=output_feature_select_file,
        method="spherize",
        samples=neg_control_query,
        spherize_center=True,
        spherize_method="ZCA-cor",
        spherize_epsilon=1e-6,
        output_file=output_spherized_file,
        output_type="parquet",
    )

    print(f"Saved feature-selected profiles to {output_feature_select_file}")
    print(f"Saved spherized profiles to {output_spherized_file}")


# In[8]:


# Check output file
test_df = pd.read_parquet(output_feature_select_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)


# In[9]:


# Check output file
test_df = pd.read_parquet(output_annotated_file)

print(test_df.shape)
print("Plate:", test_df.Metadata_Plate.unique())
print(
    "Metadata columns:", [col for col in test_df.columns if col.startswith("Metadata_")]
)
test_df.head(2)

