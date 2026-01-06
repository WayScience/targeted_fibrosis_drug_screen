#!/usr/bin/env python
# coding: utf-8

# ## Generate QC report of the whole image quality control flags
# 
# We rerun the notebook per platemap layout.

# ## Import libraries

# In[1]:


import pandas as pd
import pathlib
import re

import seaborn as sns
from upsetplot import from_indicators, plot
import matplotlib.pyplot as plt

import warnings

# Ignore upset plot warning regarding  behavior change that will occur in pandas 3.0
warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")


# # Set paths and load in data

# In[2]:


# Platemap name for processing (e.g., platemap_#)
platemap_name = "platemap_4"

# Set constants for prefixes we want to keep for the dataframes
prefixes = ("Metadata_", "ImageQuality_PercentMaximal", "ImageQuality_PowerLogLogSlope")

# Output directory for plots
output_directory = pathlib.Path(f"./qc_plots/{platemap_name}")
output_directory.mkdir(exist_ok=True)

# path to the qc results
qc_results = pathlib.Path(f"./qc_results/{platemap_name}")
# Find all Image.csv files in the directory
csv_files = list(qc_results.rglob("Image.csv"))
print(f"Found {len(csv_files)} CSV files for {platemap_name}.")


# In[3]:


dataframes = []
problem_files = []

for csv_file in csv_files:
    df = pd.read_csv(csv_file, usecols=lambda col: col.startswith(prefixes))

    unique_plates = df["Metadata_Plate"].unique()
    correct_plate_name = csv_file.parent.name

    # Check if multiple plate names exist
    if len(unique_plates) > 1:
        problem_files.append((csv_file, unique_plates))
        print(
            f"⚠️ {csv_file} has multiple plate names {unique_plates}, fixing to {correct_plate_name}"
        )

    # Overwrite the plate column to the folder name
    df["Metadata_Plate"] = correct_plate_name

    # Verify the fix
    assert (
        df["Metadata_Plate"].nunique() == 1
        and df["Metadata_Plate"].iloc[0] == correct_plate_name
    ), f"Failed to fix {csv_file}"

    dataframes.append(df)

if problem_files:
    print("\nSummary: CSV files with multiple plate names were fixed.")
else:
    print("✅ All CSV files had consistent plate names.")


# In[4]:


# Concatenate the CSV dataframes
combined_df = pd.concat(dataframes, ignore_index=True)

print(combined_df.shape)
combined_df.head()


# ## Update plate names to be easier to read

# In[5]:


# Extract numeric suffixes directly into a list
plate_names = combined_df["Metadata_Plate"].unique()
numeric_suffixes = [int(re.search(r"(\d+)$", name).group(1)) for name in plate_names]

# Zip names with their numeric part, sort by numeric value
sorted_plates = [name for _, name in sorted(zip(numeric_suffixes, plate_names))]

# Create the mapping to Plate_1, Plate_2, ...
plate_mapping = {name: f"Plate_{i+1}" for i, name in enumerate(sorted_plates)}

# Apply mapping
combined_df["Metadata_Plate_Alias"] = combined_df["Metadata_Plate"].map(plate_mapping)

print("Unique plate names:")
print(combined_df["Metadata_Plate"].unique())

print("\nUnique plate aliases:")
print(combined_df["Metadata_Plate_Alias"].unique())

print(combined_df.shape)
combined_df.head()


# ## Add flags per channel

# In[6]:


# Start by assuming no channel is flagged
combined_df["Flagged_Saturation_Channel"] = "none"

# Set the columns for percent maximal (saturation)
percent_max_cols = [
    col for col in combined_df.columns if col.startswith("ImageQuality_PercentMaximal")
]

# Boolean mask where the Saturation flag is on
saturation_flagged_mask = combined_df["Metadata_Saturation_Flag"].astype(bool)

# Iterate through each PercentMaximal column and create a new boolean column per channel
for col in percent_max_cols:
    channel = col.replace("ImageQuality_PercentMaximal_", "")
    combined_df[f"{channel}_Saturated"] = saturation_flagged_mask & (
        combined_df[col] > 0.10
    )

print(combined_df.shape)
combined_df.head()


# In[7]:


# dictionary with blur thresholds per channel
blur_thresholds = {
    "OrigActin": -1.8891791699802942,
    "OrigDNA": -2.2456075474546515,
    "OrigER": -2.2825812279725524,
    "OrigMito": -2.012531942517173,
    "OrigPM": -2.4309820530015642,
}

# Boolean mask for Blur flag
blur_flagged_mask = combined_df["Metadata_Blur_Flag"].astype(bool)

# Iterate through PowerLogLogSlope columns and apply thresholds
for col in combined_df.columns:
    if col.startswith("ImageQuality_PowerLogLogSlope_"):
        channel = col.replace("ImageQuality_PowerLogLogSlope_", "")
        if channel in blur_thresholds:
            threshold = blur_thresholds[channel]
            combined_df[f"{channel}_Blur"] = blur_flagged_mask & (
                combined_df[col] < threshold
            )

print(combined_df.shape)
combined_df.head()


# In[8]:


# Ensure Failed_Any column exists
combined_df["Failed_Any"] = combined_df[
    ["Metadata_Blur_Flag", "Metadata_Saturation_Flag"]
].any(axis=1)

# Calculate total FOVs and number of failed FOVs per plate
plate_fov_counts = (
    combined_df.groupby("Metadata_Plate_Alias")
    .agg(Total_FOVs=("Failed_Any", "size"), Failed_FOVs=("Failed_Any", "sum"))
    .reset_index()
)

print(plate_fov_counts)


# In[9]:


# Count image-sets that failed due to blur only (blur True, saturation False)
blur_only_mask = combined_df["Metadata_Blur_Flag"].astype(bool) & ~combined_df[
    "Metadata_Saturation_Flag"
].astype(bool)
blur_only = combined_df.loc[blur_only_mask].copy()

# Total number
total_blur_only = len(blur_only)
print(f"Total image-sets failing for blur only: {total_blur_only}")

# Breakdown per plate (if multiple plates exist)
blur_only_per_plate = (
    blur_only.groupby("Metadata_Plate_Alias")
    .size()
    .reset_index(name="Blur_Only_Failed")
)
print("\nPer-plate counts:")
print(blur_only_per_plate.to_string(index=False))

# If plate_fov_counts is available, compute percent of FOVs per plate that failed for blur only
if "plate_fov_counts" in globals():
    pct = pd.merge(
        blur_only_per_plate,
        plate_fov_counts[["Metadata_Plate_Alias", "Total_FOVs"]],
        on="Metadata_Plate_Alias",
        how="left",
    )
    pct["Percent_Blur_Only"] = (pct["Blur_Only_Failed"] / pct["Total_FOVs"]) * 100
    print("\nPer-plate percent of FOVs failing for blur only:")
    print(pct.to_string(index=False))


# ## Plot the percentage of failed FOVs across plates regardless of condition

# In[10]:


# Calculate percentage of rows that failed any QC check, grouped by plate
combined_df["Failed_Any"] = combined_df[
    ["Metadata_Blur_Flag", "Metadata_Saturation_Flag"]
].any(axis=1)
failed_percent_by_plate = (
    combined_df.groupby("Metadata_Plate_Alias")["Failed_Any"]
    .mean()
    .reset_index(name="Percent_Failed")
)
failed_percent_by_plate["Percent_Failed"] *= 100

# Create bar plot with a single color
plt.figure(figsize=(6, 5))
ax = sns.barplot(
    data=failed_percent_by_plate,
    x="Metadata_Plate_Alias",
    y="Percent_Failed",
    color="steelblue",
)

# Add percentage labels on top of bars
for p in ax.patches:
    height = p.get_height()
    ax.text(
        x=p.get_x() + p.get_width() / 2,
        y=height + 1,  # slightly above the bar
        s=f"{height:.1f}%",
        ha="center",
        va="bottom",
        fontsize=10,
    )

plt.ylabel("Percentage of FOVs failing QC (%)")
plt.title("Proportion of FOVs failing QC per plate")
plt.ylim(0, 100)
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(output_directory / "qc_failure_rate_per_plate.png", dpi=500)
plt.show()


# ## Plot percentage failed FOV based on blur or saturation regardless of channel

# In[11]:


# Add a column for "failed both"
combined_df["Failed_Both"] = combined_df["Metadata_Saturation_Flag"].astype(
    bool
) & combined_df["Metadata_Blur_Flag"].astype(bool)

# Melt all three flag columns
flag_counts = combined_df.melt(
    id_vars="Metadata_Plate_Alias",
    value_vars=["Metadata_Saturation_Flag", "Metadata_Blur_Flag", "Failed_Both"],
    var_name="QC_Flag",
    value_name="Flagged",
)

# Convert to boolean if needed
flag_counts["Flagged"] = flag_counts["Flagged"].astype(bool)

# Clean up the legend labels
flag_counts["QC_Flag"] = flag_counts["QC_Flag"].map(
    {
        "Metadata_Saturation_Flag": "Failed Saturation",
        "Metadata_Blur_Flag": "Failed Blur",
        "Failed_Both": "Failed Both",
    }
)

# Count total per plate per flag type
total_counts = (
    flag_counts.groupby(["Metadata_Plate_Alias", "QC_Flag"])
    .size()
    .reset_index(name="Total")
)

# Count number of failed flags
fail_counts = (
    flag_counts[flag_counts["Flagged"]]
    .groupby(["Metadata_Plate_Alias", "QC_Flag"])
    .size()
    .reset_index(name="Failed")
)

# Merge and compute percentage
qc_summary = pd.merge(
    total_counts, fail_counts, on=["Metadata_Plate_Alias", "QC_Flag"], how="left"
).fillna(0)
qc_summary["Percent_Failed"] = (qc_summary["Failed"] / qc_summary["Total"]) * 100

# Plot
sns.set_theme(style="whitegrid")
plt.figure(figsize=(8, 6))
sns.barplot(
    data=qc_summary, x="Metadata_Plate_Alias", y="Percent_Failed", hue="QC_Flag"
)
plt.xticks(rotation=45, ha="right")
plt.ylabel("Percent failed FOV (%)")
plt.title("Proportion of FOVs failing QC based on condition per plate")
plt.legend(
    title="QC flag",
    bbox_to_anchor=(0.5, -0.25),
    loc="upper center",
    ncol=3,  # number of columns for a horizontal layout
)
plt.tight_layout()
plt.subplots_adjust(bottom=0.28)  # reduce the whitespace belo
plt.savefig(output_directory / "qc_failure_rate_by_flag_type_per_plate.png", dpi=500)
plt.show()


# ## Create upset plot for all plates in combination with the breakdown of channels failing QC

# In[12]:


# Define the relevant QC columns
qc_columns = [
    "OrigActin_Saturated",
    "OrigDNA_Saturated",
    "OrigER_Saturated",
    "OrigMito_Saturated",
    "OrigPM_Saturated",
    "OrigActin_Blur",
    "OrigDNA_Blur",
    "OrigER_Blur",
    "OrigMito_Blur",
    "OrigPM_Blur",
]

# Make sure columns are boolean
qc_data = combined_df[qc_columns].astype(bool)

# Create the upset input
upset_data = from_indicators(qc_columns, qc_data)

# Create the figure with desired size
fig = plt.figure(figsize=(14, 10))  # Adjust the size as needed

# Plot the UpSet plot
plot(
    upset_data,
    fig=fig,
    element_size=None,
    show_counts=True,
    show_percentages=True,
    min_subset_size=75,
)
plt.suptitle("UpSet plot of channel-wise QC failures\nacross plates in platemap layout")
plt.savefig(
    output_directory / "upset_plot_channel_qc_failures.png",
    bbox_inches="tight",
    dpi=500,
)
plt.show()


# ## Generate heatmap to visualize which channels and conditions most impact the failing image sets

# In[13]:


# Get boolean DataFrames for saturation and blur separately
sat_cols = [col for col in qc_columns if col.endswith("_Saturated")]
blur_cols = [col for col in qc_columns if col.endswith("_Blur")]

# Extract saturation and blur DataFrames
sat_df = qc_data[sat_cols]
blur_df = qc_data[blur_cols]

# Align columns so we can compare channels directly (remove suffixes for matching)
sat_df.columns = [col.replace("_Saturated", "") for col in sat_df.columns]
blur_df.columns = [col.replace("_Blur", "") for col in blur_df.columns]

# Calculate per-channel failure rates for each condition
sat_fail_rate = sat_df.mean() * 100
blur_fail_rate = blur_df.mean() * 100

# Combine into one DataFrame for plotting
heatmap_data = pd.DataFrame(
    {
        "Saturation": sat_fail_rate,
        "Blur": blur_fail_rate,
    }
)

plt.figure(figsize=(8, 5))
sns.heatmap(
    heatmap_data.T,
    annot=True,
    fmt=".1f",
    cmap="Reds",
    cbar_kws={"label": "Percent FOVs failing QC (%)"},
)
plt.title(
    "QC failure rates per channel and condition\nacross all plates in platemap layout"
)
plt.xlabel("Channel")
plt.ylabel("QC Condition")
plt.tight_layout()
plt.savefig(output_directory / "qc_failure_rates_heatmap.png", dpi=500)
plt.show()


# In[14]:


# Use blur_thresholds dict for thresholds in the plot
blur_thresholds = {
    "OrigActin": -1.8891791699802942,
    "OrigDNA": -2.2456075474546515,
    "OrigER": -2.2825812279725524,
    "OrigMito": -2.012531942517173,
    "OrigPM": -2.4309820530015642,
}


channels = ["OrigActin", "OrigDNA", "OrigER", "OrigMito", "OrigPM"]

plt.figure(figsize=(18, 4))
for i, channel in enumerate(channels, 1):
    plt.subplot(1, len(channels), i)
    sns.histplot(
        combined_df[f"ImageQuality_PowerLogLogSlope_{channel}"],
        bins=30,
        kde=True,
        color="mediumslateblue",
    )
    plt.axvline(
        blur_thresholds[channel], color="red", linestyle="--", label="Blur threshold"
    )
    plt.xlabel("PowerLogLogSlope")
    plt.title(channel)
    if i == 1:
        plt.ylabel("Count")
    else:
        plt.ylabel("")
    plt.legend([], [], frameon=False)  # Remove legend

plt.suptitle("Distribution of PowerLogLogSlope for each channel")
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()


# In[15]:


# Find all FOVs that failed DNA blur QC
failed_dna_blur = combined_df[combined_df["OrigDNA_Blur"]]

# Sort by PowerLogLogSlope for OrigDNA (ascending: more negative = blurrier)
sorted_failed = failed_dna_blur.sort_values("ImageQuality_PowerLogLogSlope_OrigDNA")

# Get the middle index (or indices) for 3 examples
n = len(sorted_failed)
mid = n // 2
indices = [mid - 1, mid, mid + 1] if n >= 3 else list(range(n))

# Select the rows at those indices
middle_examples = sorted_failed.iloc[indices][
    [
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_Site",
        "ImageQuality_PowerLogLogSlope_OrigDNA",
    ]
]
middle_examples


# In[66]:


# Find rows that are saturated for OrigMito and have no other saturation or blur failures
saturation_cols = [
    "OrigActin_Saturated",
    "OrigDNA_Saturated",
    "OrigER_Saturated",
    "OrigMito_Saturated",
    "OrigPM_Saturated",
]
blur_cols = [
    "OrigActin_Blur",
    "OrigDNA_Blur",
    "OrigER_Blur",
    "OrigMito_Blur",
    "OrigPM_Blur",
]

# mito saturated and no other saturation/blur flags
mito_only_mask = (
    combined_df["OrigMito_Saturated"]
    & ~combined_df[[c for c in saturation_cols if c != "OrigMito_Saturated"]].any(
        axis=1
    )
    & ~combined_df[blur_cols].any(axis=1)
)

mito_only = combined_df[mito_only_mask].copy()

if mito_only.empty:
    print(
        "No rows that are saturated for OrigMito only (no other saturation or blur failures)."
    )
else:
    print(f"Found {len(mito_only)} rows failing OrigMito saturation only.")
    # sample up to 3 reproducibly
    sample_n = mito_only.sample(n=min(3, len(mito_only)))
    cols = [
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_Site",
        "ImageQuality_PercentMaximal_OrigMito",
    ]
    display(sample_n[cols].reset_index(drop=True))


# In[67]:


from pathlib import Path
import matplotlib.pyplot as plt

# Base directory and platemap
images_base_dir = Path("/media/18tbdrive/CFReT_screening_data/compound_screen")
search_base = images_base_dir / platemap_name  # search within this platemap folder

exts = [".tif", ".tiff", ".png", ".jpg", ".jpeg", ".ome"]
keywords = ["d1", "mito", "origmito"]


# Helper to check if filename matches keywords (special handling for d1)
def keyword_match(name, keywords):
    name_lower = name.lower()
    for k in keywords:
        if k.lower() == "d1":
            # match if filename ends with 'd1' before the extension
            if any(name_lower.endswith(f"d1{ext}") for ext in exts):
                return True
        elif k.lower() in name_lower:
            return True
    return False


# Helper to find image for a row ensuring correct plate
def find_image_for_row_correct_plate(row, search_base, exts, keywords):
    well = str(row["Metadata_Well"]).lower()
    site = str(row["Metadata_Site"]).lower()
    plate = str(row["Metadata_Plate"]).lower()  # ensure correct plate match

    # Walk recursively through all subfolders
    for p in search_base.rglob("*"):
        if p.is_file() and p.suffix.lower() in exts:
            name = p.name.lower()
            parent_plate = p.parent.name.lower()
            # Require exact match of plate folder + well + site + mito keyword
            if (
                parent_plate == plate
                and well in name
                and site in name
                and keyword_match(name, keywords)
            ):
                return p, parent_plate
    return None, None


# Prepare list of image paths for up to 3 sampled rows
rows = sample_n.reset_index(drop=True)
image_paths = []
plates_found = []

for idx in range(min(3, len(rows))):
    row = rows.loc[idx]
    img, plate_found = find_image_for_row_correct_plate(
        row, search_base, exts, keywords
    )
    image_paths.append(img)
    plates_found.append(plate_found)

# Display found images
found = [p for p in image_paths if p is not None]
if not found:
    print("No matching Mito (d1) images found for the sampled rows.")
else:
    n = len(found)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, p in zip(axes, found):
        try:
            img = plt.imread(p)
            ax.imshow(img, cmap="gray")
        except Exception:
            ax.text(0.5, 0.5, f"Failed to read\n{p.name}", ha="center")
        ax.set_title(p.name, fontsize=9)
        ax.axis("off")
    plt.tight_layout()
    plt.show()

# Print which files and plate folders were used
for i, (p, pl) in enumerate(zip(image_paths, plates_found), 1):
    print(f"{i}: {p if p is not None else 'NOT FOUND'} (Plate folder: {pl})")


# In[71]:


import pathlib
import pandas as pd

# Load all Image.csv under qc_results, fix plate names per-file, then compute failed FOVs
qc_root = pathlib.Path("qc_results")
all_csvs = list(qc_root.rglob("Image.csv"))
print(f"Found {len(all_csvs)} Image.csv files under {qc_root}")

loaded = []
skipped = []

prefixes = ["Metadata_"]  # adjust if needed

for csv_file in all_csvs:
    try:
        df = pd.read_csv(
            csv_file, usecols=lambda c: any(c.startswith(p) for p in prefixes)
        )
    except Exception as e:
        print(f"Skipping {csv_file}: failed to read ({e})")
        skipped.append((csv_file, "read_error"))
        continue

    plate_name = csv_file.parent.name  # folder name is authoritative per-file
    if "Metadata_Plate" in df.columns and df["Metadata_Plate"].nunique() > 1:
        print(
            f"⚠️ {csv_file} had multiple plate names {df['Metadata_Plate'].unique()}, fixing to {plate_name}"
        )

    # Overwrite to ensure single unique plate name per file
    df["Metadata_Plate"] = plate_name

    # Save the platemap layout as the parent folder of the plate folder
    df["Platemap"] = str(csv_file.parent.parent.name)

    # required columns for the downstream aggregation
    required = {
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_Site",
        "Metadata_Blur_Flag",
        "Metadata_Saturation_Flag",
        "Platemap",
    }
    missing = required - set(df.columns)
    if missing:
        print(f"Skipping {csv_file}: missing required columns {missing}")
        skipped.append((csv_file, f"missing:{missing}"))
        continue

    # normalize boolean flags
    df["Metadata_Blur_Flag"] = df["Metadata_Blur_Flag"].astype(bool)
    df["Metadata_Saturation_Flag"] = df["Metadata_Saturation_Flag"].astype(bool)

    loaded.append(df)

if not loaded:
    raise RuntimeError("No valid Image.csv files were loaded.")

combined_all = pd.concat(loaded, ignore_index=True)
combined_all["Failed_Any"] = combined_all[
    ["Metadata_Blur_Flag", "Metadata_Saturation_Flag"]
].any(axis=1)

# Aggregate failed FOVs grouped by plate / well / site
failed_by_plate_well_site = (
    combined_all.groupby(
        ["Metadata_Plate", "Platemap", "Metadata_Well", "Metadata_Site"]
    )
    .agg(Total_FOVs=("Failed_Any", "size"), Failed_FOVs=("Failed_Any", "sum"))
    .reset_index()
)

# Filter for just the control wells
wells_of_interest = ["B02", "B05", "B08", "B11", "E02", "E05", "E08", "E11"]
failed_controls = failed_by_plate_well_site[
    failed_by_plate_well_site["Metadata_Well"].isin(wells_of_interest)
]

# Aggregate per plate including platemap
failed_controls_per_plate = (
    failed_controls.groupby(["Metadata_Plate", "Platemap"])
    .agg(Total_FOVs=("Total_FOVs", "sum"), Failed_FOVs=("Failed_FOVs", "sum"))
    .reset_index()
)
failed_controls_per_plate["Percent_Failed"] = (
    failed_controls_per_plate["Failed_FOVs"] / failed_controls_per_plate["Total_FOVs"]
) * 100

# Show results
print(
    "Failed FOVs for control wells (25 sites per well) per plate with platemap layout:"
)
print(failed_controls_per_plate.to_string(index=False))


# In[ ]:




