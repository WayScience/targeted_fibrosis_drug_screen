#!/usr/bin/env python
# coding: utf-8

# # Generate IC pipeline per platemap layout with updated QC thresholds

# In[1]:


import pathlib
import json
import re


# In[2]:


threshold_dir = pathlib.Path("../1a.whole_image_qc/blur_thresholds")
blur_thresholds = {}

# Load all platemap JSONs into a dictionary
for i in range(1, 12):  # 1 through 11
    file_path = threshold_dir / f"blur_thresholds_platemap_{i}.json"
    with open(file_path, "r") as f:
        blur_thresholds[i] = json.load(f)


# In[3]:


# Directory to save new pipeline files
pipeline_dir = pathlib.Path("./pipeline/updated_pipelines")
pipeline_dir.mkdir(exist_ok=True)

# Load in the base illumination pipeline content
base_illum_file = pathlib.Path("./pipeline/base_illum.cppipe")
with open(base_illum_file, "r") as f:
    base_illum_content = f.read()

# Updated pattern — handles indentation and flexible formatting
pattern_template = (
    r"(^\s*Which measurement\?\s*:\s*ImageQuality_PowerLogLogSlope_{}"
    r"\s*\n\s*Flag images based on low values\?\s*:\s*Yes"
    r"\s*\n\s*Minimum value\s*:\s*)(-?\d+\.?\d*)"
)

for platemap_num, thresholds in blur_thresholds.items():
    new_content = base_illum_content

    for channel, min_val in thresholds.items():
        new_content, n_replacements = re.subn(
            pattern_template.format(re.escape(channel)),
            r"\g<1>{}".format(min_val),
            new_content,
            flags=re.MULTILINE,
        )

        if n_replacements == 0:
            print(f"⚠️ No match found for channel {channel} in platemap {platemap_num}")

    # Save updated pipeline file
    new_file = pipeline_dir / f"illum_platemap_{platemap_num}.cppipe"
    with open(new_file, "w") as f:
        f.write(new_content)

    print(f"✅ Saved updated pipeline for platemap {platemap_num}")

