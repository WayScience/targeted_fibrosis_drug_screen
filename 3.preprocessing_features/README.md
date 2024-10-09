# Preprocessing single-cell profiles

In this module, we extract single-cell profiles using CytoTable, perform single-cell quality control (QC) with coSMicQC, preprocess the features with pycytominer.
All files generated are in `parquet` file format.

## Preprocess morphology features

To perform preprocessing of single-cell morphology features, run the bash script [preprocess_features.sh](./preprocess_features.sh) using the command below:

```bash
# Make sure your current working dir is the 3.preprocessing_features folder
source preprocess_features.sh
```
