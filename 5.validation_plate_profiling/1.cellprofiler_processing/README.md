# Segmentation and feature extraction

In this module, we apply a modified version of the CellProfiler analysis pipeline from the [cellpainting_predicts_cardiac_fibroblasts](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis) repository.

## CellProfiler pipeline

The pipeline that we will be utilizing performs two tasks:

1. Segmentation for cellular compartments (nuclei, cells, and cytoplasm)
2. Extract features from compartments per channel.

## Perform analysis on data

To segment and extract features from the plate, run the bash script [perform_extraction.sh](./perform_extraction.sh) using the command below:

```bash
source perform_extraction.sh
```
