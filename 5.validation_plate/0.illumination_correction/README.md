# Illumination correction & whole image quality control (QC)

In this module, we apply a modified version of the CellProfiler illumination correction pipeline from the [cellpainting_predicts_cardiac_fibroblasts](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis) repository.

## CellProfiler pipeline

The pipeline that we will be utilizing performs two tasks:

1. Whole image quality control (QC) -> skip IC processing on an image set if it is too blur or over-saturated.
2. Correct and save images for illumination errors

## Perform IC on data

To calculate, apply, and save images that have been corrected, run the bash script [perform_ic.sh](./perform_ic.sh) using the command below:

```bash
source perform_ic.sh
```
