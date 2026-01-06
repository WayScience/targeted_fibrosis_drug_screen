# CellProfiler segmentation and feature extraction

In this module, we apply a modified version of the CellProfiler processing/analysis pipeline from the [cellpainting_predicts_cardiac_fibroblasts](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis) repository.
After careful manual evaluation of the segmentation parameters on this dataset, we decided to update the parameters in the dataset after downloading slightly to get better segmentation results.

## Perform analysis on data

To perform segmentation and extract morphological features, run the bash script [perform_processing.sh](./perform_processing.sh) using the command below:

```bash
# Make sure your current working dir is the 2.cellprofiler_processing folder
source perform_processing.sh
```

**It took approximately 27 hours to run segmentation and feature extraction on 4 plates at the same time with ~1,100 image sets (group of channels per FOV) per plate using a Linux-based machine running Pop_OS! LTS 22.04 with an AMD Ryzen 7 3700X 8-Core Processor.**
**There is a total of 16 CPUs with 125 GB of MEM.**
