# Validation plate profiling

This module processes the validation plate (`CARD-CelIns-CX7_260803130001`) with the same workflow as modules 2, 3, and 4, using the plate map in [`metadata/platemaps/`](./metadata/platemaps/).

| Step | Folder | What it does |
|---|---|---|
| 0 | [`0.illumination_correction`](./0.illumination_correction/) | Create the LoadData file and correct illumination with CellProfiler |
| 1 | [`1.cellprofiler_processing`](./1.cellprofiler_processing/) | Segment cells and extract morphology features with CellProfiler |
| 2 | [`2.preprocessing_features`](./2.preprocessing_features/) | Convert to parquet, run single-cell QC, process single cells, and aggregate to bulk profiles |

[`metadata/`](./metadata/) holds the validation plate map and its visualization.

Run each step from its own folder with the bash script it contains:

```bash
# for example
cd 0.illumination_correction
source perform_ic.sh
```

[`2.preprocessing_features/preprocess_features.sh`](./2.preprocessing_features/preprocess_features.sh) runs conversion, QC, single-cell processing, and bulk processing in order.
Its outputs are written to `2.preprocessing_features/data/`, which git ignores.
For details on what each preprocessing step does, see the [module 4 README](../4.image_based_profiling/README.md).
