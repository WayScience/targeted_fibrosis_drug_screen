# Whole image quality control

In this module, we run the whole image QC CellProfiler pipeline on all 44 plates and summarize which image sets pass.
Poor-quality images (over-saturated or blurry) are flagged here, and [module 2](../2.illumination_correction/) skips the flagged image sets during illumination correction.

| Step | Notebook | What it does |
|---|---|---|
| 0 | [`0.run_image_qc.ipynb`](./0.run_image_qc.ipynb) | Run the CellProfiler QC pipeline in parallel on every plate |
| 1 | [`1.find_qc_thresholds.ipynb`](./1.find_qc_thresholds.ipynb) | Derive blur and saturation thresholds (`blur_thresholds/`) |
| 2 | [`2.whole_screen_results.ipynb`](./2.whole_screen_results.ipynb) | Summarize QC results across the screen (`qc_plots/`) |

Run all steps with the CellProfiler environment (`fibrosis_cp_env`):

```bash
# Make sure your current working dir is the 1.whole_image_qc folder
source perform_image_qc.sh
```

The pipeline is in [`pipeline/`](./pipeline/).
