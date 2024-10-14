# Assess fibrosis reversal

In this module, we will apply a pre-trained binary logistic regression model trained on cardiac fibroblasts which predicts if a cell's morphology is from a nonfailing patient (healthy) or failing patient (failing) on a scale of 0 to 1.
This model comes from the [cellpainting_predicts_cardiac_fibroblasts](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis) repository.
We download the final and shuffled baseline models in the first notebook.
This is where we apply the models each plate and output a parquet file per plate containing the relevant metadata and predicted probability values for healthy and failing per single-cell.

## Apply models and evaluate fibrosis reversal

To perform apply the models to the data and evaluate, run the bash script [apply_model.sh](./apply_model.sh) using the command below:

```bash
# Make sure your current working dir is the 4.assess_fibrosis_reversal folder
source apply_model.sh
```
