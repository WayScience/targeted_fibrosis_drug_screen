# Train individual machine learning models per plate

In this module, we train plate specific logistic regression binary classifiers (44 plates * 2 models (final/shuffled) * fold split = 88 models * fold split) on DMSO control wells, 
predicting from the profiles the patient derivation cell type label(healthy vs. failing).

We apply a fold split strategy per plate to leave one whole well of healthy and a whole well of failing cells out to assess potential cross plate generalizability, 
and help with catching image batch effect led model performance instability. 

## Running

TBA
