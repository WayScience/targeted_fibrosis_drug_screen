# Train pooled machine learning models per platemap

In this module, we train platemap specific logistic regression binary classifiers (11 platemaps * 2 models (final/shuffled) * fold split = 22 models * fold split). 
Each platemap design has 4 replicate plates, 
which in this analysis module we pool together to produce training data. 
We train models the same way as we train the plate specific models, 
using DMSO control cells predicting the heart type of origin, 
either healthy or failing (`cell_type` column in metadata).

We apply a hold one plate out fold split strategy per platemap. 
This helps assess potential cross plate generalizability, 
and help with catching image batch effect led erratic model performance instability. 

## Running

TBA
