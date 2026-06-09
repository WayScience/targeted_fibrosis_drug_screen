# Train individual machine learning models per plate

In this module, we train plate specific logistic regression binary classifiers (44 plates * 2 models (final/shuffled) * fold split = 88 models * fold split).
We train models using DMSO control cells predicting the heart type of origin, 
either healthy or failing (`cell_type` column in metadata).

We apply a fold split strategy per plate. 
This method produces minority class well level replicate number folds each leaving out a different combination of whole well of healthy and a whole well of failing cells. 
This ensures single cells from the same well cannot span both the train and test split during any model training, 
which helps assess potential cross plate generalizability, 
and help with catching image batch effect led erratic model performance instability. 

## Running

TBA
