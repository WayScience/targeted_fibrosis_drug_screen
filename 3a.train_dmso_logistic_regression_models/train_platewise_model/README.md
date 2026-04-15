# Train individual machine learning models per plate

In this module, we train the plate specific logistic regression binary classifier 
models using the DMSO control cells predicting from the profiles the cell type label(healthy and failing).

We apply a pseudo-stratified fold split strategy per plate that consistently
divides a whole well of healthy and a whole well of failing cells into the
testing set, and leaves the rest for training. 
This is to ensure single cells from the same well cannot span train and test split which is one major source of erratic and hard to explain model performance drop due to problematic well(s). 

## Running

TBA
