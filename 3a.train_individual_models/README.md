# Train individual machine learning models per plate

In this module, we train a logistic regression binary classifier per plate from each batch using the DMSO control cells (healthy and failing).
We then train one classifier that combines all 4 replicate plates from the batch.
For all models, the dataset is split 70/30 for training and testing.
We apply the models to the other three plates it wasn't trained on as holdouts (except for the combined batch model).

It takes approximately **10 minutes** to train all models, including final and shuffled baseline (train datasets with features values shuffled independently) using the [`sklearn`](https://scikit-learn.org/stable/install.html) package, which utilizes CPU.
Model training is performed on a Linux-based machine running Pop_OS! LTS 22.04 with an AMD Ryzen 7 3700X 8-Core Processor, 16 CPU cores, and 64GB of RAM.

Once the models are trained, we extract precision, recall, and predicted probabilities (ranging from 0-1, 1 being 100% healthy and 0 being 100% failing).
We generate PR curves and probability distributions to evaluate the performance of the models.

## Split data, train, and evaluate models

To perform the processes described above, run the following bash [script](./train_models_and_eval.sh) in command line terminal:

```bash
source train_models_and_eval.sh
```
