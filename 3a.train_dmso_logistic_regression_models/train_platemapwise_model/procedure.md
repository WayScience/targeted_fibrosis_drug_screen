## Plate-specific logistic regression model training

Plate-specific logistic regression models were trained to classify single-cell profiles as either failing or healthy. Here, “healthy” refers to cells from healthy patient heart samples, whereas “failing” refers to cells from patients with dilated cardiomyopathy. Models were trained independently for each plate and fold split.

For each plate-specific training split, features exhibiting complete or quasi-complete separation were first identified and removed when they were perfectly predictive of class label. Recursive feature elimination was then performed iteratively to reduce the feature set until the number of retained features satisfied the prespecified event-per-variable rule, using either a 1-in-10 or 1-in-20 feature-to-minority-class-sample-size criterion depending on the split and plate. A final logistic regression model was then fit using the post-RFE feature set. Model performance was evaluated on the held-out test split, and receiver operating characteristic and precision-recall curve visualizations were generated for each fitted model.

For every non-shuffled model, a corresponding shuffled-label control model was trained using the same plate, fold split, preprocessing procedure, feature selection workflow, and model-fitting procedure, except that class labels were randomly shuffled prior to training.

## Inference on treated single-cell profiles

For treatment response inference, all non-shuffled, fold-specific plate-level logistic regression models were applied to treated single-cell profiles. The inference dataset consisted of treated, non-DMSO, failing single-cell profiles that were not used during model training.

Each trained model was applied across the treated single-cell profile dataset to produce fold-specific prediction scores for each cell. These scores estimate whether each treated failing-cell profile more closely resembled cells from healthy individuals or cells from patients with dilated cardiomyopathy. The resulting k-fold prediction scores were used as single-cell morphology rescue scores in downstream treatment-response and hit-calling analyses.

## Relative cell count computation

Relative cell abundance was computed by normalizing each treated well to its plate-matched DMSO control baseline. For each plate, the median DMSO cell count was computed from control wells and used as the normalization denominator. Each treated well’s relative cell count was then calculated as:

$\text{relative cell count}_{i} = \frac{\text{treated well cell count}_i}{\text{median DMSO control cell count on the same plate}}$

This plate-level normalization accounts for plate-to-plate differences in baseline cell abundance. A value of 1 indicates that the treated well had the same cell count as the median DMSO control well on that plate, whereas values below or above 1 indicate lower or higher cell abundance relative to plate-matched controls.

For each treatment, the final relative count was computed by averaging the plate-normalized relative cell counts across all treated wells assigned to that treatment. Therefore, the reported treatment-level relative count is the mean of within-plate normalized treated-well counts, rather than a ratio computed from pooled raw counts across plates.

$\text{relative cell count}_{\text{treatment}} = \frac{1}{n} \sum_n \text{relative cell count}_{i}$

## Mean average precision computation

[To be added.]

## Hit calling

Hits were identified by applying fixed thresholds to the treatment-level summary metrics. A treatment was called as a hit if it satisfied all three criteria: mean relative cell count ≥ 0.8, mean logit score ≥ 0.9, and mean average precision ≤ 0.6. These criteria retained treatments with sufficient cell abundance, strong healthy-like morphology rescue signal, and reduced residual disease-like classification performance. The final hit set consisted only of treatments passing all three filters simultaneously.
