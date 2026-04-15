#!/usr/bin/env python
# coding: utf-8

# ## Train plate specific logistic regression models to predict failing or healthy cell status
# 
# Each training unit consumes a fold split of a plate and does the following:
# 1. Perform complete quasi-separation check and remove features perfectly predictive of label.
# 2. Perform iterative Recursive Feature Elimination (RFE) to cut down the number of features until the 1 in 10 or 20 rule is satisified per split and plate based on the minority class size.
# 3. Fit the final logistic regression model with post RFE features. 
# 4. Evaluate on testing split and produce ROC and PRC visualizations.
# 
# Per every model fitted as described, a shuffled model is also trained following the exact same procedure and same split data except with labels shuffled.   

# ## Import libraries

import json
import pathlib
from dataclasses import dataclass
from typing import Iterable, Optional, Literal
import joblib
from joblib import Parallel, delayed

import numpy as np
import pandas as pd
from pandas.api.types import (
    is_bool_dtype,
    is_categorical_dtype,
    is_object_dtype,
)
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.discrete.discrete_model import Logit
from scipy.linalg import LinAlgError
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    PrecisionRecallDisplay, 
    RocCurveDisplay,
    average_precision_score,
    roc_auc_score,
)


# ## Pathing and global parameters

random_state = 0
np.random.seed(random_state)
metadata_prefix = "Metadata_"
label_col = "Metadata_cell_type"

datasplit_dir = pathlib.Path(".") / "datasplits"
if not datasplit_dir.exists():
    raise FileNotFoundError(f"Datasplit directory not found: {datasplit_dir}")

fitted_model_dir = pathlib.Path(".") / "models"
fitted_model_dir.mkdir(exist_ok=True)

eval_plot_dir = pathlib.Path(".") / "eval_plots"
eval_plot_dir.mkdir(exist_ok=True)


# ## Helpers (a lot of them)

# ### Helpers for the quasi-separation screen step

@dataclass
class SeparationFlag:
    feature: str
    kind: str
    detail: str


def _is_categorical_like(s: pd.Series, max_unique_numeric: int = 10) -> bool:
    if is_bool_dtype(s) or is_categorical_dtype(s) or is_object_dtype(s):
        return True
    if pd.api.types.is_numeric_dtype(s):
        nunique = s.dropna().nunique()
        return nunique <= max_unique_numeric
    return False


def screen_numeric_quasi_separation(
    X: pd.DataFrame,
    y: pd.Series | list,
    x_cols: Optional[Iterable[str]] = None,
    prob_clip: float = 1e-6,
    coef_abs_threshold: float = 10.0,
    se_threshold: float = 10.0,
    maxiter: int = 100,
) -> pd.DataFrame:
    """
    Univariate logistic screen for numeric predictors.

    Flags numeric features when the one-feature logistic model shows:
      - fitted probs extremely close to 0/1
      - very large coefficient magnitude
      - very large standard error
      - fit failure / singularity / non-convergence
    """
    if x_cols is None:
        if hasattr(y, 'name'):
            x_cols = [c for c in X.columns if c != y.name]
        else:
            x_cols = list(X.columns)

    y_name = y.name if hasattr(y, 'name') else "target"
    if set(pd.Series(y).dropna().unique()) - {0, 1}:
        raise ValueError(f"{y_name} must be binary coded as 0/1.")

    rows = []

    for col in x_cols:
        x = X[col]
        if not pd.api.types.is_numeric_dtype(x):
            continue
        if _is_categorical_like(x):
            continue

        tmp = pd.DataFrame({"x": x, "y": y}).dropna()
        if tmp["x"].nunique() < 2:
            rows.append(
                {
                    "feature": col,
                    "status": "skip_constant_or_single_value",
                    "coef": np.nan,
                    "se": np.nan,
                    "min_p": np.nan,
                    "max_p": np.nan,
                    "flag": True,
                    "detail": "fewer than 2 unique non-missing values",
                }
            )
            continue

        _X = sm.add_constant(tmp["x"], has_constant="add")

        try:
            model = sm.Logit(tmp["y"], _X)
            res = model.fit(disp=False, maxiter=maxiter)

            p = res.predict(_X)
            coef = float(res.params["x"])
            se = float(res.bse["x"]) if "x" in res.bse.index else np.nan
            min_p = float(np.min(p))
            max_p = float(np.max(p))

            reasons = []
            if min_p < prob_clip or max_p > 1 - prob_clip:
                reasons.append("extreme_fitted_probabilities")
            if abs(coef) > coef_abs_threshold:
                reasons.append("large_abs_coef")
            if np.isfinite(se) and se > se_threshold:
                reasons.append("large_standard_error")
            if not bool(res.mle_retvals.get("converged", True)):
                reasons.append("nonconverged")

            rows.append(
                {
                    "feature": col,
                    "status": "ok" if not reasons else "flagged",
                    "coef": coef,
                    "se": se,
                    "min_p": min_p,
                    "max_p": max_p,
                    "flag": bool(reasons),
                    "detail": ",".join(reasons) if reasons else "",
                }
            )

        except Exception as e:
            rows.append(
                {
                    "feature": col,
                    "status": "fit_failed",
                    "coef": np.nan,
                    "se": np.nan,
                    "min_p": np.nan,
                    "max_p": np.nan,
                    "flag": True,
                    "detail": f"{type(e).__name__}: {e}",
                }
            )

    return pd.DataFrame(rows)


# ### Data-prep helpers

def compute_n_features(
    rule: Literal["one_in_ten", "one_in_twenty"],
    labels: np.ndarray,
) -> int:

    return min(np.sum(labels == 1), np.sum(labels == 0)) // (10 if rule == "one_in_ten" else 20)

def build_feature_matrix(df: pd.DataFrame, metadata_prefix: str="Metadata_") -> pd.DataFrame:
    """
    Extracts feature columns, apply standard scalar and return back as dataframe
    """

    scaler = StandardScaler().set_output(transform="pandas")

    feat_cols = [col for col in df.columns if not col.startswith(metadata_prefix)]
    df: pd.DataFrame = df.loc[:, feat_cols].apply(pd.to_numeric, errors="coerce")
    df_norm: pd.DataFrame = scaler.fit_transform(df)

    return df_norm

def split_and_prep_data(
    df: pd.DataFrame,
    train_indices: list[int],
    test_indices: list[int],
    metadata_prefix: str = "Metadata_",
    label_col: str = "Metadata_cell_type",
    shuffle_random_state: int = 42
) -> tuple:
    """
    Splits data to train and test given indices, 
        extract features and labels,
        perform basic variance checks,
        and also produce shuffled labels    
    """
    train_data = df.iloc[train_indices, :].copy()
    test_data = df.iloc[test_indices, :].copy()

    train_profiles: pd.DataFrame = build_feature_matrix(train_data, metadata_prefix)
    test_profiles: pd.DataFrame = build_feature_matrix(test_data, metadata_prefix)

    train_labels, test_labels = train_data[label_col], test_data[label_col]
    train_labels_shuffled = shuffle_labels(
        train_labels, 
        random_state=shuffle_random_state
    )

    # force test columns to match train columns exactly
    test_profiles = test_profiles.reindex(columns=train_profiles.columns)

    # Drop zero-variance columns in train
    nonconstant_cols = train_profiles.columns[train_profiles.nunique(dropna=True) > 1]
    train_profiles = train_profiles.loc[:, nonconstant_cols]
    test_profiles = test_profiles.loc[:, nonconstant_cols]

    return train_profiles, test_profiles, train_labels, test_labels, train_labels_shuffled

def shuffle_labels(labels: list, random_state: int) -> list:
    """
    We all know what this does
    """
    rng = np.random.default_rng(seed=random_state)
    return rng.permutation(labels)

def pre_fit_selection(profiles: pd.DataFrame, labels: list) -> tuple:
    """
    Basic pre-fit selection that complements the numeric quasi-separation screen
    """
    # 1. Low Variance Analysis 
    # this is now absolutely needed because I found out current cytomining 
    # feature selection may not work with continuous features. 
    variances = profiles.var()
    low_var_threshold = 1e-4
    low_var_cols = variances[variances < low_var_threshold].index

    # 2. Covariance / Colinearity Analysis
    # stricter threshold than cytomining default because regression can be
    # rather sensitive to colinearity
    corr_matrix = profiles.corr().abs()
    # Select upper triangle of correlation matrix
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    high_corr_threshold = 0.95
    high_corr_cols = [col for col in upper.columns if any(upper[col] > high_corr_threshold)]

    # Summary of columns to inspect or drop
    problematic_cols = set(low_var_cols) | set(high_corr_cols)

    return (
        low_var_cols,
        high_corr_cols,
        problematic_cols
    )


# ### RFE, final model fit and eval helpers

def fit_rfe_l1_logit_selector(
    X: pd.DataFrame,
    y: pd.Series,
    random_state: int,
    C: float = 0.1,
    max_iter: int = 5000,
    n_features_rule: Literal["one_in_ten", "one_in_twenty"] = "one_in_twenty",
    rfe_step: int = 1,
) -> tuple[pd.Index, RFE | None]:
    """
    Fits a sklearn L1 regualrized logit for feature selection purposes
    Assumes X inputs are scaled
    Assumes y inputs are binary encoded as 0 and 1
    """

    y_arr = y.to_numpy() if hasattr(y, 'to_numpy') else np.asarray(y)
    n_features = compute_n_features(n_features_rule, y_arr)
    if n_features >= X.shape[1]:
        return X.columns, None

    est = LogisticRegression(
        penalty="l1",
        solver="saga",
        random_state=random_state,
        C=C,
        max_iter=max_iter,
        class_weight="balanced",
    )

    selector = RFE(est, n_features_to_select=n_features, step=rfe_step)
    selector = selector.fit(X, y)

    return X.columns[selector.support_], selector

def fit_statsmodels_logit(
    X_train: pd.DataFrame,
    y_train: pd.Series,
) -> tuple[object, list[str]]:
    """
    Fit statsmodels Logit with an intercept.
    Returns fitted result and the actual columns used.
    """
    X_train_const = sm.add_constant(X_train, has_constant="add")

    smt_logistic = Logit(
        y_train,
        X_train_const,
        missing="drop",
        check_rank=True,
    )
    smt_result = smt_logistic.fit(disp=0)
    return smt_result, list(X_train_const.columns)

def evaluate_model(
    smt_result,
    X_test: pd.DataFrame,
    y_test: pd.Series,
) -> tuple[np.ndarray, float, float]:
    """
    Evaluate smt model on provided ds
    """
    X_test_const = sm.add_constant(X_test, has_constant="add")
    # align exactly to training columns used by the fitted statsmodels result
    X_test_const = X_test_const.reindex(columns=smt_result.model.exog_names, fill_value=0.0)

    y_score = np.asarray(smt_result.predict(X_test_const))
    y_true = y_test.to_numpy()

    ap = average_precision_score(y_true, y_score)
    roc_auc = roc_auc_score(y_true, y_score)

    return y_score, ap, roc_auc

def save_curve_plot(
    y_true: np.ndarray,
    y_score: np.ndarray,
    out_path: pathlib.Path,
    title_prefix: str,
) -> None:
    """
    Plot PRC and ROC in the same figure and save
    """
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    RocCurveDisplay.from_predictions(
        y_true,
        y_score,
        ax=ax[0],
        name="Failing",
        pos_label=1,
    )
    ax[0].set_title(f"{title_prefix} ROC Curve")
    ax[0].grid(alpha=0.3)

    PrecisionRecallDisplay.from_predictions(
        y_true,
        y_score,
        ax=ax[1],
        name="Failing",
        pos_label=1,
    )
    ax[1].set_title(f"{title_prefix} Precision-Recall Curve")
    ax[1].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ### Main orchestrator per plate-split-shuffle 

def _empty_fold_result(
    plate_repr: str,
    fold: int | str,
    shuffle_status: str,
    labels: pd.Series | np.ndarray,
    test_labels: pd.Series | np.ndarray,
    train_profiles: pd.DataFrame,
    selected_features: list[str],
) -> pd.DataFrame:
    """
    Helper invoked by the main orchestrator for making an empty result row for 
        a fold that failed to produce a valid model, with NaN for metrics and 
        counts reflecting the data and features at the point of failure.
    Model fitting can fail due to many reasons and therefore the failure
        indicating return may be needed at multiple points in the main
        orchestrator function. This helper allows the central construction
        of the empty result row to avoid code duplication and ensure consistency.
    """
    return pd.DataFrame(
        {
            "plate": [plate_repr],
            "fold": [fold],
            "shuffled": [shuffle_status == "shuffled"],
            "n_train": [len(labels)],
            "n_test": [len(test_labels)],
            "n_input_features": [train_profiles.shape[1]],
            "n_selected_features": [len(selected_features)],
            "average_precision": [np.nan],
            "roc_auc": [np.nan],
        }
    )

def process_model_fitting(
    train_profiles: pd.DataFrame, 
    test_profiles: pd.DataFrame, 
    labels: pd.Series | np.ndarray, 
    test_labels: pd.Series | np.ndarray,
    shuffle_status: str, 
    plate_repr: str, 
    fold: int | str, 
    plate_fitted_model_dir: pathlib.Path, 
    plate_eval_plot_dir: pathlib.Path, 
    random_state: int
):
    """
    The main Orchestrator function for end to end data prep, cleaning, RFE, 
        final logit model fitting, and evaluation for a given fold of a given 
        plate-split-shuffle.
    """
    # filter for low variance, covariance, and quasi-separation issues
    feat_filter_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_feat.tsv"
    if feat_filter_ckpt.exists():
        cols2keep = feat_filter_ckpt.read_text().splitlines()       
    else:
        (
            low_var_cols,
            high_corr_cols,
            cols2rm
        ) = pre_fit_selection(train_profiles, labels)
        print(f"\tPlate {plate_repr} Fold {fold} {shuffle_status} > Low var cols: {len(low_var_cols)}, High corr cols: {len(high_corr_cols)}")

        report = screen_numeric_quasi_separation(
            X=train_profiles,
            y=labels,
            coef_abs_threshold=10,
            se_threshold=10,
        )
        report.to_csv(plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_quasi_separation_report.csv", index=False)
        flagged_features = report.loc[report["flag"], "feature"].tolist()
        print(f"\tPlate {plate_repr} Fold {fold} {shuffle_status} > Flagged {len(flagged_features)} features for quasi-separation.")

        cols2rm = set(cols2rm) | set(flagged_features)
        cols2keep = [col for col in train_profiles.columns if col not in cols2rm]

        with open(feat_filter_ckpt, "w") as f:
            f.write("\n".join(cols2keep))

    if not cols2keep:
        print(f"\tNo features left after filtering in plate {plate_repr} fold {fold}, skipping model fitting.")
        return None

    train_profiles_filtered = train_profiles.loc[:, cols2keep]

    # rfe feature selection prior to final model fitting 
    rfe_selected_feat_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_rfe_selected_features.tsv"
    rfe_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_rfe.joblib"
    if rfe_selected_feat_ckpt.exists():
        selected_features = rfe_selected_feat_ckpt.read_text().splitlines()
    else:
        selected_features, rfe = fit_rfe_l1_logit_selector(
            train_profiles_filtered,
            labels,
            random_state=random_state,
            C=0.1,
            max_iter=5000,
            n_features_rule="one_in_twenty",
            rfe_step=1,
        )

        with open(rfe_ckpt, "wb") as f:
            joblib.dump(rfe, f)
        with open(rfe_selected_feat_ckpt, "w") as f:
            f.write("\n".join(selected_features))

    if len(selected_features) == 0:
        print(f"\tNo features selected by RFE in plate {plate_repr} fold {fold}, skipping model fitting.")
        return None

    train_profiles_rfe = train_profiles_filtered.loc[:, selected_features]

    # fit final statsmodels logit and save
    model_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_statsmodels_logit.joblib"
    model_fail_to_converge_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_statsmodels_logit_fit_failed.txt"

    # failure row construction here as failure can arise from multiple paths
    empty_result = _empty_fold_result(
        plate_repr=plate_repr,
        fold=fold,
        shuffle_status=shuffle_status,
        labels=labels,
        test_labels=test_labels,
        train_profiles=train_profiles,
        selected_features=selected_features,
    )

    if model_ckpt.exists():
        smt_result = joblib.load(model_ckpt)
        if smt_result is None:
            print(f"\tLoaded an empty model from {model_ckpt}, likely from a previous failed run. Treating as failure.")
            return empty_result
    elif model_fail_to_converge_ckpt.exists():
        print(f"\tPrevious attempt to fit statsmodels Logit for plate {plate_repr} fold {fold} {shuffle_status} failed to converge, skipping.")
        smt_result = None
        train_design_cols = train_profiles_rfe.columns.tolist()
         # return empty row for this fold
        return empty_result
    else:
        try:
            smt_result, train_design_cols = fit_statsmodels_logit(
                train_profiles_rfe,
                labels,
            )
            with open(plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_smt_summary.txt", "w") as f:
                f.write(str(smt_result.summary()))            
            with open(model_ckpt, "wb") as f:
                joblib.dump(smt_result, f)                
        except LinAlgError as e:
            print(f"\tStatsmodels Logit failed to fit for plate {plate_repr} fold {fold} {shuffle_status} due to LinAlgError {e}, skipping.")
            smt_result = None
            with open(model_fail_to_converge_ckpt, "w") as f:
                f.write(f"Statsmodels Logit fit failed due to LinAlgError {e}.")
            train_design_cols = train_profiles_rfe.columns.tolist()
            return empty_result
        except Exception as e:
            print(f"\tStatsmodels Logit failed to fit for plate {plate_repr} fold {fold} {shuffle_status} due to {e}, skipping")
            with open(model_fail_to_converge_ckpt, "w") as f:
                f.write(f"Statsmodels Logit fit failed due to Error {e}.")
            train_design_cols = train_profiles_rfe.columns.tolist()
            return empty_result

    # Evaluate model on test set
    eval_ckpt = plate_eval_plot_dir / f"fold_{fold}_{shuffle_status}_eval.json"
    eval_plot_ckpt = plate_eval_plot_dir / f"fold_{fold}_{shuffle_status}_roc_pr.png"
    if eval_plot_ckpt.exists():
        metric_row = json.loads(eval_ckpt.read_text())
    else: 
        y_score, ap, roc_auc = evaluate_model(
            smt_result,
            test_profiles.loc[:, selected_features],
            test_labels,
        )

        metric_row = {
            "plate": plate_repr,
            "fold": fold,
            "shuffled": shuffle_status == "shuffled",
            "n_train": len(labels),
            "n_test": len(test_labels),
            "n_input_features": train_profiles.shape[1],
            "n_selected_features": len(selected_features),
            "average_precision": ap,
            "roc_auc": roc_auc,
        }
        with open(eval_ckpt, "w") as f:
            json.dump(metric_row, f, indent=4)

        save_curve_plot(
            test_labels.to_numpy(),
            y_score,
            eval_plot_ckpt,
            title_prefix=f"{plate_repr} fold {fold} {shuffle_status}",
        )

    return metric_row


# ## Parallelized model fitting due to RFE being really slow

encoding_path = datasplit_dir / "cell_type_encoding.json"
if not encoding_path.exists():
    raise FileNotFoundError(f"Cell type encoding file not found: {encoding_path}")
encoding_dict = json.loads(encoding_path.read_text())
print(f"Loaded cell type encoding for {len(encoding_dict)} cell types.")

plate_level_splits = [
    p for p in datasplit_dir.iterdir() if p.is_dir()
]
print(f"Found {len(plate_level_splits)} plate level splits")
if not plate_level_splits:
    raise ValueError(f"No plate level splits found in {datasplit_dir}")

all_real_metrics = []
all_shuffled_metrics = []
tasks = []

for plate_dir in plate_level_splits:
    split_json_files = list(plate_dir.glob("*.json"))
    if not split_json_files:
        continue
    dmso_parquet = plate_dir / "DMSO.parquet"
    if not dmso_parquet.exists():
        continue
    dmso_df = pd.read_parquet(dmso_parquet)
    dmso_df['Metadata_cell_type'] = dmso_df['Metadata_cell_type'].map(encoding_dict)

    plate_repr = plate_dir.name
    print(f"Queueing tasks for plate {plate_repr} with {len(split_json_files)} splits")

    plate_fitted_model_dir = fitted_model_dir / plate_repr
    plate_fitted_model_dir.mkdir(exist_ok=True)

    plate_eval_plot_dir = eval_plot_dir / plate_repr
    plate_eval_plot_dir.mkdir(exist_ok=True)

    train_splits = [None] * len(split_json_files)
    test_splits = [None] * len(split_json_files)
    for split_json in split_json_files:
        with open(split_json, "r") as f:
            split_info = json.load(f)

        train_splits[split_info['fold']] = split_info["train_index"]
        test_splits[split_info['fold']] = split_info["test_index"]

    for fold_idx, (train_idx, test_idx) in enumerate(zip(train_splits, test_splits)):

        if train_idx is None or test_idx is None:
            continue

        # basic data prep
        (
            train_profiles,
            test_profiles,
            train_labels,
            test_labels,
            train_labels_shuffled,
        ) = split_and_prep_data(
            dmso_df, 
            train_idx, 
            test_idx, 
            shuffle_random_state=random_state,
            metadata_prefix=metadata_prefix,
            label_col=label_col
        )        

        n_pos = train_labels.sum()
        n_neg = len(train_labels) - n_pos
        min_class = min(n_pos, n_neg)
        if min_class <= 50 or (n_pos + n_neg <= 100):
            print(f"\tNot enough train samples in plate {plate_repr} fold {fold_idx}")
            continue

        n_pos_test = test_labels.sum()
        n_neg_test = len(test_labels) - n_pos_test
        if n_pos_test == 0 or n_neg_test == 0:
            print(f"\tTest set missing a class in plate {plate_repr} fold {fold_idx}")
            continue

        for shuffle_status, labels in zip(
            ["original", "shuffled"],
            [train_labels, train_labels_shuffled]
        ):
            tasks.append({
                "train_profiles": train_profiles,
                "test_profiles": test_profiles,
                "labels": labels,
                "test_labels": test_labels,
                "shuffle_status": shuffle_status,
                "plate_repr": plate_repr,
                "fold": fold_idx,
                "plate_fitted_model_dir": plate_fitted_model_dir,
                "plate_eval_plot_dir": plate_eval_plot_dir,
                "random_state": random_state
            })

print(f"Executing {len(tasks)} model fitting tasks in parallel (n_jobs=8)...")
results = Parallel(n_jobs=8)(delayed(process_model_fitting)(**kwargs) for kwargs in tasks)

for result in results:
    if result is not None:
        if result["shuffled"]:
            all_shuffled_metrics.append(result)
        else:
            all_real_metrics.append(result)

print(f"Finished collecting metrics: {len(all_real_metrics)} real, {len(all_shuffled_metrics)} shuffled.")

