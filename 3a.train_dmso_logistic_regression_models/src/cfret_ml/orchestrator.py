"""
Orchestrator module for the end to end process of data prep, cleaning, RFE, 
    final logit model fitting, and evaluation for a given fold of a given 
    plate-split-shuffle.
"""

import pathlib
import json
import joblib

import numpy as np
from scipy.linalg import LinAlgError
import pandas as pd

from .preprocess import pre_fit_selection, screen_numeric_quasi_separation
from .regression_feature_selector import fit_rfe_l1_logit_selector
from .regression_model import fit_statsmodels_logit
from .eval import evaluate_model, save_curve_plot


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
        ) = pre_fit_selection(train_profiles)
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
        else:
            print(f"\tLoaded existing fitted model for plate {plate_repr} fold {fold} {shuffle_status} from checkpoint.")
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
            print(f"\tSuccessfully fitted statsmodels Logit for plate {plate_repr} fold {fold} {shuffle_status} and saved to checkpoint.")                
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


def _empty_fold_result(
    plate_repr: str,
    fold: int | str,
    shuffle_status: str,
    labels: pd.Series | np.ndarray,
    test_labels: pd.Series | np.ndarray,
    train_profiles: pd.DataFrame,
    selected_features: list[str],
) -> dict[str, any]:
    """
    Helper invoked by the main orchestrator for making an empty result row for 
        a fold that failed to produce a valid model, with NaN for metrics and 
        counts reflecting the data and features at the point of failure.
    Model fitting can fail due to many reasons and therefore the failure
        indicating return may be needed at multiple points in the main
        orchestrator function. This helper allows the central construction
        of the empty result row to avoid code duplication and ensure consistency.
    """     
    return {
            "plate": plate_repr,
            "fold": fold,
            "shuffled": shuffle_status == "shuffled",
            "n_train": len(labels),
            "n_test": len(test_labels),
            "n_input_features": train_profiles.shape[1],
            "n_selected_features": len(selected_features),
            "average_precision": np.nan,
            "roc_auc": np.nan,
        }
