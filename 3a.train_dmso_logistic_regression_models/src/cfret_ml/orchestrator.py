"""
Orchestrator module for the end to end process of data prep, cleaning, RFE,
    final logit model fitting, and evaluation for a given fold of a given
    plate-split-shuffle.
"""

import json
import pathlib

import joblib
import numpy as np
import pandas as pd
from scipy.linalg import LinAlgError

from .eval import evaluate_model, save_curve_plot
from .eval_checkpoints import EvaluationCheckpoints, evaluation_fingerprint
from .linear_sep import repair_rfe_separation
from .preprocess import pre_fit_selection, screen_numeric_quasi_separation
from .regression_feature_selector import fit_rfe_l1_logit_selector
from .regression_model import fit_statsmodels_logit


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
    random_state: int,
    filter_linear_separation: bool = True,
    force_overwrite_post_rfe: bool = False,
):
    """
    The main Orchestrator function for end to end data prep, cleaning, RFE,
        final logit model fitting, and evaluation for a given fold of a given
        plate-split-shuffle.

    Set ``force_overwrite_post_rfe`` to preserve preprocessing and RFE
        checkpoints while regenerating all final model and evaluation outputs.
    """

    ## Step 1) filter for low variance, covariance, and quasi-separation issues
    feat_filter_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_feat.tsv"
    log_prefix = f"\tPlate {plate_repr} Fold {fold} {shuffle_status} >"
    if feat_filter_ckpt.exists():
        cols2keep = feat_filter_ckpt.read_text().splitlines()
    else:
        (
            low_var_cols,
            high_corr_cols,
            cols2rm
        ) = pre_fit_selection(train_profiles)
        print(f"{log_prefix} Low var cols: {len(low_var_cols)}, High corr cols: {len(high_corr_cols)}")

        report = screen_numeric_quasi_separation(
            X=train_profiles,
            y=labels,
            coef_abs_threshold=10,
            se_threshold=10,
        )
        report.to_csv(plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_quasi_separation_report.csv", index=False)
        flagged_features = report.loc[report["flag"], "feature"].tolist()
        print(f"{log_prefix} Flagged {len(flagged_features)} features for quasi-separation.")

        cols2rm = set(cols2rm) | set(flagged_features)
        cols2keep = [col for col in train_profiles.columns if col not in cols2rm]

        with open(feat_filter_ckpt, "w") as f:
            f.write("\n".join(cols2keep))

    if not cols2keep:
        print(f"{log_prefix} No features left after filtering, skipping model fitting.")
        return None

    train_profiles_filtered = train_profiles.loc[:, cols2keep]


    ## Step 2) RFE (recursive feature elimination) feature selection

    # This step reduces the feature size to a number proportional to the 
    # number of the minority binary classes (one-in-ten vs one-in-twenty rule),
    # which is crucial for subsequent ordinary least square model fitting.
    # (if too many features are used to model too little sample size, the model
    # struggles to find a unique solution and usually will lose statistical
    # inference power). 

    rfe_selected_feat_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_rfe_selected_features.tsv"
    rfe_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_rfe.joblib"
    if rfe_selected_feat_ckpt.exists():
        selected_features = rfe_selected_feat_ckpt.read_text().splitlines()
        rfe = joblib.load(rfe_ckpt) if rfe_ckpt.exists() else None
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
        print(f"{log_prefix} No features selected by RFE, skipping model fitting.")
        return None

    ## 2.5) force re-run after RFE

    # Because the RFE is such a time consuming step, and more robust + compared to ordinary least squares, 
    # having the ability to partially re-run the faster steps following RFE is desirable. 
    # Here all subsequent checkpoints are wiped if forced to re-run.
    
    model_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_statsmodels_logit.joblib"
    model_failure_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_statsmodels_logit_fit_failed.txt"
    model_summary_ckpt = plate_fitted_model_dir / f"fold_{fold}_{shuffle_status}_smt_summary.txt"
    eval_checkpoints = EvaluationCheckpoints.for_fold(
        plate_eval_plot_dir,
        fold,
        shuffle_status,
    )

    if force_overwrite_post_rfe:
        for checkpoint in (
            model_ckpt,
            model_failure_ckpt,
            model_summary_ckpt,
        ):
            checkpoint.unlink(missing_ok=True)
        eval_checkpoints.clear()
        print(f"{log_prefix} Cleared existing post-RFE outputs.")

    ##  3) Final linear separation repair & fit final logit model

    # The filter linear separation is almost necessary given how
    # highly correlated certain image feature combinations can be vs. outcome.

    if filter_linear_separation:
        # This step repairs linear separation issues on the RFE
        # reduced feature set, drops the least predictive feature that
        # causes linear separation, and attempts to re-fill
        # the feature set with "wait-listed" (previously eliminated) features.

        repair = repair_rfe_separation(
            train_profiles_filtered,
            labels,
            rfe,
            pre_rfe_feats=train_profiles_filtered.columns.tolist()
            if rfe is not None
            else selected_features,
            target_n_features=len(selected_features),
        )
        selected_features = repair["selected_features"]
        print(
            f"{log_prefix} Separation repair summary:\n"
            f"\t  Target features: {repair['target_n_features']}\n"
            f"\t  Final features: {repair['final_n_features']}\n"
            f"\t  Dropped breakers: {len(repair['dropped_breakers'])}\n"
            f"\t  Rescued features: {len(repair['rescues'])}"
        )

    train_profiles_rfe = train_profiles_filtered.loc[:, selected_features]

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

    smt_result = None
    if model_ckpt.exists():
        smt_result = joblib.load(model_ckpt)
        if smt_result is None:
            print(f"\tLoaded an empty model from {model_ckpt}, likely from a previous failed run. Treating as failure.")
            return empty_result
        if list(smt_result.model.exog_names[1:]) == selected_features:
            print(f"\tLoaded existing fitted model for plate {plate_repr} fold {fold} {shuffle_status} from checkpoint.")
        else:
            print(f"{log_prefix} Existing model uses different features; refitting.")
            smt_result = None
    elif model_failure_ckpt.exists():
        print(f"\tPrevious attempt to fit statsmodels Logit for plate {plate_repr} fold {fold} {shuffle_status} failed to converge, skipping.")
        return empty_result

    # Fit final logit model using the repaired feautres if enabled. 
    # This can fail in many varied ways, the most common path of failure
    # due to perfect linear separation is addressed by the separation repair.

    if smt_result is None:
        try:
            smt_result, _ = fit_statsmodels_logit(
                train_profiles_rfe,
                labels,
            )
            with open(model_summary_ckpt, "w") as f:
                f.write(str(smt_result.summary()))
            with open(model_ckpt, "wb") as f:
                joblib.dump(smt_result, f)
            print(f"\tSuccessfully fitted statsmodels Logit for plate {plate_repr} fold {fold} {shuffle_status} and saved to checkpoint.")
        except LinAlgError as e:
            print(f"\tStatsmodels Logit failed to fit for plate {plate_repr} fold {fold} {shuffle_status} due to LinAlgError {e}, skipping.")
            with open(model_failure_ckpt, "w") as f:
                f.write(f"Statsmodels Logit fit failed due to LinAlgError {e}.")
            return empty_result
        except Exception as e:
            print(f"\tStatsmodels Logit failed to fit for plate {plate_repr} fold {fold} {shuffle_status} due to {e}, skipping")
            with open(model_failure_ckpt, "w") as f:
                f.write(f"Statsmodels Logit fit failed due to Error {e}.")
            return empty_result

    ## 4) Evaluate model on test set

    # Because the eval step has the least amount of clues on what has happened
    # and changed in previous step. It has a dedicated finger printing step
    # to validate whether its checkpointed results are outdated relative to the
    # full profile set and features. 

    test_profiles_selected = test_profiles.loc[:, selected_features]
    eval_fingerprint = evaluation_fingerprint(
        smt_result,
        test_profiles_selected,
        test_labels,
    )
    metric_row = eval_checkpoints.load_if_valid(eval_fingerprint)

    if metric_row is None:
        eval_checkpoints.clear()

        y_score, ap, roc_auc = evaluate_model(
            smt_result,
            test_profiles_selected,
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
        with open(eval_checkpoints.metrics, "w") as f:
            json.dump(metric_row, f, indent=4)

        save_curve_plot(
            test_labels.to_numpy(),
            y_score,
            eval_checkpoints.plot,
            title_prefix=f"{plate_repr} fold {fold} {shuffle_status}",
        )
        eval_checkpoints.mark_complete(eval_fingerprint)

    return metric_row


def _empty_fold_result(
    plate_repr: str,
    fold: int | str,
    shuffle_status: str,
    labels: pd.Series | np.ndarray,
    test_labels: pd.Series | np.ndarray,
    train_profiles: pd.DataFrame,
    selected_features: list[str],
) -> dict[str, object]:
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
