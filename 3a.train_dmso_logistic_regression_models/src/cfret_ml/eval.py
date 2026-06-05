"""
Evaluation functions for ML classifiers
"""

import pathlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    RocCurveDisplay,
    PrecisionRecallDisplay,
)
import statsmodels.api as sm
from statsmodels.discrete.discrete_model import Logit


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
