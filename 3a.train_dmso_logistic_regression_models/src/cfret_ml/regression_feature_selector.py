"""
Module for feature selection functions for logistic regression models
"""

from typing import Literal

import pandas as pd
import numpy as np
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression


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
    n_features = _compute_n_features(n_features_rule, y_arr)
    if n_features >= X.shape[1]:
        return X.columns, None

    est = LogisticRegression(
        solver="saga",
        random_state=random_state,
        C=C,
        max_iter=max_iter,
        class_weight="balanced",
        l1_ratio=1.0,
    )

    selector = RFE(est, n_features_to_select=n_features, step=rfe_step)
    selector = selector.fit(X, y)

    return X.columns[selector.support_], selector


def _compute_n_features(
    rule: Literal["one_in_ten", "one_in_twenty"],
    labels: np.ndarray,
) -> int:
    
    return min(
        np.sum(labels == 1), np.sum(labels == 0)
    ) // (10 if rule == "one_in_ten" else 20)
