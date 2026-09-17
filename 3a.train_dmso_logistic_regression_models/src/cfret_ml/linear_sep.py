"""
Detect and repair complete separation in binary feature matrices.

Statistically does something remotely similar as the quasi-separation filter in `preprocess.py`.
Both the quasi-separation filter and this module has the same aim of dropping
    overly predictive features to stabilize ordinary least square convergence.
    The quasi-separation filter is faster but weaker as it does not fit
    full linear models. This module performs a slower but much more relevant
    filtering on reduced feature sets and hence has better stabilizing effect.
"""

from dataclasses import dataclass
from typing import TypedDict

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import linprog
from sklearn.feature_selection import RFE


@dataclass(frozen=True)
class CompleteSeparationResult:
    separable: bool
    status: int
    message: str
    coefficients: pd.Series | None = None


class SeparationRepairResult(TypedDict):
    selected_features: list[str]
    dropped_breakers: list[str]
    rescues: list[str]
    skipped_rescues: list[str]
    target_n_features: int
    final_n_features: int


def repair_rfe_separation(
    X: pd.DataFrame,
    y: pd.Series | np.ndarray,
    rfe: RFE | None,
    pre_rfe_feats: list[str] | None = None,
    target_n_features: int | None = None,
    max_rescue_trials: int = 50,
) -> SeparationRepairResult:
    """
    Remove low-importance RFE features until linear separation is broken, then refill.
    - Please see `_check_complete_separation` for more details on definition of linear separation.
    
    Having a linearly separable dataset makes the fitting of ordinary least square
        very unstable as the coefficients will arbitrarily swell towards infinity.
    This helper finds problematic features and applies minimal reduction to break separability.
    It also attempts to "refill" the feature set by adding back the lower ranking features
        ranked by recursive feature elimination so the maximal amount of features
        is preserved for modelling. 

    ``rfe`` may be ``None`` when the selector retained every input feature.
    """
    feature_names = np.asarray(pre_rfe_feats or X.columns.tolist())

    if rfe is None:
        selected = feature_names.tolist()
        importance = pd.Series(0.0, index=selected)
        ranked_candidates: list[str] = []
        default_target = len(selected)
    else:
        selected = feature_names[rfe.support_].tolist()
        importance = pd.Series(
            np.abs(rfe.estimator_.coef_).ravel(),
            index=selected,
        )
        ranked_candidates = [
            feature
            for feature, rank in sorted(
                zip(feature_names, rfe.ranking_), key=lambda item: item[1]
            )
            if rank > 1
        ]
        default_target = int(rfe.n_features_)

    target = default_target if target_n_features is None else target_n_features
    dropped_breakers: list[str] = []

    while selected and _check_complete_separation(X[selected], y).separable:
        breaker = importance.loc[selected].idxmin()
        selected.remove(breaker)
        dropped_breakers.append(breaker)

    rescues: list[str] = []
    skipped_rescues: list[str] = []
    rescue_order = list(reversed(dropped_breakers)) + ranked_candidates
    for candidate in rescue_order[:max_rescue_trials]:
        if len(selected) >= target:
            break
        if not _check_complete_separation(X[selected + [candidate]], y).separable:
            selected.append(candidate)
            rescues.append(candidate)
        else:
            skipped_rescues.append(candidate)

    return {
        "selected_features": selected,
        "dropped_breakers": dropped_breakers,
        "rescues": rescues,
        "skipped_rescues": skipped_rescues,
        "target_n_features": target,
        "final_n_features": len(selected),
    }


def find_single_feature_separation_breakers(
    X: pd.DataFrame,
    y: pd.Series | np.ndarray,
) -> list[str]:
    """Return features whose individual removal breaks complete separation."""
    if not _check_complete_separation(X, y).separable:
        return []

    return [
        feature
        for feature in X.columns
        if not _check_complete_separation(X.drop(columns=feature), y).separable
    ]


def _check_complete_separation(
    X: pd.DataFrame,
    y: pd.Series | np.ndarray,
) -> CompleteSeparationResult:
    """
    Test whether coefficients exist satisfying:

        y_i* (x_i @ beta) >= 1

    for every observation, where y_i* is in {-1, +1}.
    In other words, the dataset is completely separable if there exists a
        hyperplane (linear decision boundary in n-D) that perfectly separates
        the two classes. 
    """
    design = sm.add_constant(
        X.apply(pd.to_numeric, errors="raise"), has_constant="add"
    )
    y_values = np.asarray(y)
    if set(np.unique(y_values)) != {0, 1}:
        raise ValueError("y must contain both classes coded as 0 and 1.")
    if len(y_values) != len(design):
        raise ValueError("X and y must contain the same number of observations.")

    signed_design = np.where(y_values == 1, 1.0, -1.0)[:, None] * design.to_numpy(
        dtype=float
    )

    result = linprog(
        c=np.zeros(design.shape[1]),
        A_ub=-signed_design,
        b_ub=-np.ones(len(y_values)),
        bounds=[(None, None)] * design.shape[1],
        method="highs",
    )

    return CompleteSeparationResult(
        separable=result.success,
        status=result.status,
        message=result.message,
        coefficients=pd.Series(result.x, index=design.columns)
        if result.success
        else None,
    )
