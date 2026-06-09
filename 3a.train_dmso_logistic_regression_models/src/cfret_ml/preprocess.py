"""
Preprocessing functions for logistic regression models
"""

from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from pandas.api.types import (
    is_bool_dtype,
    is_categorical_dtype,
    is_object_dtype,
)
import statsmodels.api as sm


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


def pre_fit_selection(profiles: pd.DataFrame) -> tuple:
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
