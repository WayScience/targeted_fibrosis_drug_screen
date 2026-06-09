"""
Model fitting functions for logistic regression models
"""

import pandas as pd
import statsmodels.api as sm
from statsmodels.discrete.discrete_model import Logit


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
