"""
Helper module for already split data loading and formatting to use with
    direct model training and eval. 
"""


import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler


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

    train_profiles: pd.DataFrame = _build_feature_matrix(train_data, metadata_prefix)
    test_profiles: pd.DataFrame = _build_feature_matrix(test_data, metadata_prefix)

    train_labels, test_labels = train_data[label_col], test_data[label_col]
    train_labels_shuffled = _shuffle_labels(
        train_labels, 
        random_state=shuffle_random_state
    )

    # force test columns to match train columns exactly
    test_profiles = test_profiles.reindex(columns=train_profiles.columns)

    # Drop zero-variance columns in train
    nonconstant_cols = train_profiles.columns[train_profiles.nunique(dropna=True) > 1]
    train_profiles = train_profiles.loc[:, nonconstant_cols]
    test_profiles = test_profiles.loc[:, nonconstant_cols]

    return (
        train_profiles, 
        test_profiles, 
        train_labels, 
        test_labels, 
        train_labels_shuffled
    )


def _build_feature_matrix(
    df: pd.DataFrame, 
    metadata_prefix: str="Metadata_"
) -> pd.DataFrame:
    """
    Extracts feature columns, apply standard scalar and return back as dataframe
    """

    scaler = StandardScaler().set_output(transform="pandas")

    feat_cols = [col for col in df.columns if not col.startswith(metadata_prefix)]
    df: pd.DataFrame = df.loc[:, feat_cols].apply(pd.to_numeric, errors="coerce")
    df_norm: pd.DataFrame = scaler.fit_transform(df)
    
    return df_norm


def _shuffle_labels(labels: list, random_state: int) -> list:
    """
    Shuffle labels for training null models. 
    """
    rng = np.random.default_rng(seed=random_state)
    return rng.permutation(labels)
