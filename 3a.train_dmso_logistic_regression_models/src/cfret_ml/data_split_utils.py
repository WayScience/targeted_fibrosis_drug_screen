"""
Helper module for initial data processing and splitting into train vs test sets
"""


import hashlib

import pandas as pd
import numpy as np


def string_to_int_seed(s: str, n_bytes: int = 8) -> int:
    """
    Helper hashing plate identifier to seed/salt for reproducible and varied randomization of splits
    """
    digest = hashlib.sha256(s.encode("utf-8")).digest()
    return int.from_bytes(digest[:n_bytes], "big")


def salt_seed(seed: int, salt: int) -> int:
    """
    Helper for salting the global random state with plate specific salt
    """
    # Ensure the salted seed is within the range of valid seeds for random state
    return (seed + salt) % (2**32)


def split_summary(
    df:pd.DataFrame, 
    train_index:list, 
    test_index:list,
    group_col: str = "Metadata_Well",
    label_col: str = "Metadata_cell_type",
) -> dict:
    """
    Helper for generating descriptive summary of group assignments across splits
        that is human readable and logged

    :param df: the full dataframe being split
    :param train_index: list of indices for the training split
    :param test_index: list of indices for the test split
    :param group_col: the column name in df that identifies groups (e.g. wells)
    :param label_col: the column name in df that identifies class labels (e.g. cell types)
    :return: a nested dictionary summarizing the group assignments for train and test splits
    """
    summary = {}
    
    for split_name, split_index in zip(["train", "test"], [train_index, test_index]):
        
        summary[split_name] = {
            ct: list(wells)
            for ct, wells in (
                df.iloc[split_index].groupby([group_col] + [label_col])
                .size()
                .reset_index(name="count")
            ).groupby(label_col)[group_col] 
        }

    return summary


def stratified_fold_split(
    df: pd.DataFrame,
    group_col: str = "Metadata_Well",
    class_col: str = "Metadata_cell_type",
    random_state: int | None = None,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """
    Custom datasplit function because sklearn.StratifiedGroupFold does not
    ensure representation of all classes when splitting across folds.

    This implementation produces k-fold split (more precisely k-fold heldout)
    from the input dataset, with 2 or more classes in the class label column
    and at least one group per class, keeping groups intact (not spanning
    across splits) while still ensuring representation.

    :param df: the full dataframe of metadata, with or without profiles features, 
        to be split, if training on specific condition is desired, then the input
        should be pre-filtered to only contain that condition. 
    :param group_col: the column name in df that identifies groups (e.g. wells)
    :param class_col: the column name in df that identifies class labels (e.g. cell types)
    :param random_state: optional integer seed for reproducibility of the splits
    :return: a list of tuples, where each tuple contains the train and test indices for
    """
    if group_col not in df.columns:
        raise KeyError(f"Column '{group_col}' not found in dataframe.")
    if class_col not in df.columns:
        raise KeyError(f"Column '{class_col}' not found in dataframe.")

    rng = np.random.default_rng(random_state)

    # count unique grouping labels and ensure class labels are unmixed within groups
    # in context of this project, we expect replicate wells from the same plate
    # to either contain all healthy or failing cells, never a mix of both. 
    group_label_counts = df.groupby(group_col)[class_col].nunique()
    mixed_groups = group_label_counts[group_label_counts > 1]
    if not mixed_groups.empty:
        raise ValueError(
            f"Groups must be pure (single class). Mixed groups found: "
            f"{mixed_groups.index.tolist()}"
        )

    # expect at least 2 classes and at least one group per class to be able to split
    classes = df[class_col].dropna().unique().tolist()
    if len(classes) < 2:
        raise ValueError(f"Expected at least 2 classes in {class_col}, got {len(classes)}")

    # for each class, get the unique groups that belong to that class and obtain
    # random ordering for fold assignment
    class_groups: dict[str, np.ndarray] = {}
    for label in classes:
        groups = df.loc[df[class_col] == label, group_col].unique()
        if len(groups) == 0:
            raise ValueError(f"No groups found for class '{label}' in {class_col}")
        groups = np.array(sorted(groups), dtype=object)
        rng.shuffle(groups)
        class_groups[label] = groups

    # Automatically determine the number of k fold splits based on the 
    # group replicate number from minority class. 
    # e.g. if class A has 3 groups and class B has 5 groups, this splitting
    # function will produce 3-fold split with 1 heldout group from class A and 
    # 1 heldout group from class B in each fold.
    # In our dataset, we expect 4 replicates per class from the control wells
    # so this usually results in 4-fold split. However, there can always be
    # missing wells due to viability/segmentation problems upstream so
    # we avoid hard-coding number 4 here.
    n_splits = min(len(groups) for groups in class_groups.values())
    # if group number is less than class label number, then it is highly
    # likely there is an error in metadata annotation so throw an
    # error.  
    if n_splits < 1:
        raise ValueError("Not enough groups per class to create splits.")
    
    splits: list[tuple[np.ndarray, np.ndarray]] = []
    for fold_idx in range(n_splits):
        test_groups: list[object] = []
        train_groups: list[object] = []

        for label, groups in class_groups.items():
            test_group = groups[fold_idx]
            # because this splitting is a leave one out style, testing group
            # will always be a single scalar so collection of test groups is 
            # done with append
            test_groups.append(test_group)
            # training groups will be all other groups from the same class which
            # can be one or more so the collection of training groups is done
            # with list comprehension and extend
            train_groups.extend([g for g in groups if g != test_group])

        train_idx = df[df[group_col].isin(train_groups)].index.to_numpy()
        test_idx = df[df[group_col].isin(test_groups)].index.to_numpy()
        splits.append((train_idx, test_idx))

    return splits
