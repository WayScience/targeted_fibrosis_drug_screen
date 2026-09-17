"""Utilities for writing normalized model-fitting metadata."""

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd


FIT_SUMMARY_COLUMNS = [
    "plate",
    "fold",
    "shuffled",
    "fit_status",
    "split_seed",
    "model_random_state",
    "n_train",
    "n_test",
    "n_input_features",
    "n_selected_features",
    "average_precision",
    "roc_auc",
]

SPLIT_ASSIGNMENT_COLUMNS = [
    "plate",
    "fold",
    "split",
    "cell_type",
    "well",
]


def write_model_fit_summary(
    plate_repr: str,
    tasks: Sequence[Mapping[str, Any]],
    results: Sequence[Mapping[str, Any] | None],
    split_records: Sequence[Mapping[str, Any]],
    output_dir: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Write scalar model-fit metadata and normalized well assignments.

    ``model_fit_summary.csv`` contains one row per fitted model and only scalar
    values. ``model_fit_datasplit_assignments.csv`` contains one row per
    fold/split/cell-type/well assignment. The two files can be joined on
    ``plate`` and ``fold`` without encoding Python lists inside CSV cells.

    :param plate_repr: The plate representation string.
    :param tasks: A sequence of task dictionaries, each containing the following keys:
        - ``plate_repr``: The plate representation string.
        - ``fold``: The fold number.
        - ``shuffle_status``: The shuffle status string.
        - ``random_state``: The random state used for model fitting.
        - ``labels``: The training labels as a pandas Series.
        - ``test_labels``: The test labels as a pandas Series.
        - ``train_profiles``: The training profiles as a pandas DataFrame.
    :param results: A sequence of result dictionaries or None, each containing the following keys:
        - ``n_train``: The number of training samples.
        - ``n_test``: The number of test samples.
        - ``n_input_features``: The number of input features.
        - ``n_selected_features``: The number of selected features.
        - ``average_precision``: The average precision score.
        - ``roc_auc``: The ROC AUC score.
    :param split_records: A sequence of split record dictionaries, each containing the following keys:
        - ``fold``: The fold number.
        - ``seed``: The random seed used for the split.
        - ``train``: A dictionary mapping cell types to lists of training wells.
        - ``test``: A dictionary mapping cell types to lists of test wells.
    :param output_dir: The output directory where the CSV files will be saved.
    :return: A tuple containing two pandas DataFrames:
        - The first DataFrame contains the model fit summary.
        - The second DataFrame contains the split assignments.
    """
    if len(tasks) != len(results):
        raise ValueError(
            "Expected one model-fitting result per task, got "
            f"{len(results)} results for {len(tasks)} tasks."
        )

    split_seeds = {record["fold"]: record.get("seed") for record in split_records}
    summary_rows = [
        _model_fit_row(task, result, split_seeds.get(task["fold"]))
        for task, result in zip(tasks, results)
    ]
    fit_summary_df = pd.DataFrame(summary_rows, columns=FIT_SUMMARY_COLUMNS)

    assignment_rows = []
    for record in split_records:
        for split_name in ("train", "test"):
            assignments = record.get(split_name, {})
            if not isinstance(assignments, Mapping):
                raise TypeError(
                    f"Fold {record['fold']} {split_name!r} metadata must be a mapping."
                )
            for cell_type, wells in assignments.items():
                if isinstance(wells, (str, bytes)) or not isinstance(wells, Sequence):
                    raise TypeError(
                        f"Fold {record['fold']} {split_name!r} wells for "
                        f"{cell_type!r} must be a sequence."
                    )
                assignment_rows.extend(
                    {
                        "plate": plate_repr,
                        "fold": record["fold"],
                        "split": split_name,
                        "cell_type": cell_type,
                        "well": well,
                    }
                    for well in wells
                )

    split_assignments_df = pd.DataFrame(
        assignment_rows,
        columns=SPLIT_ASSIGNMENT_COLUMNS,
    )

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    fit_summary_df.to_csv(output_path / "model_fit_summary.csv", index=False)
    split_assignments_df.to_csv(
        output_path / "model_fit_datasplit_assignments.csv",
        index=False,
    )

    return fit_summary_df, split_assignments_df


def _model_fit_row(
    task: Mapping[str, Any],
    result: Mapping[str, Any] | None,
    split_seed: Any,
) -> dict[str, Any]:
    """
    Normalize a single model-fitting result into a row for the fit summary DataFrame.

    :param task: A task dictionary containing model-fitting information.
    :param result: A result dictionary containing model-fitting metrics, or None if the task was skipped.
    :param split_seed: The random seed used for the data split.
    :return: A dictionary representing a single row for the fit summary DataFrame.
    """
    row = dict(result or {})
    metrics_missing = result is not None and (
        pd.isna(row.get("average_precision")) or pd.isna(row.get("roc_auc"))
    )

    row.update(
        {
            "plate": task["plate_repr"],
            "fold": task["fold"],
            "shuffled": task["shuffle_status"] == "shuffled",
            "fit_status": (
                "skipped" if result is None else "failed" if metrics_missing else "success"
            ),
            "split_seed": split_seed,
            "model_random_state": task.get("random_state"),
        }
    )
    row.setdefault("n_train", len(task["labels"]))
    row.setdefault("n_test", len(task["test_labels"]))
    row.setdefault("n_input_features", task["train_profiles"].shape[1])
    row.setdefault("n_selected_features", pd.NA)
    row.setdefault("average_precision", pd.NA)
    row.setdefault("roc_auc", pd.NA)
    return row
