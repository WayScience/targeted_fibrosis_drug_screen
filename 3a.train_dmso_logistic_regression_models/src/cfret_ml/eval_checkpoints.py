"""
Checkpoint module for model evaluation.
The evaluation step can be more prone to stale state than other parts
    of full logit model training and evaluation orchestration.
    This module implements a checkpoint fingerprinting module.
"""

import hashlib
import json
import pathlib
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class EvaluationCheckpoints:
    """
    Checkpoints for logit model evaluation.
    Needed because previous evaluation outputs have no clear indication of what 
        specific model and data it was generated from. As the current implementation
        allows for re-running of the model fitting process with altered parameters,
        the evaluation outputs can become stale and inconsistent with the model checkpoint.
    This class enables evaluation outputs to be fingerprinted against the model checkpoint
        and data ensuring state consistency. 
    """
    metrics: pathlib.Path
    plot: pathlib.Path
    fingerprint: pathlib.Path

    @classmethod
    def for_fold(
        cls,
        output_dir: pathlib.Path,
        fold: int | str,
        shuffle_status: str,
    ) -> "EvaluationCheckpoints":
        """
        Create a new EvaluationCheckpoints instance for a given fold and shuffle status.
        """
        prefix = output_dir / f"fold_{fold}_{shuffle_status}"
        return cls(
            metrics=prefix.with_name(f"{prefix.name}_eval.json"),
            plot=prefix.with_name(f"{prefix.name}_roc_pr.png"),
            fingerprint=prefix.with_name(f"{prefix.name}_eval.sha256"),
        )

    def clear(self) -> None:
        """
        Helper for clearing all evaluation checkpoints. 
        This is useful for clean up when re-running evaluations (upon stale).
        """
        for checkpoint in (self.metrics, self.plot, self.fingerprint):
            checkpoint.unlink(missing_ok=True)

    def load_if_valid(self, expected_fingerprint: str) -> dict[str, object] | None:
        """Load the evaluation metrics if the fingerprint matches the expected value.        """
        if not all(
            checkpoint.exists()
            for checkpoint in (self.metrics, self.plot, self.fingerprint)
        ):
            return None

        try:
            if self.fingerprint.read_text().strip() != expected_fingerprint:
                return None
            metrics = json.loads(self.metrics.read_text())
            return metrics if isinstance(metrics, dict) else None
        except (OSError, json.JSONDecodeError):
            return None

    def mark_complete(self, fingerprint: str) -> None:
        """Write the completion marker after metrics and plot are saved."""
        self.fingerprint.write_text(fingerprint)


def evaluation_fingerprint(
    fitted_model: Any,
    test_profiles: pd.DataFrame,
    test_labels: pd.Series | np.ndarray,
) -> str:
    """
    Identify the fitted model and ordered data used for evaluation.
    Namely, the fingerprint is a SHA256 hash of the model parameters, test profiles, and test labels.

    :param fitted_model: The fitted model object (e.g., a statsmodels LogitResults instance).
    :param test_profiles: The test profiles DataFrame used for evaluation.
    :param test_labels: The test labels (as a Series or ndarray) used for evaluation.
    :return: A SHA256 fingerprint string representing the evaluation state.
    """
    metadata = {
        "schema": "cfret-evaluation-v1",
        "model_columns": list(fitted_model.model.exog_names),
        "test_columns": test_profiles.columns.tolist(),
    }
    fingerprint = hashlib.sha256(
        json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode()
    )

    for values in (
        fitted_model.params,
        test_profiles.to_numpy(),
        np.asarray(test_labels),
    ):
        normalized = np.ascontiguousarray(values, dtype="<f8")
        fingerprint.update(json.dumps(normalized.shape).encode())
        fingerprint.update(normalized.tobytes())

    return fingerprint.hexdigest()
