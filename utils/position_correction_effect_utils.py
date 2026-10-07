"""Study how the plate-position correction changes the feature space.

``position_correction_utils`` estimates the plate-position effect ("tilt") and corrects
single cells with it. The correction of a well is ``amplitude * tilt(row, column)``, and
it shifts every cell of the well by the same vector. The median profile of a well after
the correction is therefore the median profile before the correction minus that vector.
The functions here use this to compare wells before and after the correction, without
reading single cells.
"""

from collections.abc import Sequence

import numpy as np
import pandas as pd
import position_correction_utils as pcu
from sklearn.decomposition import PCA

CHANNELS = ("Actin", "DNA", "ER", "Mito", "PM")
MULTI_CHANNEL = "Multiple channels"
NO_CHANNEL = "No channel"


def correct_well_medians(medians: pd.DataFrame, fit: dict) -> pd.DataFrame:
    """Apply the position correction to well median profiles.

    Parameters
    ----------
    medians : pd.DataFrame
        One row per well, with the columns ``platemap`` (number), ``well`` and the
        features of the fit.
    fit : dict
        Position correction fit (see ``position_correction_utils.load_fit``).

    Returns
    -------
    pd.DataFrame
        Copy of ``medians`` where each feature of the fit has the correction
        ``amplitude * tilt(row, column)`` of the platemap of the well subtracted.

    Raises
    ------
    KeyError
        If a platemap of the wells has no amplitude in the fit.
    """
    rows, cols = pcu.well_index(medians["well"])
    amplitudes = medians["platemap"].map(fit["amplitudes"])
    if amplitudes.isna().any():
        missing = sorted(medians.loc[amplitudes.isna(), "platemap"].unique())
        raise KeyError(f"No amplitude in the fit for platemaps {missing}")
    offsets = amplitudes.to_numpy()[:, None] * pcu.tilt(
        fit["coef"], fit["center"], rows, cols
    )
    corrected = medians.copy()
    corrected[fit["features"]] = (
        medians[fit["features"]].to_numpy(dtype=float) - offsets
    )
    return corrected


def parse_feature_family(feature: str) -> dict[str, str]:
    """Split a CellProfiler feature name into compartment, measurement and channel.

    Parameters
    ----------
    feature : str
        Feature name such as ``"Cells_Texture_SumAverage_DNA_3_00_256"``. The first
        part is the compartment and the second part is the measurement type.

    Returns
    -------
    dict[str, str]
        A dictionary with the keys:

        - "compartment": for example ``"Cells"``
        - "measurement": for example ``"Texture"``
        - "channel": the channel (``"DNA"``), ``MULTI_CHANNEL`` for features that
          combine channels (such as correlations), or ``NO_CHANNEL`` for features
          without a channel (such as area and shape)
    """
    parts = feature.split("_")
    channels = []
    for part in parts[2:]:
        if part in CHANNELS and part not in channels:
            channels.append(part)
    if not channels:
        channel = NO_CHANNEL
    elif len(channels) > 1:
        channel = MULTI_CHANNEL
    else:
        channel = channels[0]
    return {"compartment": parts[0], "measurement": parts[1], "channel": channel}

def tilt_size_by_feature(fit: dict) -> pd.DataFrame:
    """Measure the size of the tilt of every feature and label its family.

    The size of a feature is the root mean square of its tilt over all positions of the
    plate (6 rows by 10 columns), with an amplitude of 1. The features are z-scored on
    each plate, so the size is in units of the standard deviation of a single cell.

    Parameters
    ----------
    fit : dict
        Position correction fit (see ``position_correction_utils.load_fit``).

    Returns
    -------
    pd.DataFrame
        One row per feature with the columns ``feature``, ``tilt_rms``,
        ``compartment``, ``measurement`` and ``channel``.
    """
    rows, cols = np.meshgrid(
        np.arange(pcu.N_ROW_TERMS), np.arange(pcu.N_COL_TERMS + 1), indexing="ij"
    )
    tilt = pcu.tilt(fit["coef"], fit["center"], rows.ravel(), cols.ravel())
    sizes = pd.DataFrame(
        {
            "feature": fit["features"],
            "tilt_rms": np.sqrt((tilt**2).mean(axis=0)),
        }
    )
    families = pd.DataFrame([parse_feature_family(name) for name in sizes["feature"]])
    return pd.concat([sizes, families], axis=1)

def fit_pca(profiles: np.ndarray, n_components: int = 10, random_state: int = 0) -> PCA:
    """Fit a PCA on well profiles without scaling the features.

    The features are z-scored on each plate, so they are already on a similar scale.

    Parameters
    ----------
    profiles : np.ndarray
        Wells by features matrix.
    n_components : int, optional
        Number of principal components, by default 10.
    random_state : int, optional
        Seed of the PCA, by default 0.

    Returns
    -------
    PCA
        Fitted PCA. ``transform`` centers new profiles with the mean of the fitted ones.
    """
    return PCA(n_components=n_components, random_state=random_state).fit(profiles)

def variance_explained_by(values: np.ndarray, groups: Sequence) -> float:
    """Fraction of the variance of a score that the groups explain.

    This is the between-group sum of squares divided by the total sum of squares (eta
    squared) of a one-way analysis of variance.

    Parameters
    ----------
    values : np.ndarray
        Score of each well, for example its first principal component.
    groups : Sequence
        Group of each well, for example its plate row.

    Returns
    -------
    float
        Fraction between 0 and 1. It is 0 when the variance of ``values`` is 0.
    """
    values = np.asarray(values, dtype=float)
    total = ((values - values.mean()) ** 2).sum()
    if total == 0:
        return 0.0
    group_means = pd.Series(values).groupby(np.asarray(groups)).transform("mean")
    return float(((group_means.to_numpy() - values.mean()) ** 2).sum() / total)

def well_score_table(
    medians: pd.DataFrame, before: np.ndarray, after: np.ndarray
) -> pd.DataFrame:
    """Combine the principal component scores of wells before and after correction.

    Parameters
    ----------
    medians : pd.DataFrame
        One row per well, with the columns ``plate``, ``platemap`` (number), ``well``,
        ``treatment`` and ``cell_type``.
    before, after : np.ndarray
        Wells by components scores before and after the correction, in the same
        principal components and in the row order of ``medians``.

    Returns
    -------
    pd.DataFrame
        Two rows per well, one for each ``version`` (``"Before correction"`` and
        ``"After correction"``), with the columns ``plate``, ``platemap_number``,
        ``well``, ``well_row``, ``well_col``, ``treatment``, ``cell_type``,
        ``is_control`` (DMSO wells), ``version`` and one column ``PC1``, ``PC2``, ...
        for each component.
    """
    wells = medians["well"].astype(str).reset_index(drop=True)
    info = pd.DataFrame(
        {
            "plate": medians["plate"].to_numpy(),
            "platemap_number": medians["platemap"].astype(int).to_numpy(),
            "well": wells,
            "well_row": wells.str[0],
            "well_col": wells.str[1:].astype(int),
            "treatment": medians["treatment"].to_numpy(),
            "cell_type": medians["cell_type"].to_numpy(),
            "is_control": (medians["treatment"] == "DMSO").to_numpy(),
        }
    )
    components = [f"PC{i + 1}" for i in range(before.shape[1])]
    tables = [
        pd.concat([info, pd.DataFrame(scores, columns=components)], axis=1).assign(
            version=version
        )
        for version, scores in [
            ("Before correction", before),
            ("After correction", after),
        ]
    ]
    return pd.concat(tables, ignore_index=True)
