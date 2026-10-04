"""
This module provides functions to estimate and remove plate-position ("tilt") effects
from single-cell profiles.

Wells in rows B to G and columns 2 to 11 are shifted in a consistent, feature-specific
direction depending on where they sit on the plate. Because compounds were placed at
random, their own effects average out when we pool many compounds across platemaps, and
what remains is the position effect. We estimate this effect once as a shared "tilt map",
scale it per platemap (with shrinkage toward the screen-wide average), and subtract it
from every cell according to the cell's well position.

The profiles must be normalized (per plate) before correction so that all plates share
the same feature space.
"""

import pathlib

import numpy as np
import pandas as pd

ROWS = "BCDEFG"
CONTROL_TREATMENTS = ("DMSO", "TGFRi")
N_ROW_TERMS = len(ROWS)
N_COL_TERMS = 9  # columns 3 to 11, column 2 is the baseline


def well_index(
    wells: pd.Series | np.ndarray | list[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Convert well names to zero-based row and column indices.

    Parameters
    ----------
    wells : pd.Series | np.ndarray | list[str]
        Well names such as "B04".

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Row indices (B to G become 0 to 5) and column indices (2 to 11 become 0 to 9).
    """
    wells = np.asarray(wells)
    rows = np.array([ROWS.index(well[0]) for well in wells])
    cols = np.array([int(well[1:]) - 2 for well in wells])
    return rows, cols


def discover_normalized_profiles(base_dir: str | pathlib.Path) -> pd.DataFrame:
    """Find the normalized single-cell profile of every plate.

    Expects the layout ``<base_dir>/batch_*/platemap_*/single_cell_profiles/`` containing
    ``<plate>_sc_normalized.parquet`` files.

    Parameters
    ----------
    base_dir : str | pathlib.Path
        Directory containing the batch folders.

    Returns
    -------
    pd.DataFrame
        One row per plate single cell profile file with columns "plate", "platemap" (int) and "path".
    """
    records = []
    for path in sorted(
        pathlib.Path(base_dir).glob(
            "batch_*/platemap_*/single_cell_profiles/*_sc_normalized.parquet"
        )
    ):
        records.append(
            {
                "plate": path.name.removesuffix("_sc_normalized.parquet"),
                "platemap": int(path.parents[1].name.split("_")[1]),
                "path": path,
            }
        )
    return pd.DataFrame(records)


def _plate_well_medians(path: pathlib.Path, plate: str, platemap: int) -> pd.DataFrame:
    """Median profile of every well on one plate."""
    cells = pd.read_parquet(path)
    features = [c for c in cells.columns if not c.startswith("Metadata_")]
    grouped = cells.groupby("Metadata_Well", sort=True)
    table = grouped[features].median().astype("float32")
    table.insert(0, "n", grouped.size())
    table.insert(0, "cell_type", grouped["Metadata_cell_type"].first())
    table.insert(0, "treatment", grouped["Metadata_treatment"].first())
    table.insert(0, "platemap", platemap)
    table.insert(0, "plate", plate)
    return table.rename_axis("well").reset_index()


def build_well_medians(
    plate_table: pd.DataFrame, n_jobs: int = 1
) -> tuple[pd.DataFrame, list[str]]:
    """Summarize every well as the median profile of its cells.

    Parameters
    ----------
    plate_table : pd.DataFrame
        Output of ``discover_normalized_profiles``.
    n_jobs : int, optional
        Number of plates to read in parallel (requires joblib when not 1), by default 1.

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        Well-level table (one row per plate and well) with columns "plate",
        "platemap", "well", "treatment", "cell_type", "n" (cells) and the features,
        and the list of features. Features with a missing value in any well are dropped.
    """
    jobs = list(plate_table[["path", "plate", "platemap"]].itertuples(index=False))
    if n_jobs == 1:
        tables = [_plate_well_medians(*job) for job in jobs]
    else:
        from joblib import Parallel, delayed

        tables = Parallel(n_jobs=n_jobs)(
            delayed(_plate_well_medians)(*job) for job in jobs
        )
    well_medians = pd.concat(tables, ignore_index=True)
    meta = ["plate", "platemap", "well", "treatment", "cell_type", "n"]
    features = [c for c in well_medians.columns if c not in meta]
    features = [c for c in features if not well_medians[c].isna().any()]
    return well_medians[meta + features], features


def compound_table(
    well_medians: pd.DataFrame, features: list[str], min_cells: int = 50
) -> tuple[np.ndarray, pd.DataFrame]:
    """Build one profile per compound well, averaged over replicate plates.

    Each platemap is centered on its own mean so that overall platemap level does not
    enter the position estimate.

    Parameters
    ----------
    well_medians : pd.DataFrame
        Output of ``build_well_medians``.
    features : list[str]
        Feature columns to use.
    min_cells : int, optional
        Minimum cells for a well to be used, by default 50.

    Returns
    -------
    tuple[np.ndarray, pd.DataFrame]
        Compound wells by features matrix, and a table with "platemap", "well", "row"
        and "col" describing its rows.
    """
    wells = well_medians[
        ~well_medians["treatment"].isin(CONTROL_TREATMENTS)
        & (well_medians["n"] >= min_cells)
    ]
    profiles = wells.groupby(["platemap", "well"])[features].mean()
    profiles = profiles - profiles.groupby(level="platemap").transform("mean")
    meta = profiles.index.to_frame(index=False)
    meta["row"], meta["col"] = well_index(meta["well"])
    return profiles.to_numpy(float), meta


def design_matrix(rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    """Indicator design with six row terms and nine column terms."""
    design = np.zeros((len(rows), N_ROW_TERMS + N_COL_TERMS))
    design[np.arange(len(rows)), rows] = 1
    later_cols = cols > 0
    design[np.where(later_cols)[0], N_ROW_TERMS + cols[later_cols] - 1] = 1
    return design


def fit_tilt_map(
    profiles: np.ndarray, rows: np.ndarray, cols: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Fit the shared position effect with a robust additive row plus column model.

    Residuals are clipped at three scaled MADs per feature and the model is refit once
    so that a few strongly active compounds do not drive the estimate.

    Parameters
    ----------
    profiles : np.ndarray
        Compound wells by features, centered per platemap.
    rows, cols : np.ndarray
        Zero-based row and column index of each compound well.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Coefficients (15 by features) and the per-feature mean of the fitted values,
        which centers the tilt at zero over the training wells.
    """
    design = design_matrix(rows, cols)
    coef = np.linalg.lstsq(design, profiles, rcond=None)[0]
    resid = profiles - design @ coef
    mad = 1.4826 * np.median(np.abs(resid - np.median(resid, axis=0)), axis=0) + 1e-9
    winsorized = design @ coef + np.clip(resid, -3 * mad, 3 * mad)
    coef = np.linalg.lstsq(design, winsorized, rcond=None)[0]
    return coef, (design @ coef).mean(axis=0)


def tilt(
    coef: np.ndarray, center: np.ndarray, rows: np.ndarray, cols: np.ndarray
) -> np.ndarray:
    """Feature-space offset for each well position (wells by features)."""
    return design_matrix(rows, cols) @ coef - center


def amplitude(profiles: np.ndarray, tilt_values: np.ndarray) -> float:
    """Least-squares scale (through the origin) of profiles onto the tilt."""
    return float((profiles * tilt_values).sum() / (tilt_values * tilt_values).sum())


def platemap_amplitudes(
    profiles: np.ndarray,
    meta: pd.DataFrame,
    coef: np.ndarray,
    center: np.ndarray,
    n_boot: int = 200,
    seed: int = 0,
) -> dict[int, tuple[float, float]]:
    """Estimate how strongly the tilt applies to each platemap.

    Parameters
    ----------
    profiles, meta : output of ``compound_table``
    coef, center : output of ``fit_tilt_map``
    n_boot : int, optional
        Bootstrap resamples of wells for the standard error, by default 200.
    seed : int, optional
        Random seed, by default 0.

    Returns
    -------
    dict[int, tuple[float, float]]
        Platemap number mapped to (amplitude, bootstrap standard error).
    """
    rng = np.random.default_rng(seed)
    result = {}
    for platemap in sorted(meta["platemap"].unique()):
        keep = (meta["platemap"] == platemap).to_numpy()
        tilt_values = tilt(
            coef, center, meta["row"].to_numpy()[keep], meta["col"].to_numpy()[keep]
        )
        subset = profiles[keep]
        draws = (rng.integers(0, keep.sum(), keep.sum()) for _ in range(n_boot))
        boots = [amplitude(subset[i], tilt_values[i]) for i in draws]
        result[platemap] = (amplitude(subset, tilt_values), float(np.std(boots)))
    return result


def shrink_amplitude(amplitudes: dict[int, tuple[float, float]], target: int) -> float:
    """Shrink a platemap's amplitude toward the mean of the other platemaps.

    The weight on the platemap's own estimate is tau2 / (tau2 + se^2), where tau2 is the
    platemap-to-platemap variance left after removing sampling error and se is the
    platemap's standard error.

    Parameters
    ----------
    amplitudes : dict[int, tuple[float, float]]
        Output of ``platemap_amplitudes``.
    target : int
        Platemap to shrink.

    Returns
    -------
    float
        Shrunk amplitude.
    """
    others = [v for platemap, v in amplitudes.items() if platemap != target]
    amps = np.array([v[0] for v in others])
    ses = np.array([v[1] for v in others])
    tau2 = max(0.0, amps.var(ddof=1) - np.mean(ses**2))
    mean = np.average(amps, weights=1 / (ses**2 + tau2))
    own, own_se = amplitudes[target]
    weight = tau2 / (tau2 + own_se**2) if tau2 + own_se**2 > 0 else 0.0
    return float(mean + weight * (own - mean))


def fit_position_correction(
    well_medians: pd.DataFrame,
    features: list[str],
    hold_out_platemap: int | None = None,
) -> dict:
    """Fit the tilt map and per-platemap amplitudes.

    Parameters
    ----------
    well_medians : pd.DataFrame
        Output of ``build_well_medians``.
    features : list[str]
        Features to model.
    hold_out_platemap : int | None, optional
        If set, the tilt map is fit without this platemap (useful to evaluate the
        correction on data it did not see). By default all platemaps are used.

    Returns
    -------
    dict
        "coef", "center" (tilt map), "features", "amplitudes" (platemap mapped to its
        shrunk amplitude; only the held-out platemap if one is given) and
        "raw_amplitudes" (platemap mapped to the unshrunk amplitude and its standard
        error, for every platemap).
    """
    profiles, meta = compound_table(well_medians, features)
    fit = (
        np.ones(len(meta), dtype=bool)
        if hold_out_platemap is None
        else (meta["platemap"] != hold_out_platemap).to_numpy()
    )
    rows, cols = meta["row"].to_numpy(), meta["col"].to_numpy()
    coef, center = fit_tilt_map(profiles[fit], rows[fit], cols[fit])
    raw = platemap_amplitudes(profiles, meta, coef, center)
    targets = sorted(raw) if hold_out_platemap is None else [hold_out_platemap]
    return {
        "coef": coef,
        "center": center,
        "features": features,
        "amplitudes": {p: shrink_amplitude(raw, p) for p in targets},
        "raw_amplitudes": raw,
    }


def correct_cells(
    cells: pd.DataFrame,
    fit: dict,
    platemap: int,
    features: list[str] | None = None,
) -> pd.DataFrame:
    """Subtract the platemap-scaled tilt from every cell.

    Every cell in a well is shifted by the same feature vector, so differences between
    cells in the same well are unchanged. Controls are corrected like all other wells,
    even though the tilt is estimated from compound wells only. Features that are not
    part of the fit are left unchanged.

    Parameters
    ----------
    cells : pd.DataFrame
        Single-cell profiles of one plate with a "Metadata_Well" column.
    fit : dict
        Output of ``fit_position_correction``.
    platemap : int
        Platemap number of the plate.
    features : list[str] | None, optional
        Features to correct, by default all features in the fit that the cells contain.

    Returns
    -------
    pd.DataFrame
        Copy of the cells with corrected features.
    """
    position = {f: i for i, f in enumerate(fit["features"])}
    use = [f for f in (features or fit["features"]) if f in position and f in cells]
    wells = cells["Metadata_Well"].to_numpy()
    unique_wells = pd.unique(wells)
    rows, cols = well_index(unique_wells)
    offsets = (
        fit["amplitudes"][platemap]
        * tilt(fit["coef"], fit["center"], rows, cols)[:, [position[f] for f in use]]
    )
    offsets = pd.DataFrame(offsets, index=unique_wells, columns=use)
    corrected = cells.copy()
    corrected[use] = cells[use].to_numpy() - offsets.loc[wells].to_numpy()
    return corrected


def save_fit(fit: dict, path: str | pathlib.Path) -> None:
    """Save a fit (tilt map, features and platemap amplitudes) to a compressed file.

    Parameters
    ----------
    fit : dict
        Output of ``fit_position_correction``.
    path : str | pathlib.Path
        Output file, written in NumPy ``.npz`` format.
    """
    amplitudes = fit["amplitudes"]
    np.savez_compressed(
        path,
        coef=fit["coef"].astype("float32"),
        center=fit["center"].astype("float32"),
        features=np.array(fit["features"]),
        amplitude_platemaps=np.array(list(amplitudes), dtype=int),
        amplitude_values=np.array(list(amplitudes.values()), dtype=float),
    )


def load_fit(path: str | pathlib.Path) -> dict:
    """Load a fit written by ``save_fit``.

    Parameters
    ----------
    path : str | pathlib.Path
        File written by ``save_fit``.

    Returns
    -------
    dict
        Same structure as the output of ``fit_position_correction``.
    """
    stored = np.load(path, allow_pickle=False)
    return {
        "coef": stored["coef"].astype(float),
        "center": stored["center"].astype(float),
        "features": stored["features"].tolist(),
        "amplitudes": dict(
            zip(
                stored["amplitude_platemaps"].tolist(),
                stored["amplitude_values"].tolist(),
                strict=True,
            )
        ),
    }
