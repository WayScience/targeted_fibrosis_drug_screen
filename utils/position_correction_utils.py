"""Estimate and remove plate-position ("tilt") effects from single-cell profiles.

Where a well sits on the plate (rows B to G, columns 2 to 11) shifts its profile in a
consistent, feature-specific direction.
Compounds were placed at random, so their own effects average out when we pool many
compounds across platemaps, and the position effect remains.

We estimate this effect once as a shared "tilt map" (an additive row plus column effect
for each feature).
We scale the tilt map with one amplitude per platemap, shrunk toward the average of the
other platemaps, and subtract it from every cell according to the cell's well position.

Use normalized single-cell profiles (step 2 standardizes each plate).
All plates then share the same feature set, so the features line up across platemaps.
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

    The row letter maps to its position in ``ROWS`` (B to G become 0 to 5).
    The column number minus 2 gives the column index (2 to 11 become 0 to 9).
    The column is not range-checked, so a column outside 2 to 11 gives an index outside
    0 to 9.

    Parameters
    ----------
    wells : pd.Series | np.ndarray | list[str]
        Well names such as "B04".

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Row indices and column indices, each with one entry per well.

    Raises
    ------
    ValueError
        If the row letter of a well is not in ``ROWS``.
    """
    wells = np.asarray(wells)
    rows = np.array([ROWS.index(well[0]) for well in wells])
    cols = np.array([int(well[1:]) - 2 for well in wells])
    return rows, cols


def discover_normalized_profiles(base_dir: str | pathlib.Path) -> pd.DataFrame:
    """Find the normalized single-cell profile file of every plate.

    Searches ``<base_dir>/batch_*/platemap_*/single_cell_profiles/`` for files named
    ``<plate>_sc_normalized.parquet``.
    The plate name is the file name without the ``_sc_normalized.parquet`` suffix.
    The platemap number comes from the ``platemap_<N>`` folder name.

    Parameters
    ----------
    base_dir : str | pathlib.Path
        Directory that contains the ``batch_*`` folders.

    Returns
    -------
    pd.DataFrame
        One row per normalized profile file (one per plate), sorted by path, with the
        columns "plate" (str), "platemap" (int), and "path" (pathlib.Path).
        The table is empty, with no columns, when no file matches.
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
    """Summarize one plate as the median profile of each well.

    Each feature is the median over the cells in a well.
    The treatment and cell type come from the first cell of each well, so they must be
    constant within a well.

    Parameters
    ----------
    path : pathlib.Path
        Normalized single-cell parquet file of the plate.
        It needs the columns ``Metadata_Well``, ``Metadata_treatment``, and
        ``Metadata_cell_type``.
        Every column that does not start with ``Metadata_`` is a feature.
    plate : str
        Plate name to store in the ``plate`` column.
    platemap : int
        Platemap number to store in the ``platemap`` column.

    Returns
    -------
    pd.DataFrame
        One row per well, sorted by well name, with the columns ``plate``, ``platemap``,
        ``well``, ``treatment``, ``cell_type``, ``n`` (the number of cells in the well),
        and one float32 column per feature.
        Replicate plates stay separate here; ``compound_table`` combines them.
    """
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
    """Summarize every well of every plate as the median profile of its cells.

    Reads each plate listed in ``plate_table`` with ``_plate_well_medians`` and stacks
    the results.
    Drops every feature with a missing value in any well, because the tilt fit needs
    complete data.

    Parameters
    ----------
    plate_table : pd.DataFrame
        Output of ``discover_normalized_profiles``.
    n_jobs : int, optional
        Number of plates to read in parallel, by default 1.
        The default reads one plate at a time.
        Other values use joblib (``-1`` uses all cores).

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        The well table, with one row per plate and well and the columns ``plate``,
        ``platemap``, ``well``, ``treatment``, ``cell_type``, ``n`` (the number of cells
        in the well), and the retained features, and the list of retained feature names.
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
    """Build one profile per compound well, averaged over replicate plates and centered.

    Each compound sits in one well per platemap, and the replicate plates of a platemap
    share the same layout.
    Averaging a well over the replicate plates therefore gives one profile per compound.

    The function:

    1. Drops control wells (``CONTROL_TREATMENTS``) and wells with fewer than
       ``min_cells`` cells.
    2. Averages the remaining well medians over the replicate plates of each platemap.
    3. Subtracts the mean profile of each platemap (over its compound wells) from its
       profiles, separately for each feature.

    The centering removes differences in the overall level between platemaps.
    The tilt fit then learns the position effect that the platemaps share, and not other
    platemap-to-platemap variation.

    Parameters
    ----------
    well_medians : pd.DataFrame
        Output of ``build_well_medians``.
    features : list[str]
        Feature columns to use.
    min_cells : int, optional
        Minimum number of cells in a well for the well to be used, by default 50.

    Returns
    -------
    tuple[np.ndarray, pd.DataFrame]
        The matrix of compound wells by features, and a table with one row per compound
        well, in the same order as the matrix rows.
        The table has the columns "platemap", "well", "row", and "col", where "row" and
        "col" are the zero-based indices from ``well_index``.
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
    """Build the additive row plus column design matrix of the tilt model.

    Each well gets one indicator for its row.
    It also gets one indicator for its column, except for the first column
    (``cols == 0``), which is the reference and has no column term.
    The row indicators also act as the intercept, so the design has no separate
    intercept column.

    The plate has six rows and ten columns, so the design has ``N_ROW_TERMS`` = 6 row
    terms and ``N_COL_TERMS`` = 9 column terms::

        [row_0, ..., row_5, col_1, ..., col_9]

    A well in row 2 and column 0 has only ``row_2 = 1``.
    A well in row 2 and column 4 has ``row_2 = 1`` and ``col_4 = 1``.

    This represents the position effect as an additive model, with the first column as
    the reference::

        position_effect(row, col) = row_effect(row) + column_effect(col)

    It does not model row-by-column interactions or effects of a single well.

    Parameters
    ----------
    rows : np.ndarray
        Zero-based row index of each well, from 0 to ``N_ROW_TERMS - 1``.
    cols : np.ndarray
        Zero-based column index of each well, from 0 to ``N_COL_TERMS``.
        Index 0 is the reference column.

    Returns
    -------
    np.ndarray
        Indicator matrix with shape ``(n_wells, N_ROW_TERMS + N_COL_TERMS)``.
    """
    design = np.zeros((len(rows), N_ROW_TERMS + N_COL_TERMS))
    design[np.arange(len(rows)), rows] = 1
    later_cols = cols > 0
    design[np.where(later_cols)[0], N_ROW_TERMS + cols[later_cols] - 1] = 1
    return design


def fit_tilt_map(
    profiles: np.ndarray, rows: np.ndarray, cols: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Fit the shared position effect with a robust additive row plus column model.

    The fit has two passes:

    1. Fit the row and column effects of every feature by ordinary least squares.
    2. Winsorize the residuals and refit once on the winsorized data.

    Winsorizing clips each residual to within three robust standard deviations of the
    median residual of its feature, so a few strongly active compounds do not drive the
    estimate.
    The robust standard deviation is the median absolute deviation (MAD) of the
    residuals times 1.4826, which equals the standard deviation for normally distributed
    residuals.

    Parameters
    ----------
    profiles : np.ndarray
        Compound wells by features, centered per platemap (see ``compound_table``).
    rows, cols : np.ndarray
        Zero-based row and column index of each compound well (see ``well_index``).

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The coefficients, with shape ``(N_ROW_TERMS + N_COL_TERMS, n_features)``, and
        ``center``, the mean of the fitted values of each feature over the training
        wells.
        ``tilt`` subtracts ``center``, so the tilt averages to zero over the training
        wells.
    """
    design = design_matrix(rows, cols)
    coef = np.linalg.lstsq(design, profiles, rcond=None)[0]
    MAD_TO_SIGMA = 1.4826
    resid = profiles - design @ coef
    resid_center = np.median(resid, axis=0)
    mad = MAD_TO_SIGMA * np.median(np.abs(resid - resid_center), axis=0) + 1e-9
    winsorized = design @ coef + np.clip(
        resid, resid_center - 3 * mad, resid_center + 3 * mad
    )
    coef = np.linalg.lstsq(design, winsorized, rcond=None)[0]
    return coef, (design @ coef).mean(axis=0)


def tilt(
    coef: np.ndarray, center: np.ndarray, rows: np.ndarray, cols: np.ndarray
) -> np.ndarray:
    """Compute the position effect (tilt) at each well position.

    The tilt of a well is its fitted row plus column effect minus ``center``, so the
    tilt averages to zero over the compound wells that the tilt map was fit on.

    Parameters
    ----------
    coef, center : np.ndarray
        Tilt map from ``fit_tilt_map``.
    rows, cols : np.ndarray
        Zero-based row and column index of each well (see ``well_index``).

    Returns
    -------
    np.ndarray
        Tilt with shape ``(n_wells, n_features)``, in the units of the normalized
        features.
    """
    return design_matrix(rows, cols) @ coef - center


def amplitude(profiles: np.ndarray, tilt_values: np.ndarray) -> float:
    """Find the scale of the tilt that best matches a set of profiles.

    Returns the value ``a`` that minimizes ``sum((profiles - a * tilt_values) ** 2)``
    over all wells and features, which equals
    ``sum(profiles * tilt_values) / sum(tilt_values ** 2)``.
    A value of 1 means the profiles show the tilt at its estimated strength.

    Parameters
    ----------
    profiles : np.ndarray
        Wells by features.
    tilt_values : np.ndarray
        Tilt at the same wells, with the same shape as ``profiles`` (see ``tilt``).

    Returns
    -------
    float
        Amplitude of the tilt in the profiles.
    """
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

    The amplitude of a platemap is the scale of the tilt that best matches its compound
    wells (see ``amplitude``).
    We estimate its standard error by bootstrapping: resample the compound wells of the
    platemap with replacement ``n_boot`` times, recompute the amplitude each time, and
    take the standard deviation of the results.

    Parameters
    ----------
    profiles : np.ndarray
        Matrix of compound wells by features (the first output of ``compound_table``).
    meta : pd.DataFrame
        Table that describes the rows of ``profiles`` (the second output of
        ``compound_table``).
    coef, center : np.ndarray
        Tilt map from ``fit_tilt_map``.
    n_boot : int, optional
        Number of bootstrap resamples, by default 200.
    seed : int, optional
        Seed of the random number generator, by default 0.

    Returns
    -------
    dict[int, tuple[float, float]]
        For each platemap in ``meta``, its amplitude and bootstrap standard error.
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
    """Shrink the amplitude of a platemap toward the average of the other platemaps.

    The amplitude of a single platemap is noisy, so we combine it with the amplitudes of
    the other platemaps.
    Let ``a`` and ``se`` be the amplitudes and standard errors of the other platemaps.

    - ``tau2 = max(0, var(a) - mean(se ** 2))`` is the platemap-to-platemap variance
      that remains after removing the sampling error.
    - ``mean`` is the average of ``a`` weighted by ``1 / (se ** 2 + tau2)``.
    - The result is ``mean + w * (own - mean)``, where ``own`` and ``own_se`` are the
      amplitude and standard error of the target platemap and
      ``w = tau2 / (tau2 + own_se ** 2)``.

    A platemap with a large standard error moves close to the average of the others.
    A platemap with a small standard error keeps most of its own estimate.
    The target platemap does not enter ``mean`` or ``tau2``, and at least two other
    platemaps are needed.

    Parameters
    ----------
    amplitudes : dict[int, tuple[float, float]]
        Output of ``platemap_amplitudes``.
    target : int
        Platemap to shrink; a key of ``amplitudes``.

    Returns
    -------
    float
        Shrunk amplitude of the target platemap.
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
    """Fit the tilt map and the amplitude of every platemap.

    The function:

    1. Builds the compound table (see ``compound_table``).
    2. Fits the tilt map on the compound wells of all platemaps, or of all platemaps
       except ``hold_out_platemap``.
    3. Estimates the amplitude of each platemap against the tilt map (see
       ``platemap_amplitudes``) and shrinks it (see ``shrink_amplitude``).

    Parameters
    ----------
    well_medians : pd.DataFrame
        Output of ``build_well_medians``.
    features : list[str]
        Features to model.
    hold_out_platemap : int | None, optional
        Platemap to leave out of the tilt fit, to evaluate the correction on data that
        the tilt map did not see.
        By default, the fit uses all platemaps.
        The amplitude of the held-out platemap is still estimated from its own compound
        wells, so the evaluation is not fully independent.

    Returns
    -------
    dict
        A dictionary with the keys:

        - "coef" and "center": the tilt map (see ``fit_tilt_map``)
        - "features": the modeled features
        - "amplitudes": the shrunk amplitude of each platemap, or only of the held-out
          platemap when ``hold_out_platemap`` is set
        - "raw_amplitudes": for every platemap, the unshrunk amplitude and its bootstrap
          standard error
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

    The correction of a well is ``amplitude * tilt(row, column)``.
    All cells in a well shift by the same vector, so differences between cells of the
    same well do not change.
    Control wells shift like all other wells, even though the tilt map is estimated from
    compound wells only.
    Features that are not in the fit stay unchanged.

    Parameters
    ----------
    cells : pd.DataFrame
        Normalized single-cell profiles of one plate.
        It needs a "Metadata_Well" column with well names such as "B04" in rows B to G.
    fit : dict
        Output of ``fit_position_correction`` or ``load_fit``.
    platemap : int
        Platemap number of the plate; a key of ``fit["amplitudes"]``.
    features : list[str] | None, optional
        Features to correct.
        By default, every feature of the fit that is a column of ``cells``.
        Listed features that are missing from the fit or from ``cells`` are skipped.

    Returns
    -------
    pd.DataFrame
        Copy of ``cells`` with the corrected features; ``cells`` is not modified.
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
    """Save a fit to a compressed ``.npz`` file.

    The file holds the tilt map (``coef`` and ``center``, stored as float32), the
    feature names, and the shrunk amplitudes in ``fit["amplitudes"]``.
    It does not hold the raw amplitudes or their standard errors.

    Parameters
    ----------
    fit : dict
        Output of ``fit_position_correction``.
    path : str | pathlib.Path
        Output file in NumPy ``.npz`` format.
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
        A dictionary with the keys "coef" and "center" (float64 arrays), "features"
        (list), and "amplitudes" (platemap mapped to its shrunk amplitude).
        These are the keys of the output of ``fit_position_correction``, except for
        "raw_amplitudes", which ``save_fit`` does not store.
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
