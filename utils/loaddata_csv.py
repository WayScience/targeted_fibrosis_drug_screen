"""Utility functions for building LoadData CSV files from image paths and metadata."""

from pathlib import Path
import re
from typing import Any, Dict, Iterable, List, Optional, Pattern, Set, Tuple, Union

import pandas as pd

# Define specific type hint for regex expressions
PatternLike = Union[str, Pattern[str]]


def _compile_pattern(pattern: PatternLike) -> Pattern[str]:
    """Return a compiled regex pattern from either a string or compiled pattern.

    Parameters
    ----------
    pattern: PatternLike
        A regex expression that splits the image file name into metadata.
            (e.g., well, plate, site)
    """
    if isinstance(pattern, str):
        return re.compile(pattern, flags=re.IGNORECASE)
    return pattern


def collect_image_paths(
    inputs: Iterable[Union[str, Path]], image_extensions: Set[str]
) -> List[Path]:
    """Return image files from a mix of directories and file paths.

    Parameters
    ----------
    inputs : Iterable[Union[str, Path]]
        A list of file paths and/or directory paths to search.
    image_extensions : Set[str]
        Allowed file extensions (e.g. {".tif", ".png"}).

    Returns
    -------
    List[Path]
        Sorted list of unique image file paths found recursively.
    """
    # Instantiate list for image paths
    image_paths = []

    # Loop through the input directories to find images
    for input_path in inputs:
        input_path = Path(input_path).expanduser().resolve()

        # Find valid image from extensions
        if input_path.is_file():
            if input_path.suffix.lower() in image_extensions:
                image_paths.append(input_path)
            continue

        # If the path is directory, then find all paths in directory
        if input_path.is_dir():
            image_paths.extend(
                path
                for path in input_path.rglob("*")
                if path.is_file() and path.suffix.lower() in image_extensions
            )
            continue

        raise FileNotFoundError(f"Image input does not exist: {input_path}")

    return sorted(set(image_paths))


def infer_plate(
    image_path: Path,
    plate_prefix_pattern: PatternLike,
    plate_folder_pattern: PatternLike,
    default_plate: Optional[str] = None,
) -> str:
    """
    Infer plate identifier from filename or parent directories.

    Priority order:
    1. Plate encoded in filename prefix
    2. Plate-like folder name in parent directories
    3. Default plate (if provided)
    4. Fallback to parent directory names

    Parameters
    ----------
    image_path : Path
        Path to image file.
    plate_prefix_pattern : PatternLike
        Regex pattern used to infer the plate from the filename.
    plate_folder_pattern : PatternLike
        Regex pattern used to infer the plate from parent folder names.
    default_plate : Optional[str], default=None
        Fallback plate identifier if inference fails.

    Returns
    -------
    str
        Plate name metadata.
    """
    plate_prefix_pattern = _compile_pattern(plate_prefix_pattern)
    plate_folder_pattern = _compile_pattern(plate_folder_pattern)

    # Find plate name metadata for image path
    filename_match = plate_prefix_pattern.search(image_path.name)
    if filename_match:
        return filename_match.group("plate")

    # Iterate through parent folder to find plate name if necessary
    for parent in image_path.parents:
        folder_match = plate_folder_pattern.search(parent.name)
        if folder_match:
            return folder_match.group(0)

    if default_plate is not None:
        return default_plate

    return (
        image_path.parent.parent.name
        if image_path.parent.parent.name
        else image_path.parent.name
    )


def parse_image_path(
    image_path: Path,
    channel_map: Dict[str, str],
    well_site_channel_pattern: PatternLike,
    plate_prefix_pattern: PatternLike,
    plate_folder_pattern: PatternLike,
    default_plate: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Extract metadata from an image file path.

    Parameters
    ----------
    image_path : Path
        Path to image file.
    channel_map : Dict[str, str]
        Mapping from channel code (e.g. 'd0') to image name.
    well_site_channel_pattern : PatternLike
        Regex pattern used to parse well, site, and channel from the filename.
    plate_prefix_pattern : PatternLike
        Regex pattern used to infer the plate from the filename.
    plate_folder_pattern : PatternLike
        Regex pattern used to infer the plate from parent folder names.
    default_plate : Optional[str], default=None
        Fallback plate identifier if inference fails.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing plate, well, site, and file metadata.
    """
    well_site_channel_pattern = _compile_pattern(well_site_channel_pattern)

    # Find matching metadata from image file names
    match = well_site_channel_pattern.search(image_path.name)
    if match is None:
        raise ValueError(
            "Could not parse well/site/channel from image filename. "
            f"Expected a pattern like B02f01d0 in: {image_path.name}"
        )

    # Set channel for file and path name per image path
    channel = match.group("channel").lower()

    return {
        "Metadata_Plate": infer_plate(
            image_path,
            plate_prefix_pattern=plate_prefix_pattern,
            plate_folder_pattern=plate_folder_pattern,
            default_plate=default_plate,
        ),
        "Metadata_Well": match.group("well").upper(),
        "Metadata_Site": match.group("site").zfill(2),
        "ImageName": channel_map[channel],
        "FileName": image_path.name,
        "PathName": str(image_path.parent),
    }


def build_loaddata_csv(
    image_paths: List[Path],
    channel_map: Dict[str, str],
    loaddata_columns: List[str],
    well_site_channel_pattern: PatternLike,
    plate_prefix_pattern: PatternLike,
    plate_folder_pattern: PatternLike,
    default_plate: Optional[str] = None,
) -> pd.DataFrame:
    """
    Build a LoadData-style dataframe from image paths.

    Parameters
    ----------
    image_paths : List[Path]
        List of image file paths.
    channel_map : Dict[str, str]
        Mapping from channel codes to image names.
    loaddata_columns : List[str]
        Final ordered column list for output dataframe.
    well_site_channel_pattern : PatternLike
        Regex pattern used to parse well, site, and channel from filenames.
    plate_prefix_pattern : PatternLike
        Regex pattern used to infer the plate from filenames.
    plate_folder_pattern : PatternLike
        Regex pattern used to infer the plate from parent folder names.
    default_plate : Optional[str], default=None
        Fallback plate identifier.

    Returns
    -------
    pd.DataFrame
        Wide-format LoadData dataframe with file and path columns.

    Raises
    ------
    ValueError
        If duplicate plate/well/site/channel combinations are found.
    """
    loaddata_columns = [
        column.replace("Image_FileName_", "FileName_").replace(
            "Image_PathName_", "PathName_"
        )
        for column in loaddata_columns
    ]

    # Collect records of metadata and file paths for loaddata csv file
    records = [
        parse_image_path(
            path,
            channel_map=channel_map,
            well_site_channel_pattern=well_site_channel_pattern,
            plate_prefix_pattern=plate_prefix_pattern,
            plate_folder_pattern=plate_folder_pattern,
            default_plate=default_plate,
        )
        for path in image_paths
    ]
    long_df = pd.DataFrame(records)

    # Return empty dataframe
    if long_df.empty:
        return pd.DataFrame(columns=loaddata_columns)

    # Set expected metadata columns
    metadata_columns = [
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_Site",
    ]

    # Drop duplicates if found
    duplicate_mask = long_df.duplicated([*metadata_columns, "ImageName"], keep=False)
    if duplicate_mask.any():
        duplicate_rows = long_df.loc[
            duplicate_mask,
            [*metadata_columns, "ImageName", "FileName"],
        ].sort_values([*metadata_columns, "ImageName"])
        raise ValueError(
            "Found multiple files for the same plate/well/site/channel. "
            "Resolve duplicates before creating the LoadData CSV.\n"
            f"{duplicate_rows.to_string(index=False)}"
        )

    # Create file name dataframe
    filename_df = long_df.pivot(
        index=metadata_columns,
        columns="ImageName",
        values="FileName",
    ).add_prefix("FileName_")

    # Create path name dataframe
    pathname_df = long_df.pivot(
        index=metadata_columns,
        columns="ImageName",
        values="PathName",
    ).add_prefix("PathName_")

    # Create final loaddata csv file
    loaddata_df = filename_df.join(pathname_df).reset_index()

    # Ensure all columns are included but will fill with NaN
    for column in loaddata_columns:
        if column not in loaddata_df.columns:
            loaddata_df[column] = pd.NA

    return loaddata_df.loc[:, loaddata_columns].sort_values(
        ["Metadata_Plate", "Metadata_Well", "Metadata_Site"]
    )


def summarize_missing_channels(loaddata_df: pd.DataFrame, channel_map: Dict[str, str]):
    """
    Identify wells/sites missing one or more imaging channels.

    Parameters
    ----------
    loaddata_df : pd.DataFrame
        LoadData dataframe with FileName_* columns.
    channel_map : Dict[str, str]
        Mapping from channel codes to image names.

    Returns
    -------
    pd.DataFrame
        Table of missing channels per plate/well/site.
    """
    # Set file name columns
    filename_columns = [f"FileName_{image_name}" for image_name in channel_map.values()]
    if not all(column in loaddata_df.columns for column in filename_columns):
        missing_columns = [
            column for column in filename_columns if column not in loaddata_df.columns
        ]
        raise KeyError(
            "Could not find a complete set of LoadData filename columns. "
            f"Missing FileName_* columns: {missing_columns}."
        )

    # Identify any missing files
    missing_mask = loaddata_df[filename_columns].isna().any(axis=1)

    # Check if anything is missing
    if not missing_mask.any():
        return pd.DataFrame(
            columns=[
                "Metadata_Plate",
                "Metadata_Well",
                "Metadata_Site",
                "Missing_Channels",
            ]
        )
    # Subset for any missing files
    missing_df = loaddata_df.loc[
        missing_mask,
        ["Metadata_Plate", "Metadata_Well", "Metadata_Site", *filename_columns],
    ].copy()
    # Find which channels are missing to avoid downstream errors
    missing_df["Missing_Channels"] = missing_df[filename_columns].apply(
        lambda row: ", ".join(
            column.removeprefix("FileName_")
            for column, value in row.items()
            if pd.isna(value)
        ),
        axis=1,
    )

    return missing_df[
        ["Metadata_Plate", "Metadata_Well", "Metadata_Site", "Missing_Channels"]
    ]

def add_illumination_columns_to_loaddata_csvs(
    loaddata_csv_paths: Iterable[Union[str, Path]],
    illum_directory: Union[str, Path],
    channel_map: Dict[str, str],
    output_dir: Union[str, Path],
) -> Dict[Path, pd.DataFrame]:
    """
    Add illumination correction file columns to existing LoadData CSVs and save updated versions.

    Notes
    -----
    - Illumination prefix is assumed to be "Illum"
    - Image prefix is assumed to be "Orig"
    - Output CSVs are always written to output_dir (no in-place overwriting)

    Parameters
    ----------
    loaddata_csv_paths : Iterable[Union[str, Path]]
        Paths to existing LoadData CSV files.
    illum_directory : Union[str, Path]
        Directory containing illumination correction .npy files.
        Expected layout:
        - illum_directory/plate/plate_IllumChannel.npy
        - illum_directory/plate_IllumChannel.npy
    channel_map : Dict[str, str]
        Mapping of channel codes to image names (e.g. {"d0": "OrigActin"}).
        Illumination names are derived by replacing "Orig" → "Illum".
    output_dir : Union[str, Path]
        Directory where updated CSVs will be written.

    Returns
    -------
    Dict[Path, pd.DataFrame]
        Mapping of input CSV paths to updated DataFrames.
    """
    # Resolve path and set up output dir if not already created
    illum_directory = Path(illum_directory).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Raise error if no illum directory exists
    if not illum_directory.is_dir():
        raise FileNotFoundError(f"Illumination directory does not exist: {illum_directory}")

    # Set up channel names with illum functions
    illum_names = [
        f"Illum{image_name.replace('Orig', '')}"
        for image_name in channel_map.values()
    ]

    updated_loaddata_dfs = {}
    missing_files = []

    # Load in CSV file(s) and add illum columns
    for csv_path in loaddata_csv_paths:
        csv_path = Path(csv_path).expanduser().resolve()
        df = pd.read_csv(csv_path)

        # Find Metadata_Plate to use for saving updated LoadData CSV
        if "Metadata_Plate" not in df.columns:
            raise ValueError(f"Missing Metadata_Plate column in {csv_path}")

        # Iterate over unique plate IDs and create a mask for each plate
        for plate in df["Metadata_Plate"].dropna().unique():
            plate_mask = df["Metadata_Plate"] == plate

            # Find each channel illum function per plate
            for illum_name in illum_names:
                # Set the expected illum function npy naming
                illum_file = Path(f"{plate}_{illum_name}.npy")

                # Find the illum functions
                illum_path = illum_directory / str(plate) / illum_file
                if not illum_path.is_file():
                    illum_path = illum_directory / illum_file

                # Set new column names
                fname_col = f"FileName_{illum_name}"
                fpath_col = f"PathName_{illum_name}"

                # Add columns to original loaddata csv
                if illum_path.is_file():
                    df.loc[plate_mask, fname_col] = illum_path.name
                    df.loc[plate_mask, fpath_col] = str(illum_path.parent)
                # Detect any missing illum functions
                else:
                    missing_files.append((csv_path, str(plate), illum_name, illum_path))
                    df.loc[plate_mask, fname_col] = pd.NA
                    df.loc[plate_mask, fpath_col] = pd.NA

        # Append updated loaddata csvs to list
        updated_loaddata_dfs[csv_path] = df

        # Save new loaddata csvs with illum functions as columns
        out_path = output_dir / csv_path.name.replace("loaddata_", "loaddata_with_illum_")
        df.to_csv(out_path, index=False)

    # Report missing illum functions
    if missing_files:
        print("⚠️ Missing illumination files:")
        for csv_path, plate, ch, path in missing_files:
            print(f"{csv_path}: plate={plate}, channel={ch}, expected={path}")

    return updated_loaddata_dfs
