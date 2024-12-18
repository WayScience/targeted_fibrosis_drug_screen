"""
This module provides functions for loading, writing, and modifying files in this data
repository. It centralizes input/output operations, ensuring efficient and consistent
handling of data files throughout the project.
"""

import pathlib
from collections import defaultdict
import pandas as pd


def load_barcodes(barcode_path: str | pathlib.Path) -> dict:
    """Load barcode data from a given file and organize it into batches.

    This function reads barcode data, groups it by the 'platemap_file' column, and
    organizes the grouped data into batches. Each batch is assigned a unique batch ID
    (e.g., 'batch_1', 'batch_2', etc.).

    Parameters
    ----------
    barcode_path : str | pathlib.Path
        Path to the file containing barcode data. This path can be a string or a pathlib.Path
        object.

    Returns
    -------
    dict
        A dictionary where keys are batch IDs (e.g., "batch_1", "batch_2") and values are
        dictionaries.  Each inner dictionary maps a 'platemap_file' name to the
        corresponding 'plate_barcode' data.
    """

    # type checking
    if not isinstance(barcode_path, (str | pathlib.Path)):
        raise TypeError("'barcode_path' must a str or pathlib.Path")
    if isinstance(barcode_path, str):
        barcode_path = pathlib.Path(barcode_path)

    # load barcode as csv
    barcodes = pd.read_csv(barcode_path)

    # Initialize an empty dictionary to store barcode contents
    barcode_contents = defaultdict(lambda: None)

    # Group barcodes by 'platemap_file' and iterate through each group
    for idx, (platemap_name, df) in enumerate(
        barcodes.groupby("platemap_file"), start=1
    ):
        # Generate a batch ID and plate ID
        batch_id = f"batch_{idx}"

        # Collect plate barcodes for the current platemap
        batch_content = {platemap_name: df["plate_barcode"].values.tolist()}

        # Store the batch content in the main dictionary
        barcode_contents[batch_id] = batch_content

    return barcode_contents
