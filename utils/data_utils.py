"""
In this module, we have functions that has to do with providing utilts for our
"""
import pandas as pd
from pycytominer.cyto_utils import infer_cp_features


def split_meta_and_features(
    profile: pd.DataFrame,
    compartments: list[str] = ["Nuclei", "Cells", "Cytoplasm"],
    metadata_tag: bool | None = False,
) -> tuple[list[str], list[str]]:
    """Splits metadata and feature column names

    Parameters
    ----------
    profile : pd.DataFrame
        Dataframe containing image-based profile
    compartments : list, optional
        compartments used to generated image-based profiles, by default
        ["Nuclei", "Cells", "Cytoplasm"]
    metadata_tag : Optional[bool], optional
        indicating if the profiles have metadata columns tagged with 'Metadata_'
        , by default False

    Returns
    -------
    tuple[List[str], List[str]]
        Tuple containing metadata and feature column names
    """

    # identify features names
    features_cols = infer_cp_features(profile, compartments=compartments)

    # iteratively search metadata features and retain order if the Metadata tag is not added
    if metadata_tag is False:
        meta_cols = [
            colname
            for colname in profile.columns.tolist()
            if colname not in features_cols
        ]
    else:
        meta_cols = infer_cp_features(profile, metadata=metadata_tag)

    return (meta_cols, features_cols)
