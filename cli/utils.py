"""
Functions to interpret CLI commands.
"""


from pathlib import Path
from typing import Optional

import pandas as pd
import json

from cli import paths


def get_downloaded_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return paths.DOWNLOADED_PDB_DIR / f"{pdb_id}.mmtf"


def get_prepared_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return paths.PREPARED_PDB_DIR / f"{pdb_id}.mmtf"


def get_annotated_pdb_path(pdb_id: str) -> Path:
    """
    Return the path to the cleaned PDB file for this PDB ID.
    """
    return paths.ANNOTATED_PDB_DIR / f"{pdb_id}.pdb"


def is_pdb_downloaded(pdb_id: str) -> bool:
    """
    Does the downloaded PDB archive exist.
    """
    downloaded_path = get_downloaded_pdb_path(pdb_id)
    if downloaded_path is None:
        return False

    return get_downloaded_pdb_path(pdb_id).is_file()


def is_pdb_prepared(pdb_id: str) -> bool:
    """
    Has the PDB been fully processed.
    """
    return get_prepared_pdb_path(pdb_id).is_file()


def is_pdb_annotated(pdb_id: str) -> bool:
    """
    Has the PDB been fully processed.
    """
    return get_annotated_pdb_path(pdb_id).is_file()


def setup_dirs():
    """
    Make sure the base data directories exist.
    """
    paths.DOWNLOADED_PDB_DIR.mkdir(parents=True, exist_ok=True)
    paths.PREPARED_PDB_DIR.mkdir(parents=True, exist_ok=True)
    paths.ANNOTATED_PDB_DIR.mkdir(parents=True, exist_ok=True)
    paths.ANNOTATED_DF_DIR.mkdir(parents=True, exist_ok=True)
    paths.PDB_FILTERING_METRIC_DIR.mkdir(parents=True, exist_ok=True)


def save_annotation_dataframe(annotation_df: pd.DataFrame, save_file: Path):
    """
    Save the annotation dataframe.
    """
    annotation_df.to_json(save_file)


def have_annotation(file_stem: str, resolution: float = 3.0) -> bool:
    """
    If we have already completed the annotation of this file, return True.
    False otherwise.
    """
    pdb_path = paths.ANNOTATED_PDB_DIR / f"{file_stem}.{resolution}.pdb"
    df_path = paths.ANNOTATED_DF_DIR / f"{file_stem}.{resolution}.json"
    if pdb_path.is_file() and df_path.is_file():
        return True

    return False


def load_annotation_df(file_stem: str, resolution: float = 3.0) -> pd.DataFrame:
    """
    Return the annotation dataframe associated with this file-stem.
    """
    return pd.read_json(paths.ANNOTATED_DF_DIR / f"{file_stem}.{resolution}.json")


def have_pdb_size_metrics_on_file(pdb_id: str) -> bool:
    """
    Check to see if the PDB size metrics have been recorded previously.
    """
    return (paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_size.json").is_file()


def save_pdb_size_metrics(pdb_id: str, metrics: dict[str, int]) -> None:
    """
    Save the number of atoms, residues, and chains in a PDB assembly.
    """
    with open(paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_size.json", mode="w", encoding="utf-8") as out_file:
        json.dump(metrics, out_file)


def load_pdb_size_metrics(pdb_id: str) -> Optional[dict[str, int]]:
    """
    Load the number of atoms, residues, and chains for a PDB
    """
    pdb_size_metric_path = paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_size.json"
    if not pdb_size_metric_path.is_file():
        return None

    with open(pdb_size_metric_path, mode="r", encoding="utf-8") as in_file:
        return json.load(in_file)


def have_stoichiometry_on_file(pdb_id: str) -> bool:
    """
    Check to see if the stoichiometry of the PDB assembly has been recorded previously.
    """
    return (paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_stoichiometry.json").is_file()


def save_stoichiometry(pdb_id: str, metrics: dict[int, int]) -> None:
    """
    Save the stoichiometry of the PDB assembly.
    """
    with open(paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_stoichiometry.json", mode="w", encoding="utf-8") as out_file:
        json.dump(metrics, out_file)


def load_stoichiometry(pdb_id: str) -> Optional[dict[int, int]]:
    """
    Load the stoichiometry of the PDB assembly.
    """
    pdb_size_metric_path = paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_stoichiometry.json"
    if not pdb_size_metric_path.is_file():
        return None

    with open(pdb_size_metric_path, mode="r", encoding="utf-8") as in_file:
        return json.load(in_file)


def have_secondary_structure_on_file(pdb_id: str) -> bool:
    """
    Check to see if the secondary structure has been recorded previously.
    """
    return (paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_secondary_structure.json").is_file()


def save_secondary_structure(pdb_id: str, metrics: dict[int, int]) -> None:
    """
    Save the secondary structure fractions.
    """
    with open(
        paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_secondary_structure.json", mode="w", encoding="utf-8"
    ) as out_file:
        json.dump(metrics, out_file)


def load_secondary_structure(pdb_id: str) -> Optional[dict[int, int]]:
    """
    Load the secondary structure fractions.
    """
    pdb_size_metric_path = paths.PDB_FILTERING_METRIC_DIR / f"{pdb_id}_secondary_structure.json"
    if not pdb_size_metric_path.is_file():
        return None

    with open(pdb_size_metric_path, mode="r", encoding="utf-8") as in_file:
        return json.load(in_file)


def guess_input_type(input: str) -> Optional[str]:
    """
    Based on the input string, guess if this input is:

    1. a single PDB ID
    2. a single PDB file
    3. a text file containing multiple PDB IDs
    4. a directory containing multiple PDB files.
    """

    if Path(input).is_file():
        if "pdb" in Path(input).suffix:
            return "pdb_file"
        else:
            return "id_file"
    elif Path(input).is_dir():
        return "pdb_dir"
    elif len(input) == 4:
        return "pdb_id"
    else:
        return None


def guess_analysis_input_type(input: str) -> Optional[str]:
    """
    Based on the input string, guess if this input is:

    1. a text containing PDB IDs
    2. a directory containing Dataframes
    """

    if Path(input).is_file():
        return "file"
    elif Path(input).is_dir():
        return "dir"
    else:
        return None