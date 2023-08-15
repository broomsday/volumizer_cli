"""
Functions for using the RCSB.
"""

from pathlib import Path

import biotite.structure as bts
from biotite.structure.io import mmtf
from biotite.database import rcsb as biotite_rcsb

from cli import utils
from cli.constants import PDB_ID_LENGTH


def parse_cluster_file(lines: list[str]) -> set[str]:
    """
    Take the lines from an RCSB cluster file and return a list of all the PDB IDs in the file.
    """
    return {
            f"{line.split()[0].split('_')[0]}\n"
            for line in lines
            if len(line.split()[0].split('_')[0]) == PDB_ID_LENGTH
        }


def build_pdb_set(cluster_file: Path) -> set[str]:
    """
    Get a set of all the PDB IDs we want to download and process.
    """
    with open(cluster_file, mode="r", encoding="utf-8") as fi:
        return parse_cluster_file(fi.readlines())


def download_pdb_file(pdb_id: str) -> bool:
    """
    Download the biological assembly from the RCSB.
    Unzip and save the PDB.
    """
    download_path = utils.get_downloaded_pdb_path(pdb_id)
    if not download_path.is_file():
        try:
            mmtf_file = mmtf.MMTFFile.read(biotite_rcsb.fetch(pdb_id, "mmtf"))
            mmtf_file.write(download_path)
        except (ConnectionError):
            return False

    return True


def get_biological_assembly(pdb_id: str) -> bts.AtomArray:
    """
    Load the biological assembly of a PDB.
    """
    download_path = utils.get_downloaded_pdb_path(pdb_id)
    mmtf_file = mmtf.MMTFFile.read(download_path)

    try:
        biological_assembly = mmtf.get_assembly(mmtf_file, assembly_id="1", model=1)
    except (ValueError, NotImplementedError, IndexError):
        biological_assembly = mmtf.get_structure(mmtf_file, model=1)

    return biological_assembly


def get_pdb_size_metrics(pdb_id: str) -> dict[str, int] | None:
    """
    Given a PDB ID from the RCSB, get PDB size metrics.
    """
    if not download_pdb_file(pdb_id):
        return None

    biological_assembly = get_biological_assembly(pdb_id)

    return {
        "atoms":len(biological_assembly),
        "residues": bts.get_residue_count(biological_assembly),
        "chains": bts.get_chain_count(biological_assembly),
    }
