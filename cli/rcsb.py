"""
Functions for using the RCSB.
"""

from pathlib import Path
from urllib import request
from urllib.error import HTTPError, URLError
import requests
from time import sleep
import ast

from cli import utils
from cli.paths import DOWNLOADED_PDB_DIR
from cli.constants import (
    PDB_ID_LENGTH,
    RCSB_BIOUNIT_URL,
    RCSB_BUNDLE_URL,
    RCSB_STRUCTURE_URL,
    RCSB_GENERAL_INFO_URL,
    RCSB_ASSEMBLY_INFO_URL,
    RCSB_CONTACT_RETRIES,
)


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


def download_biological_assembly(pdb_id: str, retries: int = RCSB_CONTACT_RETRIES) -> bool:
    """
    Check to see if the biological assembly is available at the RCSB.
    If so, download it.

    If not, check to see if the PDB-like bundle is available.
    If not, just down the standard PDB.
    """
    biounit_url = RCSB_BIOUNIT_URL + f"{pdb_id[1:3].lower()}/{pdb_id.lower()}.pdb1.gz"
    bundle_url = RCSB_BUNDLE_URL + f"{pdb_id[1:3].lower()}/{pdb_id.lower()}/{pdb_id.lower()}-pdb-bundle.tar.gz"
    structure_url = RCSB_STRUCTURE_URL + f"{pdb_id[1:3].lower()}/pdb{pdb_id.lower()}.ent.gz"

    for _ in range(retries):
        try:
            request.urlretrieve(biounit_url, DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb1.gz")
            return True
        except HTTPError:
            try:
                request.urlretrieve(bundle_url, DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb1.tar.gz")
                return True
            except HTTPError:
                try:
                    request.urlretrieve(structure_url, DOWNLOADED_PDB_DIR / f"{pdb_id}.pdb1.gz")
                    return True
                except HTTPError:
                    return False
        except URLError:
            sleep(1)

    return False


def get_pdb_size_metrics(pdb_id: str, retries: int = RCSB_CONTACT_RETRIES) -> dict[str, int]:
    """
    Given a PDB ID from the RCSB, get PDB size metrics.
    """
    for _ in range(retries):
        general_info = ast.literal_eval(requests.get(RCSB_GENERAL_INFO_URL + pdb_id).text)
        assembly_info = ast.literal_eval(requests.get(RCSB_ASSEMBLY_INFO_URL + pdb_id + "/1").text)

        if ("status" in assembly_info) and (assembly_info["status"] == 500):
            continue
        else:
            break

    if (("status" in assembly_info) and (assembly_info["status"] == 500)) or (
        ("status" in general_info) and (general_info["status"] == 500)
    ):
        raise requests.ConnectionError("Couldn't connect to server")

    # if we couldn't get either set of info, just return 0s
    if (("status" in assembly_info) and (assembly_info["status"] == 404)) and (
        ("status" in general_info) and (general_info["status"] == 404)
    ):
        return {
            "atoms": 0,
            "residues": 0,
            "chains": 0,
        }
    # if we got the general info but not the assembly info (e.g. an NMR structure) use the general info
    elif ("status" in assembly_info) and (assembly_info["status"] == 404):
        return {
            "atoms": int(general_info["rcsb_entry_info"].get("deposited_atom_count", 0)),
            "residues": int(general_info["rcsb_entry_info"].get("deposited_modeled_polymer_monomer_count", 0)),
            "chains": int(general_info["rcsb_entry_info"].get("deposited_polymer_entity_instance_count", 0)),
        }

    # otherwise use the assembly info
    return {
        "atoms": int(assembly_info["rcsb_assembly_info"].get("atom_count", 0)),
        "residues": int(assembly_info["rcsb_assembly_info"].get("modeled_polymer_monomer_count", 0)),
        "chains": int(assembly_info["rcsb_assembly_info"].get("polymer_entity_instance_count", 0)),
        # "chains": int(assembly_info["pdbx_struct_assembly"]["oligomeric_count"]),
    }


def download_pdb_file(pdb_id: str) -> Path:
    """
    Download the biological assembly from the RCSB.
    Unzip and save the PDB.
    """
    if not utils.is_pdb_downloaded(pdb_id):
        if download_biological_assembly(pdb_id):
            utils.decompress_pdb(pdb_id)

    return utils.get_downloaded_pdb_path(pdb_id)
