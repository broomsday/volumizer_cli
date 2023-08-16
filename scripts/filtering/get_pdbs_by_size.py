"""
Helper script to take a list of PDB IDs and return only those whose size or related metrics
fall within given ranges.
"""


from pathlib import Path
import warnings

import typer
from tqdm import tqdm
from biotite.structure.io import load_structure, save_structure, mmtf

from volumizer.pdb import clean_structure
from cli import analysis, pdb, rcsb, utils
from cli.constants import PDB_ID_LENGTH, MAX_RESOLUTION


def main(
    input_list: Path = typer.Argument(..., help="List of PDB IDs to search"),
    output_list: Path = typer.Argument(..., help="Resulting list of IDs that pass search criteria"),
    min_atoms: int = typer.Option(1, help="Minimum number of atoms"),
    max_atoms: int = typer.Option(None, help="Maximum number of atoms"),
    min_residues: int = typer.Option(1, help="Minimum number of residues"),
    max_residues: int = typer.Option(None, help="Maximum number of residues"),
    min_chains: int = typer.Option(1, help="Minimum number of chains"),
    max_chains: int = typer.Option(None, help="Maximum number of chains"),
):
    """
    Subset a PDB list based on some size metrics
    """
    # we'll be saving some data so make sure directories are available
    utils.setup_dirs()

    metrics = {
        "min_atoms": min_atoms,
        "max_atoms": max_atoms,
        "min_residues": min_residues,
        "max_residues": max_residues,
        "min_chains": min_chains,
        "max_chains": max_chains,
    }
    preparation_metrics = {
        "min_atoms": min_atoms,
        "max_atoms": max_atoms * 2 if max_atoms is not None else None,
        "min_residues": min_residues,
        "max_residues": max_residues * 2 if max_residues is not None else None,
        "min_chains": min_chains,
        "max_chains": max_chains * 2 if max_chains is not None else None,
    }

    # get the list of PDBs we need to check
    with open(input_list, mode="r", encoding="utf-8") as in_file:
        # NOTE: split around '.' to ignore any resolution suffixes
        pdb_ids = [line.rstrip().split(".")[0][:PDB_ID_LENGTH] for line in in_file.readlines()]

    # check the PDBs
    satisfied_pdb_ids = []
    for pdb_id in tqdm(pdb_ids):
        if utils.have_pdb_size_metrics_on_file(pdb_id):
            pdb_size_metrics = utils.load_pdb_size_metrics(pdb_id)
        else:
            if not utils.is_pdb_downloaded(pdb_id):
                if not rcsb.download_pdb_file(pdb_id):
                    continue

            if not utils.is_pdb_prepared(pdb_id):
                resolution = rcsb.get_resolution(utils.get_downloaded_pdb_path(pdb_id))
                if resolution is not None and resolution > MAX_RESOLUTION:
                    continue

                # early exit and don't clean if the PDB file is too big
                biological_assembly = rcsb.get_biological_assembly(pdb_id)
                pdb_size_metrics = pdb.get_pdb_size_metrics(biological_assembly)
                if not analysis.pdb_satisfies_metrics(pdb_size_metrics, preparation_metrics):
                    continue

                # size is worth cleaning
                structure = rcsb.get_biological_assembly(pdb_id)
                prepared_structure = clean_structure(structure)
                save_structure(utils.get_prepared_pdb_path(pdb_id), prepared_structure)
            
            pdb_size_metrics = pdb.get_pdb_size_metrics(load_structure(utils.get_prepared_pdb_path(pdb_id)))
            utils.save_pdb_size_metrics(pdb_id, pdb_size_metrics)

        if pdb_size_metrics is None:
            warnings.warn(f"No PDB size metrics for: {pdb_id}")
        elif analysis.pdb_satisfies_metrics(pdb_size_metrics, metrics):
            satisfied_pdb_ids.append(pdb_id)

    with open(output_list, mode="w", encoding="utf-8") as out_file:
        out_file.writelines([f"{pdb}\n" for pdb in satisfied_pdb_ids])

    print(f"Original number of PDBs: {len(pdb_ids)}")
    print(f"Final number of PDBs: {len(satisfied_pdb_ids)}")


if "__main__" in __name__:
    typer.run(main)