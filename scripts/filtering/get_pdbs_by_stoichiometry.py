"""
For each PDB file in a list, read in the file, determine the sequence of all chains
and then guess the stoichiometry based on sequence alignment.
"""


from pathlib import Path
import warnings

import typer
from tqdm import tqdm
from biotite.structure.io import load_structure, save_structure

from volumizer.pdb import clean_structure
from cli import utils, pdb, rcsb, analysis
from cli.constants import MAX_RESOLUTION


def main(
    input_list: Path = typer.Argument(..., help=""),
    output_list: Path = typer.Argument(..., help=""),
    min_chain_repeats: int = typer.Option(1, help=""),
    max_chain_repeats: int = typer.Option(None, help=""),
    min_unique_chains: int = typer.Option(1, help=""),
    max_unique_chains: int = typer.Option(None, help=""),
    stoichiometry_factorable: bool = typer.Option(
        False,
        help="If True, only structures where stoichiometry for all chains are factors of one another can pass, e.g. 8-4-2",
    ),
):
    """
    For each PDB file in a list, read in the file, determine the sequence of all chains
    and then guess the stoichiometry based on sequence alignment.
    """
    # we'll be saving some data so make sure directories are available
    utils.setup_dirs()

    metrics = {
        "min_chain_repeats": min_chain_repeats,
        "max_chain_repeats": max_chain_repeats,
        "min_unique_chains": min_unique_chains,
        "max_unique_chains": max_unique_chains,
        "stoichiometry_factorable": stoichiometry_factorable,
    }

    # get the list of PDBs we need to check
    with open(input_list, mode="r", encoding="utf-8") as in_file:
        # NOTE: split around '.' to ignore any resolution suffixes
        pdb_ids = [line.rstrip().split(".")[0] for line in in_file.readlines()]

    # check the PDBS
    satisfied_pdb_ids = []
    for pdb_id in tqdm(pdb_ids):
        if utils.have_stoichiometry_on_file(pdb_id):
            stoichiometry = utils.load_stoichiometry(pdb_id)
        else:
            if not utils.is_pdb_downloaded(pdb_id):
                if not rcsb.download_pdb_file(pdb_id):
                    continue

            if not utils.is_pdb_prepared(pdb_id):
                resolution = rcsb.get_resolution(utils.get_downloaded_pdb_path(pdb_id))
                if resolution is not None and resolution > MAX_RESOLUTION:
                    continue

                structure = rcsb.get_biological_assembly(pdb_id)
                prepared_structure = clean_structure(structure)
                save_structure(utils.get_prepared_pdb_path(pdb_id), prepared_structure)
            
            stoichiometry = pdb.get_stoichiometry(load_structure(utils.get_prepared_pdb_path(pdb_id)))
            utils.save_stoichiometry(pdb_id, stoichiometry)

        if stoichiometry is None:
            warnings.warn(f"No stoichiometry: {pdb_id}")
        elif analysis.pdb_satisfies_stoichiometry(stoichiometry, metrics):
            satisfied_pdb_ids.append(pdb_id)

    with open(output_list, mode="w", encoding="utf-8") as out_file:
        out_file.writelines([f"{pdb}\n" for pdb in satisfied_pdb_ids])

    print(f"Original number of PDBs: {len(pdb_ids)}")
    print(f"Final number of PDBs: {len(satisfied_pdb_ids)}")


if "__main__" in __name__:
    typer.run(main)