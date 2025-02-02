"""
For each PDB file in a list, read in the file, determine the secondary structure
and return a list of only those matching given secondary structure ranges.
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
    min_helix: float = typer.Option(0.0, help=""),
    max_helix: float = typer.Option(1.0, help=""),
    min_strand: float = typer.Option(0.0, help=""),
    max_strand: float = typer.Option(1.0, help=""),
    min_coil: float = typer.Option(0.0, help=""),
    max_coil: float = typer.Option(1.0, help=""),
):
    """
    For each PDB file in a list, read in the file, determine the secondary structure
    and return a list of only those matching given secondary structure ranges.
    """
    # we'll be saving some data so make sure directories are available
    utils.setup_dirs()

    metrics = {
        "min_helix": min_helix,
        "max_helix": max_helix,
        "min_strand": min_strand,
        "max_strand": max_strand,
        "min_coil": min_coil,
        "max_coil": max_coil,
    }

    # get the list of PDBs we need to check
    with open(input_list, mode="r", encoding="utf-8") as in_file:
        # NOTE: split around '.' to ignore any resolution suffixes
        pdb_ids = [line.rstrip().split(".")[0] for line in in_file.readlines()]

    # check the PDBs
    satisfied_pdb_ids = []
    for pdb_id in tqdm(pdb_ids):
        if utils.have_secondary_structure_on_file(pdb_id):
            secondary_structure = utils.load_secondary_structure(pdb_id)
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

            secondary_structure = pdb.get_secondary_structure(load_structure(utils.get_prepared_pdb_path(pdb_id)))
            utils.save_secondary_structure(pdb_id, secondary_structure)

        if secondary_structure is None:
            warnings.warn(f"No secondary structure: {pdb_id}")
        elif analysis.pdb_satisfies_secondary_structure(secondary_structure, metrics):
            satisfied_pdb_ids.append(pdb_id)

    with open(output_list, mode="w", encoding="utf-8") as out_file:
        out_file.writelines([f"{pdb}\n" for pdb in satisfied_pdb_ids])

    print(f"Original number of PDBs: {len(pdb_ids)}")
    print(f"Final number of PDBs: {len(satisfied_pdb_ids)}")


if "__main__" in __name__:
    typer.run(main)