"""
Command-line entry-point to find pores and cavities in PDBs.
"""

from pathlib import Path

import typer
import multiprocessing
import pandas as pd

from volumizer import volumizer
from volumizer import utils as volumizer_utils
from volumizer.constants import VOXEL_SIZE
from cli import utils as cli_utils
from cli import rcsb
from cli.paths import ANNOTATED_DF_DIR, ANNOTATED_PDB_DIR


def volumize_pdb_id(pdb_id: str) -> None:
    """
    Download the given PDB ID and then volumize it.
    """
    if not cli_utils.have_annotation(pdb_id):
        print(f"Working on: {pdb_id}")

        if not rcsb.download_pdb_file(pdb_id):
            print(f"Skipping {pdb_id}, cannot download")
            return None

        downloaded_pdb_path = cli_utils.get_downloaded_pdb_path(pdb_id)        
        annotated_pdb_path = cli_utils.get_annotated_pdb_path(pdb_id)
        annotated_df_path = annotated_pdb_path.with_suffix(".json")

        volumizer.volumize_pdb_and_save(downloaded_pdb_path, annotated_pdb_path, annotated_df_path)

        print(f"Annotation dataframe saved as: {annotated_df_path}")
        print(f"Annotated PDB saved as: {annotated_pdb_path}")
        print(f"Quick annotation output:")
        print(pd.read_json(annotated_df_path))
    else:
        print(pd.read_json(cli_utils.get_annotated_pdb_path(pdb_id).with_suffix(".json")))        


def volumize_pdb_file(pdb_file: Path) -> None:
    """
    Operate directly on the given PDB file.
    """
    if not cli_utils.have_annotation(pdb_file.stem):
        print(f"Working on: {pdb_file}")

        annotated_pdb_path = (ANNOTATED_PDB_DIR / pdb_file.stem).with_suffix(".pdb")
        annotated_df_path = (ANNOTATED_DF_DIR / pdb_file.stem).with_suffix(".json")

        volumizer.volumize_pdb_and_save(pdb_file, annotated_pdb_path, annotated_df_path)

        print(f"Annotation dataframe saved as: {annotated_df_path}")
        print(f"Annotated PDB saved as: {annotated_pdb_path}")
        print(f"Quick annotation output:")
        print(pd.read_json(annotated_df_path))
    else:
        print(pd.read_json(ANNOTATED_DF_DIR / pdb_file.stem).with_suffix(".json"))


def main(
    volumize_input: str = typer.Argument(
        ..., help="PDB ID, PDB file, file with one PDB ID per line, or folder containing PDB files"
    ),
    resolution: float = typer.Option(VOXEL_SIZE, help="Edge-length of voxels used to discretize the structure."),
    jobs: int = typer.Option(1, help="Number of threads to use."),
):
    """
    Find pores and cavities in the supplied PDB files.
    """
    cli_utils.setup_dirs()
    volumizer_utils.set_resolution(resolution)

    input_type = cli_utils.guess_input_type(volumize_input)

    if input_type == "pdb_id":
        volumize_pdb_id(volumize_input)
    elif input_type == "pdb_file":
        pdb_file = Path(volumize_input)
        volumize_pdb_file(pdb_file)
    elif input_type == "id_file":
        with open(volumize_input, mode="r", encoding="utf-8") as id_file:
            pdb_ids = [line.strip() for line in id_file.readlines()]
        tasks = [[pdb_id] for pdb_id in pdb_ids]
        with multiprocessing.Pool(processes=jobs) as pool:
            pool.starmap(volumize_pdb_id, tasks)
    elif input_type == "pdb_dir":
        pdb_files = Path(volumize_input).glob("*.pdb")
        tasks = [[pdb_file] for pdb_file in pdb_files]
        with multiprocessing.Pool(processes=jobs) as pool:
            pool.starmap(volumize_pdb_file, tasks)
    else:
        raise RuntimeError("File mode not implemented")


if "__main__" in __name__:
    typer.run(main)
