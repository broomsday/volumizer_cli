"""
Open an RCSB cluster file and generate a text file with one ID per line,
where each ID is the first ID in the cluster.
"""


from pathlib import Path

import typer

from cli.rcsb import parse_cluster_file


def main(
    cluster_file: Path= typer.Argument(..., help="Sequence ID cluster file downloaded frome RCSB"),
    output_file: Path= typer.Argument(..., help="Output text file listing RCSB IDs matching first entry of each cluster"),
) -> None:
    """
    Open an RCSB cluster file and generate a text file with one ID per line,
    where each ID is the first ID in the cluster.
    """

    with open(cluster_file, mode="r", encoding="utf-8") as file_in:
        pdbs = parse_cluster_file(file_in.readlines())

    with open(output_file, mode="w", encoding="utf-8") as file_out:
        file_out.writelines(pdbs)


if "__main__" in __name__:
    typer.run(main)