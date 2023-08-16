"""
Functions for parsing, cleaning, and modifying PDBs.
"""

import biotite.structure as bts
from biotite.sequence import ProteinSequence
from biotite.sequence.align import SubstitutionMatrix, align_optimal

from cli.constants import (
    RESIDUE_LETTER_CONVERSION,
    SEQUENCE_IDENTITY_CUTOFF,
    BIOTITE_SSE_CODES,
)


def get_pdb_size_metrics(structure: bts.AtomArray) -> dict[str, int]:
    """
    Given a PDB ID from the RCSB, get PDB size metrics.
    """
    return {
        "atoms":len(structure),
        "residues": bts.get_residue_count(structure),
        "chains": bts.get_chain_count(structure),
    }


def get_stoichiometry(structure: bts.AtomArray, match_cutoff: float = SEQUENCE_IDENTITY_CUTOFF) -> dict[int, int]:
    """
    Example return = {1: 10, 2: 5} for a heteromultimer with 10 of one chain and 5 of the other
    """
    sequences = []
    for chain in bts.chain_iter(structure):
        sequence = "".join([RESIDUE_LETTER_CONVERSION.get(residue.res_name[0], "X") for residue in bts.residue_iter(chain)])
        sequences.append(sequence)

    # compute pairwise sequence identity using biopython to assign sequence clusters
    cluster_count = 0
    clusters = {cluster_count: [sequences.pop()]}
    while len(sequences) > 0:
        query_sequence = sequences.pop()
        joined_cluster = False
        for cluster_id, cluster_sequences in clusters.items():
            alignment = align_optimal(ProteinSequence(query_sequence), ProteinSequence(cluster_sequences[0]), SubstitutionMatrix.std_protein_matrix(), max_number=1)
            identity = sum([1 for position in alignment[0].trace if position[-1] != -1]) / alignment[0].trace.shape[0]

            # if this is a match, add it to this cluster
            if identity >= match_cutoff:
                clusters[cluster_id].append(query_sequence)
                joined_cluster = True
                break

        # if we don't find any matches, start a new cluster
        if not joined_cluster:
            cluster_count += 1
            clusters[cluster_count] = [query_sequence]

    # compute the stoichiometry
    return {cluster_id: len(cluster_sequences) for cluster_id, cluster_sequences in clusters.items()}


def get_secondary_structure(structure: bts.AtomArray) -> dict[str, float]:
    """
    Compute the fraction of basic secondary structures, helix, strand, loop
    """
    sse = bts.annotate_sse(structure)

    helix_residues = len([res for res in sse if res in BIOTITE_SSE_CODES["helix"]])
    strand_residues = len([res for res in sse if res in BIOTITE_SSE_CODES["strand"]])
    total_sse_residues = len([res for res in sse if res != ""])

    if total_sse_residues == 0:
        return {
            "helix": 0.0,
            "strand": 0.0,
            "coil": 1.0,
        }

    return {
        "helix": helix_residues / total_sse_residues,
        "strand": strand_residues / total_sse_residues,
        "coil": (total_sse_residues - helix_residues - strand_residues) / total_sse_residues,
    }
