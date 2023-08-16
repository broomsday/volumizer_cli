"""
Constants used to avoid magic numbers.
"""

PDB_ID_LENGTH = 4
MAX_RESOLUTION = 5.0

RESIDUE_LETTER_CONVERSION = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}
SEQUENCE_IDENTITY_CUTOFF = 0.90

BIOTITE_SSE_CODES = {
    "helix": frozenset(["a"]),
    "strand": frozenset(["b"]),
}
