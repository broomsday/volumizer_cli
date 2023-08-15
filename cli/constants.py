"""
Constants used to avoid magic numbers.
"""

PDB_ID_LENGTH = 4

RCSB_CLUSTER_URL = "https://cdn.rcsb.org/resources/sequence/clusters/bc-90.out"
RCSB_CCD_URL = "https://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif"
RCSB_BIOUNIT_URL = "https://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/divided/"
RCSB_BUNDLE_URL = "https://files.rcsb.org/pub/pdb/compatible/pdb_bundle/"
RCSB_STRUCTURE_URL = "https://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"
RCSB_GENERAL_INFO_URL = "https://data.rcsb.org/rest/v1/core/entry/"
RCSB_ASSEMBLY_INFO_URL = "https://data.rcsb.org/rest/v1/core/assembly/"
RCSB_CONTACT_RETRIES = 10

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
