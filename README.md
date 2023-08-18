# Purpose
Allows processing of many PDB files.

## Current use case
1. Use `scripts/utils/rcsb_cluster_to_ids.py` to get a list of unique IDS
2. Use `scripts/filtering/get_pdbs_by_size.py` to get only PDBs within a size range
3. Use `scripts/filtering/get_pdbs_by_stoichiometry.py` to get only PDBs that are 8+ chain homo-oligomers
4. Use `scripts/filtering/get_pdbs_by_secondary_structure.py` to only PDBs that are 50% helical
5. Run `scripts/volumize.py` on the resulting list of IDs above
6. Run `scripts/filtering/get_pdbs_by_metrics.py` on the same list as above to get the winners