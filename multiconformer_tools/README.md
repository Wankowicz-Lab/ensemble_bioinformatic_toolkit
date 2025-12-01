This folder contains scripts to compare information across groups of PDBs. They mainly use outputs gathered from the individual metrics files. 

The majority of scripts also expect a csv identifying PDB or structure IDs to Cluster. PDBs can be clustered in any way you want. 



## Comparing PDBs
1) alpha_carbon_rmsd.py
       This script calculated the alpha-carbon RMSD via Pymol. 
       **How to run**: python pymol_rmsd.py structure_1.pdb structure_2.pdb

2) FGW Pocket Distance Calculator

This script builds ligand-centered binding pockets from a folder of PDB files and computes pairwise **Fused Gromov–Wasserstein (FGW)** distances between them using the POT (`ot`) library.

Builds:
   - `C` : pairwise distance matrix between pocket atoms (geometry)
   - `F` : feature matrix (one-hot element type + ligand/protein flags + residue class)
   - `w` : uniform mass vector over pocket atoms
For each pair of pockets, computes the squared FGW distance:
   - Balancing geometry vs. chemistry via `--alpha`
   - Using entropic regularization (`--epsilon`) via POT’s `fused_gromov_wasserstein`

  
3) cooccurence_network.py
   This script analyzes **residue co-occurrence** across multiple **PDB structures**, considering **occupancy-weighted distances** between residues. It constructs a **network graph** of residue         interactions and saves the results as:
    - A **visual network plot** (`series1_2_residue_network.png`)
    - A **CSV file** containing residue interaction frequencies (`series1_2_residue_network_weights.csv`).
    **How to run**: python compute_residue_network.py <pdb_directory>
