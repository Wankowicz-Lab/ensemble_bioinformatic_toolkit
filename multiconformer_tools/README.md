This folder contains scripts that can be used on multiconformer models. These models encode conformational heterogeneity using altlocs. If the PDB or mmCIF have multiple models, please use the ensemble scripts.

## PDB Metrics
1) SASA - by atom 
    This script calculates the **Solvent Accessible Surface Area (SASA)** for all atoms in a **PDB file** using the **Shrake-Rupley algorithm** implemented in **BioPython**. It processes all **atoms except hydrogen**, accounts for **alternate conformations (alt locs)**, and saves the results to a CSV file.
   **How to run**: python compute_sasa.py <pdb_file>



Internally built tools:
1) Hydrogen Bond
2) Protein-ligand interactions
3) Voroni (volume) analysis 



Externally built tools: 
1) fpocket: https://github.com/Discngine/fpocket
    - get pocket size of protein
2) Order Parameters: https://github.com/ExcitedStates/qfit-3.0
    - residue level metric of flexibility

## Comparing PDBs
1) alpha_carbon_rmsd.py
       This script calculated the alpha-carbon RMSD via Pymol. 
       **How to run**: python pymol_rmsd.py structure_1.pdb structure_2.pdb


## Hydrogen Bond Analysis

This will take in a PDB and provide a CSV of all hydrogen bonds weighted by the relative occupancies of the two alternative conformers. 
