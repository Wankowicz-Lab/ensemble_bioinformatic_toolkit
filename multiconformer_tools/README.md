This folder contains scripts that can be used on multiconformer models. These models encode conformational heterogeneity using altlocs. If the PDB or mmCIF have multiple models, please use the ensemble scripts.


Internally built tools:
1) Hydrogen Bond
2) Protein-ligand interactions 



Externally built tools: 
1) Order Parameters: https://github.com/ExcitedStates/qfit-3.0
    - residue level metric of flexibility

## Comparing PDBs
1) alpha_carbon_rmsd.py
       This script calculated the alpha-carbon RMSD via Pymol. 
       **How to run**: python pymol_rmsd.py structure_1.pdb structure_2.pdb
   
3) cooccurence_network.py
   This script analyzes **residue co-occurrence** across multiple **PDB structures**, considering **occupancy-weighted distances** between residues. It constructs a **network graph** of residue         interactions and saves the results as:
    - A **visual network plot** (`series1_2_residue_network.png`)
    - A **CSV file** containing residue interaction frequencies (`series1_2_residue_network_weights.csv`).
    **How to run**: python compute_residue_network.py <pdb_directory>
