This folder contains scripts that can be used on multiconformer models. These models encode conformational heterogeneity using altlocs. If the PDB or mmCIF have multiple models, please use the ensemble scripts.

Internally built tools:
1) Hydrogen Bond
2) Protein-ligand interactions
3) Voroni (volume) analysis 



Externally built tools: 
1) fpocket: https://github.com/Discngine/fpocket
    - get pocket size of protein
2) Order Parameters: https://github.com/ExcitedStates/qfit-3.0
    - residue level metric of flexibility

Comparing PDBs:
1) Alpha-carbon RMSD via Pymol
       python pymol_rmsd.py structure_1.pdb structure_2.pdb


## Hydrogen Bond Analysis

This will take in a PDB and provide a CSV of all hydrogen bonds weighted by the relative occupancies of the two alternative conformers. 
