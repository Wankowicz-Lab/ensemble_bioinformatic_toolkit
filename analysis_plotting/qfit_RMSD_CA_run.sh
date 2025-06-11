#!/bin/bash

## Calculate the RMSD between CA atoms of the PDB pairs.
## Need to have the bash script to activate the SBGrid environment to import pymol.
# source /programs/sbgrid.shrc

# Install conda environment for pymol-open-source first. Activate the conda environment in terminal.
# "conda create -n pymol_env python=3.10 pymol-open-source -c schrodinger -c conda-forge"
conda init
# source ~/.bashrc
conda activate pymol_env

# Specify the paths to list of PDB, python script and output files.
FILENAME="/Users/mingbin/Desktop/Wankowicz_lab_Mac/Projects/Catalysis_serine/Ser_bioinfo/list_PDBs/Subtilisin_qfit.txt" # path to pdb list.

PDBPATH="/Users/mingbin/Desktop/Wankowicz_lab_Mac/Projects/Catalysis_serine/Ser_bioinfo/Subtilisin/qFit_061025"

SCRIPTNAME="/Users/mingbin/Desktop/Wankowicz_lab_Mac/Projects/Catalysis_serine/Ser_bioinfo/scripts/analysis_plots/qfit_RMSD_CA.py"

CHAIN="E" # label of catalytic chain. For subtilisin it is chain E. 
## NEED TO RENAME THE CHAIN FOR 1v5i from chain A to become chain E (use `rename_pdb_chain.py` script) ##

OUTNAME="chain_${CHAIN}_CA_RMSD.txt" # output file.

# Run the python script.
python $SCRIPTNAME --file $FILENAME --PDBpath $PDBPATH --chain $CHAIN --output $OUTNAME

echo "Calculations of CA RMSD finished. Results saved to ${OUTNAME}."