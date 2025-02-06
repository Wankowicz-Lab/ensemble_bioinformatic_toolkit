#!/usr/bin/env python
import sys
from pymol import cmd
import os

# Ensure two arguments are provided (the paths to the two PDB files)
if len(sys.argv) != 3:
    print("Usage: pymol_rmsd.py <structure_1.pdb> <structure_2.pdb>")
    sys.exit(1)

# Load the structures into PyMOL
structure_1 = sys.argv[1]
structure_2 = sys.argv[2]

cmd.load(structure_1, "structure_1")
cmd.load(structure_2, "structure_2")

# Select the alpha carbons (CA) from chain A of both structures
cmd.select("chain_A_1", "structure_1 and chain A and name CA")
cmd.select("chain_A_2", "structure_2 and chain A and name CA")

# Ensure that the selection was made successfully
if cmd.count_atoms("chain_A_1") == 0 or cmd.count_atoms("chain_A_2") == 0:
    print("Error: No alpha carbons found in chain A of one or both structures.")
    sys.exit(1)

# Calculate RMSD between the selected atoms (CA atoms from chain A)
rmsd = cmd.align("chain_A_1", "chain_A_2")[0]

# Prepare the output file name
output_file = f"{os.path.splitext(structure_1)[0]}_RSMD.txt"

# Write the PDB filename and RMSD to the text file
with open(output_file, "w") as f:
    f.write(f"{structure_1}\t{rmsd:.3f}\n")
