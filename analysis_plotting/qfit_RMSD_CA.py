#!/usr/bin/env python
import sys
from pymol import cmd
import os
import argparse
import itertools

# Given the list of PDB files, generate pairs of PDB for comparison. Calculate the RMSD of alpha C of the pair.
# Example in the bash script: "python $SCRIPTNAME --file $FILENAME --PDBpath $PDBPATH --chain $CHAIN --output"

def get_args():
    """
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(description="Calculate the alpha C RMSD between two chains of the two PDBs in the given file.")
    parser.add_argument("--file", type=str, required=True, help="Collection of all PDBs.")
    parser.add_argument("--PDBpath", type=str, required=True, help="Directory containing PDB outputs.")
    parser.add_argument("--chain", type=str, required=True, help="Label of catalytic chain.")
    parser.add_argument("--output", type=str, required=True, help="Output file.")
    
    args = parser.parse_args()
    
    return args.file, args.PDBpath, args.chain, args.output

def load_data(filepath):
    """
    Load the CSV file containing the list of PDBs.
    """
    with open(filepath, 'r') as f:
        pdbs = f.readlines()
    pdbs = [a.strip().lower() for a in pdbs] # Generate the list of PDBs in lower case.
    return pdbs

def calc_CA_RMSD(chain, structure_1, structure_2):
    """
    Calculate the RMSD of alpha C atoms of the two given structures.

    Input:
    --------
    chain: string, label of the catalytic chain.
    structure_1: string, PDB ID.
    structure_2: string, PDB ID.

    Output:
    --------
    rmsd value: float
    """
    # Delete the previous structures.
    cmd.reinitialize()
    # Load structures
    cmd.load(structure_1, "structure_1")
    cmd.load(structure_2, "structure_2")

    # Select alpha carbon atoms of the catalytic residues. Could also select by chain.
    cmd.select(f"chain_{chain}_1", f"structure_1 and chain {chain} and name CA")
    cmd.select(f"chain_{chain}_2", f"structure_2 and chain {chain} and name CA")

    # Check if atom selection is valid
    if cmd.count_atoms(f"chain_{chain}_1") == 0:
        print(f"Error: Chain {chain} not found in {structure_1} or contains no CA atoms.")
        sys.exit(1)
    if cmd.count_atoms(f"chain_{chain}_2") == 0:
        print(f"Error: Chain {chain} not found in {structure_2} or contains no CA atoms.")
        sys.exit(1)

    # Compute RMSD
    rmsd = cmd.align(f"chain_{chain}_1", f"chain_{chain}_2")[0]
    return rmsd

# ---------- main block ---------- #
if __name__ == "__main__":
    # Step 1. Load the file of PDBs, label of catalytic chain, and output path.
    filepath, pdb_path, chain, output = get_args()
    pdbs = load_data(filepath)

    # Step 2. Generate unique pairs of PDBs. All need to be in lower case.
    pdb_pairs = list(itertools.combinations(pdbs, 2))
    print("Number of PDB pairs: ", len(pdb_pairs))

    # Step 3. Load the pbd file. For subtilisin, the catalytic chain is chain E.
    strings = [] # output strings to write into output file.
    for pdb_1, pdb_2 in pdb_pairs: # loop through the pairs of pdb files.
        # Point to .pdb file.
        structure_1 = f"{pdb_path}/{pdb_1}_final_qFit.pdb" 
        structure_2 = f"{pdb_path}/{pdb_2}_final_qFit.pdb" 
        # Calculate RMSD value between two each pair.
        try:
            rmsd = calc_CA_RMSD(chain, structure_1, structure_2)
        except Exception as e:
            print(f"{pdb_1} and {pdb_2} pair has error.")
        strings.append(f"{pdb_1}_chain_{chain} {pdb_2}_chain_{chain} {rmsd:.4f}\n") # Write out PDB1, PDB2, rmsd value.

    # Write output file.
    with open(output, 'w') as f:
        f.write("Chain_1 Chain_2 RMSD(Å)\n")
        for line in strings:
            f.write(line) # Write out PDB1, PDB2, rmsd value.
