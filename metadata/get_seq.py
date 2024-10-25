#!/usr/bin/env python
import Bio.PDB
from Bio.PDB import *
from Bio import PDB
import sys

def get_chain_sequences(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('PDB_structure', pdb_file)

    chain_info = []

    for model in structure:
        for chain in model:
            sequence = ''
            start_residue_number = None

            for residue in chain:
                if PDB.is_aa(residue):
                    if start_residue_number is None:
                        start_residue_number = residue.id[1]
                    sequence += residue.get_resname()

            chain_id = chain.id
            # Append formatted string with chain ID, sequence, and starting residue number
            chain_info.append(f"{chain_id}\t{sequence}\t{start_residue_number}")

    return chain_info

if __name__ == "__main__":
    pdb_file = sys.argv[1]
    chain_sequences = get_chain_sequences(pdb_file)

    # Print results to stdout
    for info in chain_sequences:
        print(info)
