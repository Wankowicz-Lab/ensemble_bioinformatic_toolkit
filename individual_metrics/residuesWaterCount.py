# watersAtResidue.py <inputPDB> <outputCSVPath> [cutoff = 2.0]
# Returns CSV of the number of waters within specified cutoff distance of each residue.
# CSV -> 'chain_id', 'residue_number', 'residue_name', 'water_count'


import os
import sys
import numpy as np
import argparse
import csv
from collections import defaultdict

import biotite.structure as struc
import biotite.structure.io.pdb as pdb


# Reutrns the distance between two 3d positions
def calculate_distance(pos1, pos2):
    return np.linalg.norm(np.array(pos1) - np.array(pos2))


# Consolidates atoms into residue groups that we need
def group_atoms_by_residue(protein_atoms):
    residue_groups = defaultdict(list)
    
    for i in range(len(protein_atoms)):
        chain_id = protein_atoms.chain_id[i]
        res_id = protein_atoms.res_id[i]
        res_name = protein_atoms.res_name[i]
        
        residue_key = (chain_id, res_id, res_name)
        residue_groups[residue_key].append(i)
    
    return residue_groups


def count_waters_near_residue(residue_atom_indices, protein_atoms, water_atoms, max_distance):
    unique_waters = set()
    
    # Iterates all atoms residue-wise to determine the unique waters surrounding them
    for atom_idx in residue_atom_indices:
        atom_pos = protein_atoms.coord[atom_idx]
        
        for water_idx in range(len(water_atoms)):
            water_pos = water_atoms.coord[water_idx]
            distance = calculate_distance(atom_pos, water_pos)
            
            if distance <= max_distance:
                water_res_id = water_atoms.res_id[water_idx]
                unique_waters.add(water_res_id)
    
    return len(unique_waters)


def get_waters_near_protein(pdb_path, output_path, cutoff=2.0):

    # Read PDB file
    pdb_file = pdb.PDBFile.read(pdb_path)
    structure = pdb_file.get_structure()[0]
    
    protein_mask = struc.filter_amino_acids(structure) # model with nothing but protein
    water_mask = struc.filter_solvent(structure) # model with nothing but water
    
    protein_atoms = structure[protein_mask]
    water_atoms = structure[water_mask]
    
    print(f"Found {len(protein_atoms)} protein atoms and {len(water_atoms)} water atoms")
    
    # Group
    residue_groups = group_atoms_by_residue(protein_atoms)
    
    # Count waters 
    results = []
    for residue_key, atom_indices in residue_groups.items():
        chain_id, res_id, res_name = residue_key
        
        water_count = count_waters_near_residue(
            atom_indices, protein_atoms, water_atoms, cutoff
        )
        
        # insert result to output table
        results.append({
            'chain_id': chain_id,
            'residue_number': res_id,
            'residue_name': res_name,
            'water_count': water_count
        })
    
    # sort the result by chain and resi for consistency
    results.sort(key=lambda x: (x['chain_id'], x['residue_number']))
    
    with open(output_path, 'w', newline='') as csvfile:
        fieldnames = ['chain_id', 'residue_number', 'residue_name', 'water_count']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Results written to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Count water molecules near each protein residue from PDB File'
    )
    parser.add_argument('input_pdb', help='Input PDB file path')
    parser.add_argument('output_csv', help='Output CSV file path')
    parser.add_argument(
        '--cutoff', 
        type=float, 
        default=2.0, 
        help='Distance cutoff of potential waters in Angstroms (default = 2.0)'
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_pdb):
        print(f"Error: Input file '{args.input_pdb}' not found")
        sys.exit(1)
    
    get_waters_near_protein(args.input_pdb, args.output_csv, args.cutoff)


if __name__ == '__main__':
    main()
