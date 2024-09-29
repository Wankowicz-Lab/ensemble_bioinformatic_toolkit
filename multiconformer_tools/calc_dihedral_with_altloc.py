
import Bio.PDB
import numpy as np
import pandas as pd
import math
from pathlib import Path
from tqdm import tqdm
from glob import glob

# Function to calculate dihedrals, accounting for alternate locations and occupancy
def calculate_dihedrals_with_altloc(pdb_file, output_csv):
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    chi_angles = []
    backbone_dihedrals = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.resname
                resnum = residue.id[1]
                
                # Get all atoms in the residue, grouped by altloc
                atom_dict = {}
                for atom in residue:
                    altloc = atom.altloc if atom.altloc != ' ' else ''
                    if altloc not in atom_dict:
                        atom_dict[altloc] = []
                    atom_dict[altloc].append(atom)

                # Calculate chi angles for each altloc
                for altloc, atoms in atom_dict.items():
                    # Only calculate if chi angles are defined for this residue type
                    if resname in CHI1:
                        for chi_atoms in CHI1[resname]:
                            chi_angle = calculate_dihedral_angle(residue, chi_atoms, altloc)
                            occupancy = min([residue[atom].get_occupancy() for atom in chi_atoms if residue.has_id(atom)])
                            chi_angles.append((resname, resnum, altloc, chi_angle, occupancy))

    # Save chi angles and occupancies to a CSV file
    chi_df = pd.DataFrame(chi_angles, columns=["Residue", "Residue Number", "AltLoc", "Chi Angle", "Occupancy"])
    chi_df.to_csv(output_csv, index=False)
    
def calculate_dihedral_angle(residue, chi_atoms, altloc):
    # Helper function to calculate the dihedral angle for a given set of atoms
    atoms = [residue[atom].get_vector() for atom in chi_atoms if residue.has_id(atom)]
    if len(atoms) == 4:
        return Bio.PDB.calc_dihedral(*atoms)
    else:
        return np.nan  # Return NaN if we don't have enough atoms for the dihedral calculation
