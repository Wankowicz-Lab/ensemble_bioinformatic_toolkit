import os
import pandas as pd
from Bio import PDB
from Bio.PDB import NeighborSearch, Vector
import argparse
import math
import numpy as np

# D-A Distance: The donor-acceptor distance should generally be between 2.7 to 3.6 Å.
#The D-H-A angle should be close to linear, typically greater than 120°, with an ideal angle being around 180°.

def process_pdb(pdb_file):
    # Parse PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    # Split ligand and protein
    ligand = []
    protein = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0].strip():  # Hetero atoms (ligands)
                    ligand.extend(residue.get_atoms())
                else:  # Protein atoms
                    protein.extend(residue.get_atoms())

    # Check for hydrogens
    has_hydrogens = any(atom.element == 'H' for atom in protein + ligand)

    if not has_hydrogens:
        return
    # Identify hydrogen bonds
    ns = NeighborSearch(protein + ligand)
    hbonds = []
  
    for donor_atom in protein + ligand:
        if donor_atom.element in ['N', 'O']:  # Potential donor 
            donor_altlocs = get_altlocs(donor_atom)
            for donor_altloc in donor_altlocs:
                potential_partners = ns.search(donor_altloc.coord, 3.6)  # 3.6 Angstrom cutoff
                for partner in potential_partners:
                    if partner.element in ['N', 'O', 'F']: #potential acceptor
                        acceptor_altlocs = get_altlocs(partner)
                        for acceptor_altloc in acceptor_altlocs:
                            distance = donor_altloc - acceptor_altloc
                            if 2.7 <= distance <= 3.6:
                                print('distance is good.')
                                # Find hydrogens attached to the donor atom
                                hydrogens = [atom for atom in donor_altloc.get_parent().get_atoms() if atom.element == 'H']

                                # Sort hydrogens by distance to donor atom
                                hydrogens.sort(key=lambda h: np.linalg.norm(donor_altloc.coord - h.coord))

                                # Consider up to two closest hydrogens
                                close_hydrogens = hydrogens[:3]

                                for h_atom in close_hydrogens:
                                    v1 = Vector(donor_altloc.coord)
                                    v2 = Vector(h_atom.coord)
                                    v3 = Vector(acceptor_altloc.coord)
                                    v1_v2 = v1 - v2
                                    v3_v2 = v3 - v2
                                    angle = v1_v2.angle(v3_v2)
                                    angle_degrees = math.degrees(angle)

                                    if angle_degrees > 120:
                                        donor_parent = donor_altloc.get_parent()
                                        acceptor_parent = acceptor_altloc.get_parent()
                                        hbonds.append({
                                            'donor_chain': donor_parent.get_parent().id,
                                            'donor_residue_number': donor_parent.id[1],
                                            'donor_residue_name': donor_parent.resname,
                                            'donor_atom': donor_altloc.name,
                                            'donor_altloc': donor_altloc.altloc,
                                            'donor_occupancy': donor_altloc.occupancy,
                                            'donor_bfactor': donor_altloc.bfactor,
                                            'hydrogen_atom': h_atom.name,
                                            'acceptor_chain': acceptor_parent.get_parent().id,
                                            'acceptor_residue_number': acceptor_parent.id[1],
                                            'acceptor_residue_name': acceptor_parent.resname,
                                            'acceptor_atom': acceptor_altloc.name,
                                            'acceptor_altloc': acceptor_altloc.altloc,
                                            'acceptor_occupancy': acceptor_altloc.occupancy,
                                            'acceptor_bfactor': acceptor_altloc.bfactor,
                                            'distance': distance,
                                            'angle': angle_degrees
                                        })
                                        break  # We've found a valid hydrogen bond, no need to check other hydrogens

    # Create DataFrame and save to CSV
    df = pd.DataFrame(hbonds)
    output_file = os.path.splitext(os.path.basename(pdb_file))[0] + '_hbonds.csv'
    df.to_csv(output_file, index=False)
    print(f"Hydrogen bond information saved to {output_file}")


def get_altlocs(atom):
    """Return a list of atom objects including the original and its alternate locations."""
    if not atom.is_disordered():
        return [atom]
    return [atom.child_dict[altloc] for altloc in atom.disordered_get_id_list()]

parser = argparse.ArgumentParser(description="Process a PDB file to identify hydrogen bonds.")
parser.add_argument("pdb_file", help="Path to the PDB file")
args = parser.parse_args()

process_pdb(args.pdb_file)
