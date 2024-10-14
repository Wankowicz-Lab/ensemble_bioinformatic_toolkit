import argparse
from Bio.PDB import PDBParser, NeighborSearch, is_aa
from Bio.PDB.Selection import unfold_entities
import os
import numpy as np
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze protein-ligand interactions from a PDB file.")
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("ligand", help="Ligand of interest (e.g., 'HOH', 'ATP')")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Distance cutoff in angstroms")
    parser.add_argument("--output", default=None, help="Output CSV file")
    return parser.parse_args()

three_to_one_dict = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XAA': 'X', 'UNK': 'X'
}

def three_to_one(residue_name):
    """Convert three-letter amino acid code to one-letter code."""
    return three_to_one_dict.get(residue_name.upper(), 'X')  # 'X' for unknown residues

def get_ligand_atoms(structure, ligand_name):
    ligand_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " " and residue.resname == ligand_name:
                    ligand_atoms.extend(residue.get_atoms())
    return ligand_atoms

def get_residues_within_cutoff(structure, ligand_atoms, cutoff):
    all_atoms = unfold_entities(structure, "A")  # 'A' for atoms
    ns = NeighborSearch(all_atoms)

    nearby_residues = set()
    for ligand_atom in ligand_atoms:
        nearby_atoms = ns.search(ligand_atom.coord, cutoff)
        for atom in nearby_atoms:
            if is_aa(atom.get_parent()):  # Check if atom belongs to an amino acid residue
                nearby_residues.add(atom.get_parent())
    return list(nearby_residues)

def handle_alt_locs(residue):
    """Handles alternate locations, returning the atom with the highest occupancy"""
    atoms = residue.get_unpacked_list()
    alt_locs = {}

    for atom in atoms:
        alt_loc = atom.get_altloc()
        if alt_loc not in alt_locs:
            alt_locs[alt_loc] = []
        alt_locs[alt_loc].append(atom)

    highest_occupancy_atoms = []
    for alt_loc, atom_list in alt_locs.items():
        atom_list.sort(key=lambda a: a.get_occupancy(), reverse=True)
        highest_occupancy_atoms.append(atom_list[0])

    return highest_occupancy_atoms

def is_hydrogen_bond(ligand_atom, residue_atom, distance):
    donor_atoms = ['N', 'O']
    acceptor_atoms = ['N', 'O']
    if ligand_atom.element in donor_atoms and residue_atom.element in acceptor_atoms:
        return 2.2 <= distance <= 3.5
    return False

def is_vdw_interaction(ligand_atom, residue_atom, distance):
    return 3.0 <= distance <= 5.0

def is_halogen_interaction(ligand_atom, residue_atom, distance):
    halogens = ['CL', 'BR', 'I']
    nucleophilic_atoms = ['O', 'N']
    return ligand_atom.element in halogens and residue_atom.element in nucleophilic_atoms and 3.0 <= distance <= 4.0

def is_pi_pi_interaction(ligand_atom, residue_atom, distance):
    def is_benzene_ring(atoms):
        return len(atoms) == 6 and all(atom.element == 'C' for atom in atoms)

    def get_ring_center(atoms):
        return np.mean([atom.coord for atom in atoms], axis=0)

    def get_ring_normal(atoms):
        center = get_ring_center(atoms)
        vectors = [atom.coord - center for atom in atoms]
        return np.cross(vectors[0], vectors[1])

    ligand_residue = ligand_atom.get_parent()
    ligand_ring = [atom for atom in ligand_residue if atom.element == 'C']
    residue_ring = [atom for atom in residue_atom.get_parent() if atom.element == 'C']

    if is_benzene_ring(ligand_ring) and is_benzene_ring(residue_ring):
        ligand_center = get_ring_center(ligand_ring)
        residue_center = get_ring_center(residue_ring)
        center_distance = np.linalg.norm(ligand_center - residue_center)

        if 3.5 <= center_distance <= 6.0:
            ligand_normal = get_ring_normal(ligand_ring)
            residue_normal = get_ring_normal(residue_ring)
            angle = np.degrees(np.arccos(np.dot(ligand_normal, residue_normal) / (np.linalg.norm(ligand_normal) * np.linalg.norm(residue_normal))))

            return 0 <= angle <= 30 or 150 <= angle <= 180

    return False

def classify_interaction(ligand_atom, residue_atom, distance):
    if is_hydrogen_bond(ligand_atom, residue_atom, distance):
        return "Hydrogen Bond"
    elif is_vdw_interaction(ligand_atom, residue_atom, distance):
        return "van der Waals"
    elif is_halogen_interaction(ligand_atom, residue_atom, distance):
        return "Halogen Bond"
    elif is_pi_pi_interaction(ligand_atom, residue_atom, distance):
        return "Pi-Pi Interaction"
    else:
        return "Other"

def analyze_interactions(ligand_atoms, residues):
    interaction_data = []
    for residue in residues:
        residue_atoms = handle_alt_locs(residue)
        for ligand_atom in ligand_atoms:
            for residue_atom in residue_atoms:
                distance = np.linalg.norm(ligand_atom.coord - residue_atom.coord)
                interaction_type = classify_interaction(ligand_atom, residue_atom, distance)
                interaction_data.append({
                    "Residue": residue.get_resname() + str(residue.get_id()[1]),
                    "Residue Atom": residue_atom.get_name(),
                    "Ligand Atom": ligand_atom.get_name(),
                    "Distance (Å)": round(distance, 2),
                    "Interaction": interaction_type
                })
    return interaction_data

def main():
    args = parse_arguments()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", args.pdb_file)

    ligand_atoms = get_ligand_atoms(structure, args.ligand)
    if not ligand_atoms:
        print(f"No ligand named {args.ligand} found in the PDB file.")
        return

    residues = get_residues_within_cutoff(structure, ligand_atoms, args.cutoff)

    interactions = analyze_interactions(ligand_atoms, residues)
    pdb_id = os.path.basename(args.pdb_file).rsplit('.', 1)[0]

    # Save the interaction data to a CSV file
    df = pd.DataFrame(interactions)
    output_filename = os.path.join(os.getcwd(), f"{pdb_id}_interactions.csv")
    df.to_csv(output_filename, index=False)

if __name__ == "__main__":
    main()
~          
