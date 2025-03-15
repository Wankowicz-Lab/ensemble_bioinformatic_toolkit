import torch
import glob
import os
import numpy as np
import networkx as nx
from scipy.spatial import KDTree
import string

ROOT_DIR = "/dors/wankowicz_lab/stephanie/macrodomain/set1/all_norm/"
OUTPUT_DIR = "/dors/wankowicz_lab/stephanie/macrodomain/set1/all_norm/adj"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === One-Hot Encoding Mappings ===
ELEMENTS = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I', 'Se', 'Zn', 'Fe', 'Cu', 'Mg', 'Mn', 'Ca', 'Na', 'K']
ELEMENT_MAP = {elem: i for i, elem in enumerate(ELEMENTS)}

RESIDUE_NAMES = ['LIG', 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
                 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
RESIDUE_MAP = {res: i for i, res in enumerate(RESIDUE_NAMES)}

CHAIN_CHARS = list(string.ascii_uppercase + string.ascii_lowercase + "0123456789")
CHAIN_MAP = {ch: i for i, ch in enumerate(CHAIN_CHARS)}

ATOM_NAMES = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CG1', 'CG2', 'CD', 'CD1', 'CD2', 'CE', 'CE1',
              'CE2', 'CE3', 'CZ', 'CZ2', 'CZ3', 'CH2', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NZ',
              'NH1', 'NH2', 'SD', 'SG', 'OG', 'OG1', 'OD1', 'OD2', 'OE1', 'OE2', 'OH', 'OXT']
ATOM_MAP = {atom: i for i, atom in enumerate(ATOM_NAMES)}

# === One-Hot Encoding Functions ===
def one_hot_encoding(value, mapping):
    vector = torch.zeros(len(mapping))  # Create zero vector
    if value in mapping:
        vector[mapping[value]] = 1  # Set the correct index to 1
    return vector

# Feature extraction functions
def getElementOneHot(element):
    return one_hot_encoding(element, ELEMENT_MAP)

def getResidueOneHot(residue):
    return one_hot_encoding(residue, RESIDUE_MAP)

def getAtomOneHot(atom):
    return one_hot_encoding(atom, ATOM_MAP)

def getChainOneHot(chain):
    return one_hot_encoding(chain, CHAIN_MAP)

def get_files():
    return glob.glob(os.path.join(ROOT_DIR, "*.pdb"))


def parse_pdb_line(pdbline):
    atom = pdbline[0:6].strip()
    atom_num = int(pdbline[6:11].strip())
    atom_name = pdbline[12:16].strip()
    residue_name = pdbline[17:20].strip()
    chain = pdbline[21].strip()
    residue_num = int(pdbline[22:26].strip())
    x = float(pdbline[30:38].strip())
    y = float(pdbline[38:46].strip())
    z = float(pdbline[46:54].strip())
    occupancy = float(pdbline[54:60].strip())
    b_factor = float(pdbline[60:66].strip())
    alt_loc = pdbline[16].strip()
    atom_type = pdbline[77:78].strip()

    return {
        "atom": atom,
        "atom_num": atom_num,
        "atom_name": atom_name,
        "atom_element": atom_type,
        "residue_name": residue_name,
        "chain": chain,
        "residue_num": residue_num,
        "x": x,
        "y": y,
        "z": z,
        "occupancy": occupancy,
        "b_factor": b_factor,
        "alt_loc": alt_loc
    }

def read_pdb(file):
    dataset_name = os.path.basename(file).split("_")[0]
    with open(file, "r") as f:
        lines = f.readlines()
        protein_atoms = []
        ligand_atoms = []
        residues = {}

        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if "HOH" in line:
                    continue
                atom = parse_pdb_line(line)

                # Handle alternate locations by keeping the highest occupancy
                res_key = (atom["chain"], atom["residue_num"], atom["atom_name"])
                if atom["alt_loc"] and atom["alt_loc"] != " ":
                    if res_key in residues:
                        if atom["occupancy"] > residues[res_key]["occupancy"]:
                            residues[res_key] = atom
                    else:
                        residues[res_key] = atom
                else:
                    residues[res_key] = atom

        # Separate ligand and protein atoms
        for atom in residues.values():
            if atom["residue_name"] == "LIG":
                ligand_atoms.append(atom)
            elif atom["residue_name"] != "HOH" or "LIG":  # Exclude water molecules
                protein_atoms.append(atom)

    return dataset_name, protein_atoms, ligand_atoms

def create_graph(protein_coords, ligand_coords, threshold=5.0):
    """
    Creates a NetworkX graph where nodes are protein residues and ligand atoms.
    Edges are formed if within threshold.
    """
    G = nx.Graph()

    # Add protein nodes
    for i, coord in enumerate(protein_coords):
        G.add_node(f"P_{i}", type="protein", pos=tuple(coord))

    # Add ligand nodes
    for j, coord in enumerate(ligand_coords):
        G.add_node(f"L_{j}", type="ligand", pos=tuple(coord))

    # Add edges based on threshold
    for i, p_coord in enumerate(protein_coords):
        for j, l_coord in enumerate(ligand_coords):
            if np.linalg.norm(p_coord - l_coord) < threshold:
                G.add_edge(f"P_{i}", f"L_{j}", weight=1)

    # Create adjacency matrix
    adjacency_matrix = nx.adjacency_matrix(G).todense()

    return G, adjacency_matrix

def save_adj_mx_to_json(adj, dsName):
    adj = adj.tolist()  # Directly convert the NumPy array to a list
    with open(os.path.join(OUTPUT_DIR, f"lig_adj_{dsName}.json"), "w") as f:
        f.write(str(adj))

def sort_by_atom_count(files):
    return sorted(files, key=lambda x: len(read_pdb(x)[2]), reverse=True)

if __name__ == "__main__":
    allLigFiles = get_files()
    allLigFiles = sort_by_atom_count(allLigFiles)

    for idx, testFile in enumerate(allLigFiles):
        if idx >= 225:
            break

        dataset_name, protein_atoms, ligand_atoms = read_pdb(testFile)

        print(f"Processing: {testFile}")
        print(f"  Ligand atom count: {len(ligand_atoms)}")
        print(f"  Protein atom count: {len(protein_atoms)}")

        # Extract coordinates for graph creation
        protein_coords = np.array([[atom["x"], atom["y"], atom["z"]] for atom in protein_atoms])
        ligand_coords = np.array([[atom["x"], atom["y"], atom["z"]] for atom in ligand_atoms])

        # Create graph and adjacency matrix
        G, new_adj_matrix = create_graph(protein_coords, ligand_coords)

        # Create nodes and existing adjacency matrix
        nodes, existing_adj_matrix = create_graph(protein_coords, ligand_coords)

        if nodes is not None:
            # Save node tensors
            torch.save(nodes, os.path.join(OUTPUT_DIR, f"lig_graph_{dataset_name}.pt"))

            # Save existing adjacency matrix
            save_adj_mx_to_json(existing_adj_matrix, dataset_name)

            # Save the new adjacency matrix
            with open(os.path.join(OUTPUT_DIR, f"graph_adj_matrix_{dataset_name}.json"), "w") as f:
                f.write(str(new_adj_matrix.tolist()))
