import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import glob
import os
import argparse
from qfit.structure import Structure  # Import qfit.structure

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Process PDB files to compute residue co-occurrence.")
    parser.add_argument("pdb_directory", type=str, help="Directory containing PDB files.")
    return parser.parse_args()

# Load multiple PDBs using qfit.structure
def load_structures(pdb_files):
    structures = [Structure.fromfile(file) for file in pdb_files]
    return structures

# Compute co-occurrence over multiple structures with occupancy consideration
def compute_cooccurrence_multi(structures, distance_threshold=4.0):
    cooccurrence = {}

    for structure in structures:
        residues = list(structure.get_residues())

        for i, res1 in enumerate(residues):
            for j, res2 in enumerate(residues):
                if i >= j:
                    continue  # Avoid duplicate pairs

                # Compute minimum distance between residues considering occupancy
                min_dist = min(
                    (atom1 - atom2) * atom1.occupancy * atom2.occupancy
                    for atom1 in res1 if atom1.element != 'H'
                    for atom2 in res2 if atom2.element != 'H'
                )

                if min_dist <= distance_threshold:
                    key = tuple(sorted([res1.get_resname() + str(res1.get_id()[1]),
                                        res2.get_resname() + str(res2.get_id()[1])]))
                    cooccurrence[key] = cooccurrence.get(key, 0) + 1  # Increment count

    return cooccurrence

# Normalize co-occurrence (frequency over total structures)
def normalize_cooccurrence(cooccurrence, num_structures):
    return {pair: count / num_structures for pair, count in cooccurrence.items()}

# Build network graph
def build_network(cooccurrence):
    G = nx.Graph()
    for (res1, res2), freq in cooccurrence.items():
        G.add_edge(res1, res2, weight=freq)
    return G

# Plot the network
def plot_network(G):
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42)
    edges = G.edges(data=True)
    edge_weights = [data["weight"] for _, _, data in edges]

    nx.draw(G, pos, with_labels=True, node_color="skyblue", edge_color=edge_weights,
            edge_cmap=plt.cm.Blues, node_size=500, font_size=10, width=2)

    plt.savefig("series1_2_residue_network.png", dpi=300)
    # Extract edge weights from the network graph
    edge_data = [(res1, res2, data["weight"]) for res1, res2, data in G.edges(data=True)]

    # Convert to DataFrame
    edge_df = pd.DataFrame(edge_data, columns=["Residue1", "Residue2", "Weight"])

    # Save as CSV
    csv_output_path = os.path.join("series1_2_residue_network_weights.csv")
    edge_df.to_csv(csv_output_path, index=False)

def main():
    args = parse_arguments()
    pdb_files = glob.glob(os.path.join(args.pdb_directory, "*.pdb"))
    structures = load_structures(pdb_files)

    # Compute co-occurrence and normalize
    cooccurrence = compute_cooccurrence_multi(structures)
    normalized_cooccurrence = normalize_cooccurrence(cooccurrence, len(pdb_files))

    # Build and plot network
    G = build_network(normalized_cooccurrence)
    plot_network(G)

if __name__ == "__main__":
    main()
