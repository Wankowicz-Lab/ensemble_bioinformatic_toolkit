import json
import numpy as np
import pandas as pd
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt
import torch
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.cluster import SpectralClustering

# Load adjacency matrices and node information from JSON and PT files
def load_adjacency_matrices_and_nodes(json_folder):
    matrices = {}
    nodes = {}
    json_files = glob.glob(os.path.join(json_folder, "*.json"))
    pt_files = glob.glob(os.path.join(json_folder, "*.pt"))
    
    for file in json_files:
        with open(file, "r") as f:
            data = json.load(f)
            matrix = np.array(data)  # Assuming JSON contains an adjacency matrix as a nested list
            matrices[os.path.basename(file)] = matrix
    
    for file in pt_files:
        node_data = torch.load(file)
        nodes[os.path.basename(file)] = node_data
    
    # Concatenate node data with matrices
    combined_data = {}
    for name in matrices:
        if name in nodes:
            # Assuming node_data is a 1D array and needs to be concatenated as a new column
            node_array = np.array(nodes[name])
            if node_array.ndim == 1:
                node_array = node_array[:, np.newaxis]  # Convert to 2D for concatenation
            combined_matrix = np.concatenate((matrices[name], node_array), axis=1)
            combined_data[name] = combined_matrix
        else:
            combined_data[name] = matrices[name]  # If no node data, use the matrix as is
    
    return combined_data, nodes

def pad_matrices(matrices):
    max_shape = max((matrix.shape for matrix in matrices.values()), key=lambda x: x[0] * x[1])
    padded_matrices = {}
    for name, matrix in matrices.items():
        padded_matrix = np.zeros(max_shape)
        padded_matrix[:matrix.shape[0], :matrix.shape[1]] = matrix
        padded_matrices[name] = padded_matrix
    return padded_matrices

# Compute pairwise similarity using different metrics
def compute_similarity(matrices, metric="cosine"):
    matrices = pad_matrices(matrices)  # Ensure all matrices are the same shape
    matrix_names = list(matrices.keys())
    flattened_matrices = np.array([matrices[name].flatten() for name in matrix_names])

    if metric == "cosine":
        similarity_matrix = cosine_similarity(flattened_matrices)
    elif metric == "euclidean":
        similarity_matrix = -squareform(pdist(flattened_matrices, metric="euclidean"))  # Negative for similarity
    elif metric == "frobenius":
        similarity_matrix = -squareform(pdist(flattened_matrices, metric=lambda x, y: np.linalg.norm(x - y)))  # Frobenius norm

    return matrix_names, similarity_matrix

# Plot similarity heatmap
def plot_similarity_heatmap(matrix_names, similarity_matrix, metric):
    plt.figure(figsize=(10, 8))
    sns.heatmap(similarity_matrix, xticklabels=matrix_names, yticklabels=matrix_names, cmap="coolwarm", fmt=".2f")
    plt.title(f"Adjacency Matrix Similarity ({metric.capitalize()} Metric)")
    plt.savefig("similarity_heatmap_frobenius.png")

    plt.figure(figsize=(10, 8))
    sns.clustermap(similarity_matrix, xticklabels=matrix_names, yticklabels=matrix_names, cmap="coolwarm", fmt=".2f")
    plt.title(f"Adjacency Matrix Similarity ({metric.capitalize()} Metric)")
    plt.savefig("similarity_clustermap_frobenius.png")

# Output table corresponding to clustermap and identify groups of values in the tightest clusters
def output_clustermap_table(matrix_names, similarity_matrix, metric):
    # Perform hierarchical clustering
    linkage_matrix = linkage(similarity_matrix, method='ward')
    
    # Create a dendrogram to identify clusters
    dendro = dendrogram(linkage_matrix, labels=matrix_names, no_plot=True)
    
    # Extract the cluster assignments
    cluster_assignments = fcluster(linkage_matrix, t=1.2, criterion='distance')
    
    # Create a DataFrame to store the cluster assignments
    cluster_df = pd.DataFrame({
        'Structure': [matrix_names[i] for i in cluster_assignments],
        'Cluster': pd.cut(range(len(cluster_assignments)), bins=len(set(cluster_assignments)), labels=False)
    })
    
    # Save the cluster assignments to a CSV file
    cluster_df.to_csv(f"clustermap_table_{metric}.csv", index=False)
    
    # Print the cluster assignments
    print(f"Cluster assignments for {metric} metric:")
    print(cluster_df)

# Folder containing JSON adjacency matrices and PT node files
json_folder = "/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/adj/"

# Load adjacency matrices and node information
combined_data, nodes = load_adjacency_matrices_and_nodes(json_folder)


# Choose similarity metric: "cosine", "euclidean", "frobenius"
metric = "cosine"
matrix_names, similarity_matrix = compute_similarity(combined_data, metric)

# Output the clustermap table and identify clusters
output_clustermap_table(matrix_names, similarity_matrix, metric)

# Plot similarity
plot_similarity_heatmap(matrix_names, similarity_matrix, metric)

def spectral_clustering_ligands(similarity_matrix, matrix_names, n_clusters=8):
    """Perform spectral clustering on the ligands based on the similarity matrix."""
    clustering = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', random_state=42)
    cluster_labels = clustering.fit_predict(similarity_matrix)
    
    # Create a DataFrame to store the cluster assignments
    cluster_df = pd.DataFrame({
        'Structure': matrix_names,
        'Cluster': cluster_labels
    })
    
    # Save the cluster assignments to a CSV file
    cluster_df.to_csv(f"spectral_clustering_{n_clusters}_clusters.csv", index=False)
    
    # Print the cluster assignments
    print(f"Spectral clustering assignments for {n_clusters} clusters:")
    print(cluster_df)

# Perform spectral clustering on the ligands
spectral_clustering_ligands(similarity_matrix, matrix_names, n_clusters=8)

# Get the top 10 pairs of structures that are similar
def get_top_similar_pairs(matrix_names, similarity_matrix, top_n=10):
    # Get the indices of the upper triangle of the similarity matrix, excluding the diagonal
    upper_triangle_indices = np.triu_indices_from(similarity_matrix, k=1)
    
    # Get the similarity values for the upper triangle
    similarity_values = similarity_matrix[upper_triangle_indices]
    
    # Get the indices of the top N similarity values
    top_indices = np.argsort(similarity_values)[-top_n:]
    
    # Get the corresponding pairs of structure names and their similarity values
    top_pairs = [(matrix_names[i], matrix_names[j], similarity_values[idx]) 
                 for idx, (i, j) in enumerate(zip(upper_triangle_indices[0][top_indices], upper_triangle_indices[1][top_indices]))]
    
    return top_pairs

# Get the top 10 similar pairs
top_similar_pairs = get_top_similar_pairs(matrix_names, similarity_matrix, top_n=20)

# Print the top 10 similar pairs
print("Top 10 similar pairs of structures:")
for pair in top_similar_pairs:
    print(f"Structures: {pair[0]} and {pair[1]}, Similarity: {pair[2]:.4f}")

