import pandas as pd
from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import AllChem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read in the clustered embeddings CSV file
clustered_embeddings_df = pd.read_csv("/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/HBDScan_clusters_3ormore.csv")
# Rename column 'Unnamed: 0' to 'PDB'
clustered_embeddings_df.rename(columns={'Unnamed: 0': 'PDB'}, inplace=True)
# Remove '_qFit_norm_lig' from all items in the 'PDB' column
clustered_embeddings_df['PDB'] = clustered_embeddings_df['PDB'].str.replace('_qFit_H.pdb', '')

print(clustered_embeddings_df.head())

# Read in the ligand smiles CSV file
ligand_smiles_df = pd.read_csv("ligand_smiles.csv")
# Remove '6_qFit_H.pdb' from pdb_file
ligand_smiles_df['pdb_file'] = ligand_smiles_df['pdb_file'].str.replace('_qFit_H.pdb', '')

# Merge the ligand smiles with the clustered embeddings
merged_df = pd.merge(clustered_embeddings_df, ligand_smiles_df, left_on='PDB', right_on='pdb_file')

# Function to calculate Tanimoto similarity
def calculate_tanimoto(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 and mol2:
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
        return TanimotoSimilarity(fp1, fp2)
    else:
        return None

merged_df.to_csv('merged_ligand_smiles_clusters.csv', index=False)

# Calculate and print the average and range of Tanimoto similarity within each cluster
overall_within_cluster_similarities = []
overall_out_of_cluster_similarities = []

for cluster in merged_df['Cluster'].unique():
    cluster_df = merged_df[merged_df['Cluster'] == cluster]
    cluster_similarities = []

    # Calculate within-cluster similarities
    for i, row1 in cluster_df.iterrows():
        for j, row2 in cluster_df.iterrows():
            if i < j:
                similarity = calculate_tanimoto(row1['smiles'], row2['smiles'])
                if similarity is not None:
                    cluster_similarities.append(similarity)
                    overall_within_cluster_similarities.append(similarity)

    # Calculate similarities with chemicals not in the cluster
    other_clusters_df = merged_df[merged_df['Cluster'] != cluster]
    other_clusters_df.to_csv(f'other_clusters_{cluster}.csv', index=False)
    other_cluster_similarities = []
    for i, row1 in cluster_df.iterrows(): 
        for j, row2 in other_clusters_df.iterrows():
            if row1['smiles'] != row2['smiles']:
                similarity = calculate_tanimoto(row1['smiles'], row2['smiles'])
            else:
                similarity = None
            if similarity is not None and similarity < 0.8:
                other_cluster_similarities.append(similarity)
                overall_out_of_cluster_similarities.append(similarity)

    # Filter out None values before calculating statistics
    filtered_other_cluster_similarities = [sim for sim in other_cluster_similarities if sim is not None]

    num_chemicals = len(cluster_df)
    if cluster_similarities:
        avg_similarity = np.mean(cluster_similarities)
        min_similarity = np.min(cluster_similarities)
        max_similarity = np.max(cluster_similarities)
        print(f"Cluster {cluster}: Number of chemicals = {num_chemicals}, Average Tanimoto similarity = {avg_similarity}, Range = ({min_similarity}, {max_similarity})")
    else:
        print(f"Cluster {cluster}: Number of chemicals = {num_chemicals}, No valid similarities found.")

    num_other_chemicals = len(other_clusters_df)
    if filtered_other_cluster_similarities:
        avg_other_similarity = np.mean(filtered_other_cluster_similarities)
        min_other_similarity = np.min(filtered_other_cluster_similarities)
        max_other_similarity = np.max(filtered_other_cluster_similarities)
        print(f"Cluster {cluster}: Number of chemicals in other clusters = {num_other_chemicals}, Average Tanimoto similarity with other clusters = {avg_other_similarity}, Range = ({min_other_similarity}, {max_other_similarity})")
    else:
        print(f"Cluster {cluster}: Number of chemicals in other clusters = {num_other_chemicals}, No valid similarities with other clusters found.")

# Calculate overall averages
overall_avg_within_cluster = np.median(overall_within_cluster_similarities) if overall_within_cluster_similarities else None
overall_avg_out_of_cluster = np.median(overall_out_of_cluster_similarities) if overall_out_of_cluster_similarities else None

print(f"Overall Average Tanimoto similarity within clusters: {overall_avg_within_cluster}")
print(f"Overall Average Tanimoto similarity out of clusters: {overall_avg_out_of_cluster}")

# Prepare data for plotting
plot_data = []

for cluster in merged_df['Cluster'].unique():
    cluster_df = merged_df[merged_df['Cluster'] == cluster]
    cluster_similarities = []

    # Calculate within-cluster similarities
    for i, row1 in cluster_df.iterrows():
        for j, row2 in cluster_df.iterrows():
            if i < j:
                similarity = calculate_tanimoto(row1['smiles'], row2['smiles'])
                if similarity is not None:
                    cluster_similarities.append(similarity)

    # Calculate similarities with chemicals not in the cluster
    other_clusters_df = merged_df[merged_df['Cluster'] != cluster]
    other_cluster_similarities = []
    for i, row1 in cluster_df.iterrows():
        for j, row2 in other_clusters_df.iterrows():
            if row1['smiles'] != row2['smiles']:
                similarity = calculate_tanimoto(row1['smiles'], row2['smiles'])
            else:
                similarity = None
            if similarity is not None and similarity < 0.8:
                other_cluster_similarities.append(similarity)

    # Add data to plot_data
    plot_data.extend([(cluster, 'Same Cluster', sim) for sim in cluster_similarities])
    plot_data.extend([(cluster, 'Different Cluster', sim) for sim in other_cluster_similarities])

# Convert plot_data to a DataFrame
plot_df = pd.DataFrame(plot_data, columns=['Cluster', 'Type', 'Tanimoto Similarity'])
plot_df.to_csv('tanimoto_by_cluster.csv')
# Plot using seaborn
plt.figure(figsize=(12, 6))
ax = sns.violinplot(
    x='Cluster',
    y='Tanimoto Similarity',
    hue='Type',
    data=plot_df,
    palette={'Same Cluster': '#470659', 'Different Cluster': '#FF8C00'},
    inner=None
)
plt.xlabel('Cluster', fontsize=20, fontweight='bold')
plt.ylabel('Tanimoto Similarity', fontsize=20, fontweight='bold')
plt.xticks(rotation=45, fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(bottom=0)
plt.ylim(top=1)
plt.tight_layout()
# Move the legend outside the plot area (upper left, outside right)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, title='Type', fontsize=16, title_fontsize=18, loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0., frameon=True)
plt.savefig('tanimoto_similarity_distribution.png', bbox_inches='tight')

# Filter the plot_data for the specified clusters
selected_clusters = [0, 3, 5, 6, 9, 10, 13]

# Prepare data for plotting
selected_plot_data = []

# Calculate Tanimoto similarities within selected clusters
for cluster in selected_clusters:
    cluster_df = merged_df[merged_df['Cluster'] == cluster]
    for i, row1 in cluster_df.iterrows():
        for j, row2 in cluster_df.iterrows():
            if i < j:
                similarity = calculate_tanimoto(row1['smiles'], row2['smiles'])
                if similarity is not None:
                    selected_plot_data.append(('Within Selected Clusters', similarity))

# Calculate Tanimoto similarities between ligands in selected clusters and all other ligands
for cluster in selected_clusters:
    cluster_df = merged_df[merged_df['Cluster'] == cluster]
    other_clusters_df = merged_df[~merged_df['Cluster'].isin(selected_clusters)]
    for i, row1 in cluster_df.iterrows():
        for j, row2 in other_clusters_df.iterrows():
            similarity = calculate_tanimoto(row1['smiles'], row2['smiles'])
            if similarity is not None:
                selected_plot_data.append(('Selected vs Others', similarity))

# Convert selected_plot_data to a DataFrame
selected_plot_df = pd.DataFrame(selected_plot_data, columns=['Type', 'Tanimoto Similarity'])

selected_plot_df['Type'] = selected_plot_df['Type'].replace({
    'Within Selected Clusters': 'Ligands with Alt Locs',
    'Selected vs Others': 'Ligands without Alt Locs'
})
print(selected_plot_df.head())

# Plot using seaborn
plt.figure(figsize=(12, 6))
sns.violinplot(x='Type', y='Tanimoto Similarity', data=selected_plot_df, inner=None)
plt.xlabel('', fontsize=20, fontweight='bold')
plt.ylabel('Tanimoto Similarity', fontsize=20, fontweight='bold')
plt.xticks(rotation=45, fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig('selected_clusters_tanimoto_similarity_distribution.png')

# Find closest non-cluster ligands for each selected cluster
print("\nClosest non-cluster ligands for each selected cluster:")
for cluster in selected_clusters:
    cluster_df = merged_df[merged_df['Cluster'] == cluster]
    other_clusters_df = merged_df[~merged_df['Cluster'].isin([cluster])]
    
    closest_matches = []
    for _, cluster_ligand in cluster_df.iterrows():
        max_similarity = -1
        closest_match = None
        
        for _, other_ligand in other_clusters_df.iterrows():
            similarity = calculate_tanimoto(cluster_ligand['smiles'], other_ligand['smiles'])
            if similarity is not None and similarity > max_similarity:
                max_similarity = similarity
                closest_match = (other_ligand['PDB'], similarity)
        
        if closest_match is not None:
            closest_matches.append({
                'cluster_pdb': cluster_ligand['PDB'],
                'closest_pdb': closest_match[0],
                'similarity': closest_match[1]
            })
    
    # Sort matches by similarity and print results
    closest_matches.sort(key=lambda x: x['similarity'], reverse=True)
    print(f"\nCluster {cluster}:")
    for match in closest_matches:
        print(f"PDB {match['cluster_pdb']} -> closest non-cluster match: {match['closest_pdb']} "
              f"(Tanimoto similarity: {match['similarity']:.3f})")

# Calculate average ligand size per cluster
print("\nAverage ligand size per cluster:")
for cluster in merged_df['Cluster'].unique():
    cluster_ligands = merged_df[merged_df['Cluster'] == cluster]
    
    # Calculate size of each ligand by counting atoms in SMILES
    sizes = []
    for smiles in cluster_ligands['smiles']: 
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                size = mol.GetNumAtoms()
                sizes.append(size)
        except:
            continue
            
    if sizes:
        avg_size = sum(sizes) / len(sizes)
        print(f"Cluster {cluster}: {avg_size:.1f} atoms (from {len(sizes)} ligands)")
