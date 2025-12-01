import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns

# Set font and figure params
plt.rcParams.update({
    'font.size': 24,
    'axes.labelsize': 24,
    'axes.labelweight': 'bold',
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'legend.fontsize': 18
})

boxprops = {'edgecolor': 'k', 'linewidth': 2}
lineprops = {'color': 'k', 'linewidth': 2}

boxplot_kwargs = dict({'boxprops': boxprops, 'medianprops': lineprops,
                       'whiskerprops': lineprops, 'capprops': lineprops,
                       'width': 0.75})

# Define the pattern for completed CSV files
completed_csv_pattern = 'PDBs/PDBs_hold/*_packing.csv'

# Read and concatenate all completed CSV files
all_dataframes = []
for file in glob.glob(completed_csv_pattern):
    df = pd.read_csv(file)
    # Extract PDB name from file path - take 5 chars after last slash
    df['PDB_Name'] = file.split('/')[-1][:5]
    all_dataframes.append(df)

# Concatenate all dataframes into a single dataframe
combined_df = pd.concat(all_dataframes, ignore_index=True)

# Calculate the average packing score for each residue by structure
average_packing_by_PDB = combined_df.groupby(['PDB_Name'])['contact_density'].mean().reset_index()

# Plot the distribution of average packing scores
plt.figure(figsize=(12, 8))
plt.hist(average_packing_by_PDB['contact_density'], bins=30, alpha=0.5, edgecolor='black')
plt.xlabel('Packing Score')
plt.ylabel('Frequency')
plt.tight_layout()
plt.savefig('average_packing_scores_distribution.png', dpi=300)

# Read the HBDScan_clusters_3ormore.csv file
hbdscan_clusters_df = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/HBDScan_clusters_3ormore.csv')
hbdscan_clusters_df['PDB_Name'] = hbdscan_clusters_df['PDB'].str.replace('_qFit_H.pdb', '', regex=False)

# Merge the average packing scores with the cluster information
merged_df = pd.merge(average_packing_by_PDB, hbdscan_clusters_df, on='PDB_Name', how='inner')

# Calculate median scores per cluster for ordering
cluster_medians = merged_df.groupby('Cluster')['contact_density'].median().sort_values()
cluster_order = cluster_medians.index.tolist()

# Plot the distribution of average packing scores for each cluster
plt.figure(figsize=(14, 10))
sns.boxplot(data=merged_df, x='Cluster', y='contact_density', palette='magma_r', order=cluster_order)
plt.xlabel('Cluster')
plt.ylabel('Packing Score')
plt.xticks(rotation=45)
plt.yticks()
plt.tight_layout()
plt.savefig('average_packing_scores_by_cluster.png', dpi=300)

# Define residue clusters
resi_clusters = {
    1: [104, 105, 103],
    2: [157, 162, 165], 
    3: [29, 30, 31],
    4: [36, 78, 79, 80, 81, 82, 83, 93, 109, 110, 114, 115, 116, 117, 118, 119, 142, 143],
    5: [37, 39, 45, 59, 60, 62, 94, 95],
    6: [13, 15, 16, 114, 146, 147, 148, 149, 151, 152]
}

# Create a dataframe with resi_cluster assignments for each residue
resi_cluster_mapping = []
for cluster_id, residues in resi_clusters.items():
    for resi in residues:
        resi_cluster_mapping.append({'resi': resi, 'resi_cluster': cluster_id})
resi_cluster_df = pd.DataFrame(resi_cluster_mapping)

# Merge the combined_df with resi_cluster assignments
combined_with_clusters = pd.merge(combined_df, resi_cluster_df, on='resi', how='inner')

# Plot distribution of packing across all PDBs for each resi_cluster
plt.figure(figsize=(14, 8))
# Calculate median for ordering
cluster_medians = combined_with_clusters.groupby('resi_cluster')['contact_density'].median().sort_values()
cluster_order = cluster_medians.index.tolist()

sns.boxplot(data=combined_with_clusters, x='resi_cluster', y='contact_density', 
            palette='magma_r', order=cluster_order)
plt.xlabel('Residue Cluster')
plt.ylabel('Packing Score')
plt.xticks()
plt.yticks()
plt.tight_layout()
plt.savefig('packing_distribution_by_resi_cluster.png', dpi=300)
plt.close()

# Print summary statistics for each residue cluster
print("\nSummary Statistics by Residue Cluster:")
print("="*70)
for cluster_id in cluster_order:
    cluster_data = combined_with_clusters[combined_with_clusters['resi_cluster'] == cluster_id]['contact_density']
    print(f"\nResidue Cluster {cluster_id}:")
    print(f"  Residues: {resi_clusters[cluster_id]}")
    print(f"  N measurements: {len(cluster_data)}")
    print(f"  Mean: {cluster_data.mean():.3f}")
    print(f"  Median: {cluster_data.median():.3f}")
    print(f"  Std: {cluster_data.std():.3f}")
    print(f"  Min: {cluster_data.min():.3f}")
    print(f"  Max: {cluster_data.max():.3f}")

# Calculate average packing scores for each residue cluster x HBDScan cluster
cluster_stats = []
for resi_cluster_id, residues in resi_clusters.items():
    for hbd_cluster in hbdscan_clusters_df['Cluster'].unique():
        # Get PDBs in this HBDScan cluster
        cluster_pdbs = hbdscan_clusters_df[hbdscan_clusters_df['Cluster'] == hbd_cluster]['PDB_Name'].tolist()
        
        # Get data for these residues in these PDBs
        cluster_data = combined_df[
            (combined_df['PDB_Name'].isin(cluster_pdbs)) & 
            (combined_df['resi'].isin(residues))
        ]
        
        if not cluster_data.empty:
            stats = {
                'Resi_Cluster': resi_cluster_id,
                'HBD_Cluster': hbd_cluster,
                'Mean_Packing': cluster_data['contact_density'].mean(),
                'Median_Packing': cluster_data['contact_density'].median(),
                'Std_Packing': cluster_data['contact_density'].std(),
                'Num_Structures': len(cluster_pdbs),
                'Residues': residues
            }
            cluster_stats.append(stats)

cluster_stats_df = pd.DataFrame(cluster_stats)

