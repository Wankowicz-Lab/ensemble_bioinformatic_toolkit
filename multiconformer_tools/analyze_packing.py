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

# Define the pattern for completed CSV files (or point to correct folder)
completed_csv_pattern = '*_packing.csv'

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

# Read cluster file
hbdscan_clusters_df = pd.read_csv('HBDScan_clusters.csv')

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

