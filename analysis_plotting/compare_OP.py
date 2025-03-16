#packages
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
from Bio.PDB import *
from analysis_functions import *
from scipy import stats

from sklearn.cluster import DBSCAN
from sklearn.feature_selection import mutual_info_regression
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram


from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.manifold import TSNE


plt.rcParams.update({'font.size': 10})

#IMPORT FILES
mac1_affinity = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/mac1_affinity.csv')
apo_op = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/7kqo_OP.out')
all_files = glob.glob("./OP/*_OP.out") #read in full protein files
li = []
pdb_remove =[]

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    df['PDB'] = filename[5:10]
    li.append(df)

order_all = pd.concat(li, axis=0, ignore_index=True)
crystal_contacts = ['48', '55', '101', '102', '159', '11', '17', '158', '166', '169', '58', '87', '156']
order_all_A = order_all[order_all['chain'] =='A'] 
order_all_A = order_all[~order_all['resi'].isin(crystal_contacts)]

all_files = glob.glob("/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/closeres/*_5.0_closeres.csv") #read in full protein files
li = []
for filename in all_files:
    df = pd.read_csv(filename)
    li.append(df)

close_resi_5 = pd.concat(li, axis=0, ignore_index=True)
close_resi_5 = close_resi_5.rename(columns={
    'PDB': 'atom_name',
    'atom_name': 'PDB'
})

close_resi_5 = close_resi_5[close_resi_5['chain'] !='S']

all_files = glob.glob("/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/closeres/*_10.0_closeres.csv") #read in full protein files
li = []

for filename in all_files:
    df = pd.read_csv(filename)
    li.append(df)

close_resi_10 = pd.concat(li, axis=0, ignore_index=True)

close_resi_10 = close_resi_10.rename(columns={
    'PDB': 'atom_name',
    'atom_name': 'PDB'
})
close_resi_10 = close_resi_10[close_resi_10['chain'] !='S']

close_resi_5['PDB'] = close_resi_5['PDB'].str[:5]
close_resi_10['PDB'] = close_resi_10['PDB'].str[:5]

#determine which atoms should be included
close_resi_5 = close_resi_5[close_resi_5['distance'] < 4]
close_resi_10 = close_resi_10[close_resi_10['distance'] > 10]


#merge apo and all OP
df_merged_order = pd.merge(order_all_A, apo_op, on=['resi', 'chain'], how='left', suffixes=('', '_apo'))
df_merged_order['s2calc_diff'] = df_merged_order['s2calc'] - df_merged_order['s2calc_apo']
order_all_A = df_merged_order
print(order_all_A.head())

order_all_A.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/order_all_A.csv', index=False)
close_resi_5.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/close_resi_5.csv', index=False)
close_resi_10.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/close_resi_10.csv', index=False)

#create close and distal 
close_OP = order_all_A.merge(
    close_resi_5[['resi', 'chain', 'PDB']],    # columns used to match
    on=['resi','chain','PDB'],            # merge keys
    how='inner'                            # only keep matches in both
)

distal_OP = order_all_A.merge(
    close_resi_10[['resi', 'chain', 'PDB']],    # columns used to match
    on=['resi','chain','PDB'],            # merge keys
    how='inner'                            # only keep matches in both
)

close_OP = close_OP.drop_duplicates()
distal_OP = distal_OP.drop_duplicates()

#create mean value for all OP
avg_close_OP = close_OP.groupby('PDB')['s2calc_diff'].mean().reset_index()
avg_distal_OP = distal_OP.groupby('PDB')['s2calc_diff'].mean().reset_index()
avg_OP = order_all_A.groupby('PDB')['s2calc_diff'].mean().reset_index()

avg_OP.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/avg_OP.csv', index=False)
avg_close_OP.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/avg_close_OP.csv', index=False)
avg_distal_OP.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/avg_distal_OP.csv', index=False)

plt.figure()
sns.kdeplot(data=close_OP, x='s2calc_diff', label='Close')
sns.kdeplot(data=distal_OP, x='s2calc_diff', label='Distal')
sns.kdeplot(data=order_all_A, x='s2calc_diff', label='All')
plt.legend()
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/OP_distribution_close_far_all.png')

plt.figure()
sns.kdeplot(data=avg_close_OP, x='s2calc_diff', label='Close')
sns.kdeplot(data=avg_distal_OP, x='s2calc_diff', label='Distal')
sns.kdeplot(data=avg_OP, x='s2calc_diff', label='All')
plt.legend()
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/OP_distribution_close_far.png')

#heatmap of OP
generate_op_heatmap(order_all_A, 'OP_heatmap_all')


# Subset pivot_op to only include the specified PDBs
openers = [ 'x3631',  'x3460', 'x4393',
    'x3458', 'x3439', 'x3466', 'x3476', 'x3623', 'x3422', 'x3924', 'x3406',
    'x3410', 'x3459', 'x3481', 'x3233', 'x3844', 'x4037', 'x3402', 'x3465'
]

series2 = [ 'x3677', 'x3675', 'x3707', 'x3694' ]

pdbs_to_include = series2
order_A_subset = order_all_A[order_all_A['PDB'].isin(openers)]
generate_op_heatmap(order_A_subset, 'OP_heatmap_openers')

plot_op_distribution(order_all_A, 'OP_distribution')
plot_stddev_s2calc(order_all_A, 'OP_sd')

# Calculate the average s2calc for each residue
avg_s2calc = order_all_A.groupby('resi')['s2calc'].mean().reset_index()

# Merge the average s2calc with the standard deviation data
merged_s2calc_sd = pd.merge(avg_s2calc, s2_sd, on='resi', suffixes=('_mean', '_std'))

# Create a scatterplot of the average s2calc vs. the standard deviation
plt.figure(figsize=(10, 8))
sns.scatterplot(x='s2calc_mean', y='s2calc_std', data=merged_s2calc_sd)

# Add a best fit line
sns.regplot(x='s2calc_mean', y='s2calc_std', data=merged_s2calc_sd, scatter=False, color='r')

# Calculate and print the Pearson correlation
pearson_corr, _ = stats.pearsonr(merged_s2calc_sd['s2calc_mean'], merged_s2calc_sd['s2calc_std'])
print(f'Pearson correlation between average s2calc and standard deviation of s2calc: {pearson_corr:.2f}')

plt.xlabel('Average s2calc')
plt.ylabel('Standard Deviation of s2calc')
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/scatter_avg_s2calc_vs_stddev.png')
plt.close()

# Prepare the data for PCA
pca_data = order_all_A.pivot(index='PDB', columns='resi', values='s2calc_diff').fillna(0)

# Perform PCA
pca = PCA(n_components=5)
pca_result = pca.fit_transform(pca_data)
print(pca.explained_variance_ratio_)
# Identify the residues contributing the most to each component of the PCA
pca_components = pd.DataFrame(pca.components_, columns=pca_data.columns, index=[f'PC{i+1}' for i in range(pca.n_components_)])

# Create a DataFrame with the PCA results
pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])
pca_df['PDB'] = pca_data.index

# Prepare the data for t-SNE
tsne_data = order_all_A.pivot(index='PDB', columns='resi', values='s2calc').fillna(0)

# Perform t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_result = tsne.fit_transform(tsne_data)

# Create a DataFrame with the t-SNE results
tsne_df = pd.DataFrame(data=tsne_result, columns=['t-SNE1', 't-SNE2'])
tsne_df['PDB'] = tsne_data.index

# Plot the t-SNE results
plt.figure(figsize=(10, 8))
sns.scatterplot(x='t-SNE1', y='t-SNE2', data=tsne_df)
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/tsne_s2calc.png')
plt.close()

# Plot the PCA results
plt.figure(figsize=(12, 10))
sns.scatterplot(x='PC1', y='PC2', data=pca_df, s=100, palette='viridis', edgecolor='w', alpha=0.7)
plt.xlabel('Principal Component 1', fontsize=20)
plt.ylabel('Principal Component 2', fontsize=20)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/pca_residue_conformational_heterogeneity.png', dpi=300)
plt.close()

plt.figure(figsize=(20, 6))
sns.boxplot(x='resi', y='s2calc', data=order_all_A, flierprops=dict(marker='o', color='gray', markersize=1))
plt.xlabel('Residue')
plt.ylabel('Order Parameter')
plt.xticks([])
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/dist_residue_conformational_heterogeneity.png')

merged_s2calc_sd.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/s2calc_sd.csv')
s2_sd_diff.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/s2calc_diff_avg.csv')

# Read in the DSSP secondary structure data
dssp_secondary_structure = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/dssp_secondary_structure.csv')
# Merge the DSSP data with the s2calc data
merged_data = pd.merge(order_all_A, dssp_secondary_structure, left_on='resi', right_on='Residue')

# Calculate the correlation between different Secondary Structure Descriptions and s2calc
correlation_stats = merged_data.groupby('SecondaryStructureDescription')['s2calc'].describe()
print(correlation_stats)

# Plot the correlation using a violin plot
plt.figure(figsize=(12, 8))
sns.violinplot(x='SecondaryStructureDescription', y='s2calc', data=merged_data)
plt.xlabel('Secondary Structure Description')
plt.ylabel('Average Order Parameter')
plt.xticks(rotation=45)
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/violinplot_secondary_structure_s2calc.png')
plt.tight_layout()
plt.close()

op_vectors = order_all_A.pivot(index="resi", columns="PDB", values="s2calc_diff").fillna(0)
op_vectors = op_vectors.T

correlation_matrix = op_vectors.corr(method="pearson")

# Compute Euclidean distance
euclidean_distances = pd.DataFrame(squareform(pdist(op_vectors, metric='euclidean')),
                                   index=op_vectors.index, columns=op_vectors.index)

# Perform hierarchical clustering
linkage_matrix = linkage(op_vectors, method="ward")

# Assign clusters
from scipy.cluster.hierarchy import fcluster
clusters = fcluster(linkage_matrix, t=1.5, criterion='distance')

# Add cluster information to the DataFrame
clustered_data = op_vectors.copy()
clustered_data['Cluster'] = clusters

# Print the clusters
print(clustered_data['Cluster'])

# Create a color palette for the clusters
import seaborn as sns
import matplotlib.pyplot as plt

# Generate a color palette with a unique color for each cluster
palette = sns.color_palette("hsv", len(set(clusters)))

# Map each cluster to a color
row_colors = clustered_data['Cluster'].map(lambda x: palette[x-1])

# Plot the clustermap with row colors
sns.clustermap(euclidean_distances, cmap="viridis", method="ward", row_colors=row_colors)
plt.title("Euclidean Distances Heatmap with Clusters")
plt.xlabel("PDB Structures")
plt.ylabel("PDB Structures")
plt.tight_layout()
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/euclidean_distances_heatmap_with_clusters_colored.png')
plt.close()

# Print out PDBs in the same cluster
for cluster_id in set(clusters):
    pdbs_in_cluster = clustered_data[clustered_data['Cluster'] == cluster_id].index.tolist()
    print(f"Cluster {cluster_id}: {pdbs_in_cluster}")


# Hierarchical Clustering
linkage_matrix = linkage(op_vectors, method="ward")

# Plot dendrogram
plt.figure(figsize=(8, 5))
dendrogram(linkage_matrix, labels=op_vectors.index, leaf_rotation=90)
plt.title("Clustering of OP Vectors")
plt.xlabel("PDB Structures")
plt.ylabel("Distance")
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/dendrogram_op_vectors.png')

