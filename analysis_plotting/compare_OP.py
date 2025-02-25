#import packages
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
from Bio.PDB import *
from scipy import stats


apo_op = pd.read_csv('7kqo_OP.out')
all_files = glob.glob("./OP/*_OP.out") #read in full protein files
li = []
pdb_remove =[]

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    df['PDB'] = filename[5:10]
    li.append(df)

order_all = pd.concat(li, axis=0, ignore_index=True)

#import files
all_files = glob.glob("*_10.0_closeres.csv") #read in close residues 
li = []

for filename in all_files:
    df = pd.read_csv(filename)
    li.append(df)

close_resi_10 = pd.concat(li, axis=0, ignore_index=True)

close_resi_10 = close_resi_10.rename(columns={
    'PDB': 'atom_name',
    'atom_name': 'PDB'
})

#clean up data
close_resi_10['PDB'] = close_resi_10['PDB'].str[:5]
close_resi_4 = close_resi_10[close_resi_10['distance'] < 4]
distal_resi_10 = close_resi_10[close_resi_10['distance'] > 10]

#find s2calc differences
df_merged_order = pd.merge(order_all_A, apo_op, on=['resi', 'chain'], how='left', suffixes=('', '_apo'))
df_merged_order['s2calc_diff'] = df_merged_order['s2calc'] - df_merged_order['s2calc_apo']

close_OP = order_all_A.merge(
    close_resi_4[['resi', 'chain', 'PDB']],    # columns used to match
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
avg_OP = df_merged_order.groupby('PDB')['s2calc_diff'].mean().reset_index()

#histogram
plt.figure()
sns.kdeplot(data=avg_close_OP, x='s2calc_diff', label='Close')
sns.kdeplot(data=avg_distal_OP, x='s2calc_diff', label='Distal')
sns.kdeplot(data=avg_OP, x='s2calc_diff', label='All')
plt.legend()
plt.savefig('OP_distribution_close_far.png')

#histogram
plt.figure()
sns.kdeplot(data=close_OP, x='s2calc_diff', label='Close')
sns.kdeplot(data=distal_OP, x='s2calc_diff', label='Distal')
sns.kdeplot(data=df_merged_order, x='s2calc_diff', label='All')
plt.legend()
plt.savefig('OP_distribution_close_far_all.png')

#heatmap of OP
s2_sd = order_all_A.groupby(['resi', 'chain'])['s2calc'].std().reset_index()
s2_sd_diff = order_all_A.groupby(['resi', 'chain'])['s2calc_diff'].std().reset_index()
s2_sd_clean = s2_sd.dropna(subset=['s2calc'])
plt.figure()
sns.color_palette("coolwarm", as_cmap=True)
# Create a pivot table for the clustermap
pivot_op = order_all_A.pivot_table(index=['PDB'], columns='resi', values='s2calc_diff').fillna(0)
pivot_op = pivot_op.drop(columns=[150], errors='ignore')
g = sns.clustermap(pivot_op,
               cmap='coolwarm_r',
               figsize=(10, 8), 
               col_cluster=True)
plt.savefig('OP_heatmap_all.png')

# Extract the row and column order from the clustermap
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

# Reorder the pivot_op DataFrame based on the clustermap order
pivot_op_reordered = pivot_op.iloc[row_order, col_order]

# Output the reordered pivot_op DataFrame
pivot_op_reordered.to_csv('pivot_op_diff_reordered.csv')



#______________________PCA__________________________________
pca_data = df_merged_order.pivot(index='PDB', columns='resi', values='s2calc_diff').fillna(0)

# Perform PCA
pca = PCA(n_components=5)
pca_result = pca.fit_transform(pca_data)
print(pca.explained_variance_ratio_)
# Identify the residues contributing the most to each component of the PCA
pca_components = pd.DataFrame(pca.components_, columns=pca_data.columns, index=[f'PC{i+1}' for i in range(pca.n_components_)])


# Plot the PCA results
plt.figure(figsize=(12, 10))
sns.scatterplot(x='PC1', y='PC2', data=pca_df, s=100, palette='viridis', edgecolor='w', alpha=0.7)
plt.xlabel('Principal Component 1', fontsize=20)
plt.ylabel('Principal Component 2', fontsize=20)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig('pca_residue_conformational_heterogeneity.png', dpi=300)
plt.close()



