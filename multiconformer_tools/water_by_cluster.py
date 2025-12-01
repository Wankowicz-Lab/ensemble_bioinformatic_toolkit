import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import PDB
import numpy as np
import pandas as pd
from itertools import combinations
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu
from sklearn.metrics import silhouette_score

import warnings
warnings.filterwarnings("ignore")
# Set font sizes
plt.rcParams.update({
    'font.size': 24,
    'axes.labelweight': 'bold',
    'axes.labelsize': 24,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'legend.fontsize': 24
})

boxprops = {'edgecolor': 'k', 'linewidth': 2}
lineprops = {'color': 'k', 'linewidth': 2}

boxplot_kwargs = dict({'boxprops': boxprops, 'medianprops': lineprops,
                       'whiskerprops': lineprops, 'capprops': lineprops,
                       'width': 0.75})

# Read in the data
water_changes_df = pd.read_csv('water_changes_by_residue.csv')
pdb_clusters_df = pd.read_csv('HBDScan_clusters.csv')

order_all = pd.read_csv('order_all_A.csv')


heatmap_data = water_changes_df.drop(columns=['DistanceToLigand', 'ResidueType', 'residueNumber'])

# Plot the overall heatmap for all PDBs
plt.figure(figsize=(12, 10))
g = sns.clustermap(heatmap_data, cmap='magma_r', center=0, 
                   xticklabels=False, yticklabels=False)
plt.tight_layout()
plt.savefig('water_changes_heatmap_all.png')
plt.close()

# Calculate the average number of water changes for each residue across all PDBs
average_water_changes = water_changes_df.drop(columns=['DistanceToLigand', 'ResidueType', 'residueNumber']).mean(axis=1)
water_changes_df['AverageWaterChanges'] = average_water_changes

# Calculate the median change in the number of water molecules for each cluster
water_changes_df = water_changes_df.set_index('residueNumber')
median_changes_by_cluster = pd.DataFrame(columns=['Cluster', 'PDB', 'MedianChange'])

# Iterate over each cluster
for cluster_id in pdb_clusters_df['Cluster'].unique():
    print(cluster_id)
    cluster_pdbs = pdb_clusters_df[pdb_clusters_df['Cluster'] == cluster_id]['PDB']
    cluster_water_changes = water_changes_df[cluster_pdbs]
    rows = []
    for pdb in cluster_pdbs:
        median_water_changes = cluster_water_changes[pdb].mean()
        rows.append({'Cluster': cluster_id, 'PDB': pdb, 'MedianChange': median_water_changes})

    median_changes_by_cluster = pd.concat([median_changes_by_cluster, pd.DataFrame(rows)], ignore_index=True)

mean_changes_by_cluster = median_changes_by_cluster.groupby('Cluster')['MedianChange'].mean().reset_index()
sorted_clusters = mean_changes_by_cluster.sort_values(by='MedianChange')['Cluster']

plt.figure(figsize=(12, 8))
sns.boxplot(data=median_changes_by_cluster, x='Cluster', y='MedianChange', palette='magma_r', order=sorted_clusters, **boxplot_kwargs)
plt.xlabel('Cluster')
plt.ylabel('Median Change in Water Molecules')
plt.xticks(rotation=45)
plt.grid(axis='y')  # Add horizontal grid lines only
plt.tight_layout()
plt.savefig('median_water_changes_distribution_by_cluster.png')
plt.close()


# Iterate over each cluster
for cluster_id in pdb_clusters_df['Cluster'].unique():
    cluster_pdbs = pdb_clusters_df[pdb_clusters_df['Cluster'] == cluster_id]['PDB'].str.replace('_qFit_H.pdb', '', regex=False)
    cluster_water_changes = water_changes_df[cluster_pdbs.tolist()]
    # Create a clustermap for each cluster
    plt.figure(figsize=(12, 10))
    sns.clustermap(
        cluster_water_changes,
        cmap='magma_r',
        center=0,  # Center the colormap at 0
    )
    plt.xlabel('PDB')
    plt.ylabel('Residue Number')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(f'water_changes_heatmap_cluster_{cluster_id}.png')
    plt.close()

    cluster_OP = order_all[order_all['PDB'].isin(cluster_pdbs)]

    # Create a heatmap for water changes
    plt.figure(figsize=(10, 8))
    cluster_water_changes.to_csv(f'water_change_vs_distance_cluster_{cluster_id}.csv')
    #merge_df = cluster_water_changes.merge(cluster_OP, on='resi')
    # Plot water changes for x3430 and x3713 against the distance to the ligand
    if 'DistanceToLigand' in cluster_water_changes.columns:
        plt.figure(figsize=(12, 8))
        for cluster in cluster_pdbs:
            plt.scatter(
                cluster_water_changes.index,
                cluster_water_changes['DistanceToLigand'],
                c=cluster_water_changes[cluster],
                cmap='magma_r',
                vmin=-6,
                vmax=6
            )

        plt.xlabel('Residue', fontsize=16)
        plt.ylabel('Distance to Ligand', fontsize=16)
        plt.colorbar(label='Water Changes')  # Add a color bar to label
        plt.legend()
        plt.tight_layout()
        plt.xticks(rotation=90, fontsize=16)
        plt.yticks(fontsize=16)
        plt.savefig(f'water_change_vs_distance_cluster_{cluster_id}.png')
        plt.close()

use_binding_site_only = False
binding_site_residues = [156, 157, 49, 22, 23, 130, 48, 52, 131, 129, 132, 160, 126, 125, 38, 128, 24, 154, 21]
# Distance metric for PDB-by-residue vectors: 'correlation' (1 - Pearson r) is a good default.
distance_metric = 'eucidean' 

# Make a residue-indexed matrix of shape (n_residues x n_PDBs)
W = water_changes_df.copy()

# Optionally restrict to binding-site residues
if use_binding_site_only:
    W = W.loc[W.index.intersection(binding_site_residues)]

#Transpose matix
X = W.T 
pdb_list = X.index.tolist()

# Map each PDB to its cluster label
cluster_map = dict(zip(
    pdb_clusters_df['PDB'].str.replace('_qFit', '', regex=False),
    pdb_clusters_df['Cluster']
))
labels = np.array([cluster_map[p] for p in pdb_list])

# ----- PAIRWISE DISTANCES -----
# Compute condensed distance vector and squareform matrix
D_condensed = pdist(X.values, metric=distance_metric)
D = squareform(D_condensed)

# Collect within- and between-cluster distances (upper triangle, i<j)
within = []
between = []
for i, j in combinations(range(len(pdb_list)), 2):
    if labels[i] == labels[j]:
        within.append(D[i, j])
    else:
        between.append(D[i, j])

within = np.array(within, dtype=float)
between = np.array(between, dtype=float)

u_stat, mw_p = mannwhitneyu(between, within)  


# ----- PERMUTATION TEST (robust, label-shuffle) -----
rng = np.random.default_rng(42)
obs_diff = between.mean() - within.mean()
n_perm = 5000
perm_diffs = []
for _ in range(n_perm):
    perm_labels = rng.permutation(labels)
    wtmp, btmp = [], []
    for i, j in combinations(range(len(pdb_list)), 2):
        if perm_labels[i] == perm_labels[j]:
            wtmp.append(D[i, j])
        else:
            btmp.append(D[i, j])
    perm_diffs.append(np.mean(btmp) - np.mean(wtmp))
perm_diffs = np.array(perm_diffs)
perm_p = (np.sum(perm_diffs >= obs_diff) + 1) / (n_perm + 1)

print("\n=== Between vs Within (water-pattern distances) ===")
print(f"Metric: {distance_metric}  |  Residues: {'binding site' if use_binding_site_only else 'all'}")
print(f"Pairs: within={len(within)}, between={len(between)}")
print(f"Means: within={within.mean():.3f}, between={between.mean():.3f}, diff={obs_diff:.3f}")
print(f"Mann–Whitney U (between > within): U={u_stat:.0f}, p={mw_p:.3g}")
print(f"Permutation test (label shuffle): p={perm_p:.3g} (n={n_perm})")

# ----- PLOTS -----
plt.figure(figsize=(8,6))
sns.kdeplot(within, label='Same Cluster', linewidth=10, color='#470659')
sns.kdeplot(between, label='Different Cluster', linewidth=10, color='#FF8C00')
plt.xlabel(f'(1 - Pearson)')
plt.ylabel('Density')
plt.legend(fontsize=22, title_fontsize='22', loc='upper right', frameon=True, shadow=True)
plt.tight_layout()
plt.savefig('water_patterns_within_vs_between_kde_binding.png', dpi=300)
plt.close()

#Calculate dispersion
dispersion = []
for c in np.unique(labels):
    idx = np.where(labels == c)[0]
    if len(idx) >= 2:
        Dc = D[np.ix_(idx, idx)]
        tri = Dc[np.triu_indices_from(Dc, k=1)]
        dispersion.append({'Cluster': c, 'MeanWithinDist': tri.mean(), 'N': len(idx)})
dispersion_df = pd.DataFrame(dispersion).sort_values('MeanWithinDist')

plt.figure(figsize=(8,5))
sns.barplot(data=dispersion_df, x='Cluster', y='MeanWithinDist')
plt.ylabel(f'Mean within-cluster distance ({distance_metric})')
plt.tight_layout()
plt.savefig('water_patterns_cluster_dispersion.png', dpi=300)
plt.close()


# ----- CORRELATION BETWEEN ORDER PARAMETERS AND WATER CHANGES -----
print("\nAnalyzing correlation between order parameters and water changes...")

# Prepare data for correlation analysis
correlation_data = []

# Get the water changes data excluding non-numeric columns
water_changes_numeric = water_changes_df.drop(columns=['DistanceToLigand', 'ResidueType', 'AverageWaterChanges'], errors='ignore')

# Iterate through each PDB and calculate correlation with s2calc_diff
for pdb in water_changes_numeric.columns:
    # Get corresponding PDB data from order_all
    order_data = order_all[order_all['PDB'] == pdb]
    
    if not order_data.empty:
        # Remove duplicate resi values by keeping the first occurrence
        order_data_unique = order_data.drop_duplicates(subset='resi', keep='first')
        
        # Merge water changes and order parameter data
        merged_data = pd.DataFrame({
            'water_changes': water_changes_numeric[pdb],
            's2calc_diff': order_data_unique.set_index('resi')['s2calc_diff']
        }).dropna()
        
        if not merged_data.empty:
            # Calculate correlation
            correlation = merged_data['water_changes'].corr(merged_data['s2calc_diff'])
            correlation_data.append({
                'PDB': pdb,
                'correlation': correlation,
                'n_residues': len(merged_data)
            })
            
            # Create scatter plot for each PDB
            plt.figure(figsize=(8, 6))
            sns.scatterplot(data=merged_data, x='s2calc_diff', y='water_changes', alpha=0.6)
            plt.xlabel('Order Parameter Change (s2calc_diff)')
            plt.ylabel('Water Changes')
            plt.title(f'Correlation for {pdb}\nr = {correlation:.3f}, n = {len(merged_data)}')
            plt.tight_layout()
            plt.savefig(f'water_vs_order_correlation_{pdb}.png', dpi=300)
            plt.close()

