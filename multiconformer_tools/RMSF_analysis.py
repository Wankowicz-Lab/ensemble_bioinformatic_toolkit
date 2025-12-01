import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ttest_ind, mannwhitneyu
from sklearn.metrics import silhouette_score
from itertools import combinations

import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({
    'font.size': 24,
    'axes.labelweight': 'bold',
    'axes.labelsize': 24,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'legend.fontsize': 24
})

def calculate_rmsf(residue_coords):
    mean_coords = np.mean(residue_coords, axis=0)
    rmsf = np.sqrt(np.mean(np.sum((residue_coords - mean_coords) ** 2, axis=1)))
    return rmsf


#READ IN RMSF FILE of All Ligand Structures
concatenated_rmsf_df = pd.read_csv('rmsf_output.csv')

# Read in the apo RMSF file
apo_rmsf_df = pd.read_csv('rmsf_output_7KQO.updated_refine_001_ensemble.csv')

#read in cluster csv
clustered_embeddings_df = pd.read_csv("HBDScan_clusters.csv")

apo_rmsf_df.rename(columns={'Resi': 'resi'}, inplace=True)
concatenated_rmsf_df.rename(columns={'Resi': 'resi'}, inplace=True)
merged_rmsf_df = pd.merge(concatenated_rmsf_df, apo_rmsf_df, on='resi', suffixes=('_ensemble', '_apo'))
merged_rmsf_df['delta_RMSF'] = merged_rmsf_df['RMSF_apo'] - merged_rmsf_df['RMSF_ensemble']


def pairwise_similarity(
    df: pd.DataFrame,
    id_col: str = "PDB_ensemble",
    x_col: str = "resi",
    y_col: str = "delta_RMSF",
) -> pd.DataFrame:
    """
    Compute pairwise similarity across all IDs (e.g., PDBs), aligning each pair on the
    intersection of x-values (e.g., residues). Returns an N x N DataFrame.
    """
    series_by_id = {
        pid: g.set_index(x_col)[y_col].sort_index()
        for pid, g in df.groupby(id_col, sort=True)
    }
    ids = list(series_by_id.keys())
    n = len(ids)
    mat = np.full((n, n), np.nan, dtype=float)

    for i in range(n):
        mat[i, i] = 1.0
        si = series_by_id[ids[i]]
        for j in range(i + 1, n):
            sj = series_by_id[ids[j]]
            common = si.index.intersection(sj.index)
            if len(common) == 0:
                continue  # no overlap → leave NaN
            v1 = si.loc[common].to_numpy().reshape(1, -1)
            v2 = sj.loc[common].to_numpy().reshape(1, -1)

            val = cosine_similarity(v1, v2)[0, 0]

            mat[i, j] = mat[j, i] = val

    return pd.DataFrame(mat, index=ids, columns=ids)

# Calculate cosine similarity for PDBs in the same cluster and different clusters
same_cluster_similarities = []
different_cluster_similarities = []

# Iterate over each unique cluster
for cluster_id in clustered_embeddings_df['Cluster'].unique():
    # Get PDBs in the current cluster
    cluster_pdbs = clustered_embeddings_df[clustered_embeddings_df['Cluster'] == cluster_id]['PDB']
    
    # Calculate pairwise similarities within the same cluster
    for i, pdb1 in enumerate(cluster_pdbs):
        for pdb2 in cluster_pdbs[i+1:]:
            if pdb1 in cos_sim_df.index and pdb2 in cos_sim_df.columns:
                same_cluster_similarities.append(cos_sim_df.loc[pdb1, pdb2])

# Calculate pairwise similarities for PDBs in different clusters
for i, pdb1 in enumerate(clustered_embeddings_df['PDB']):
    for pdb2 in clustered_embeddings_df['PDB'][i+1:]:
        cluster1 = clustered_embeddings_df[clustered_embeddings_df['PDB'] == pdb1]['Cluster'].values[0]
        cluster2 = clustered_embeddings_df[clustered_embeddings_df['PDB'] == pdb2]['Cluster'].values[0]
        if cluster1 != cluster2:
            if pdb1 in cos_sim_df.index and pdb2 in cos_sim_df.columns:
                different_cluster_similarities.append(cos_sim_df.loc[pdb1, pdb2])

# Bootstrap analysis
n_bootstrap = 1000
bootstrap_stats = {'same_cluster': [], 'different_cluster': []}

for _ in range(n_bootstrap):
    # Sample with replacement
    same_cluster_sample = np.random.choice(same_cluster_similarities, 
                                         size=len(same_cluster_similarities), 
                                         replace=True)
    different_cluster_sample = np.random.choice(different_cluster_similarities,
                                              size=len(different_cluster_similarities),
                                              replace=True)
    
    # Calculate statistics for this bootstrap sample
    bootstrap_stats['same_cluster'].append({
        'mean': np.mean(same_cluster_sample),
        'median': np.median(same_cluster_sample),
        'std': np.std(same_cluster_sample)
    })
    
    bootstrap_stats['different_cluster'].append({
        'mean': np.mean(different_cluster_sample),
        'median': np.median(different_cluster_sample),
        'std': np.std(different_cluster_sample)
    })

# Calculate confidence intervals
confidence_intervals = {
    'same_cluster': {
        'mean': np.percentile([x['mean'] for x in bootstrap_stats['same_cluster']], [2.5, 97.5]),
        'median': np.percentile([x['median'] for x in bootstrap_stats['same_cluster']], [2.5, 97.5]),
        'std': np.percentile([x['std'] for x in bootstrap_stats['same_cluster']], [2.5, 97.5])
    },
    'different_cluster': {
        'mean': np.percentile([x['mean'] for x in bootstrap_stats['different_cluster']], [2.5, 97.5]),
        'median': np.percentile([x['median'] for x in bootstrap_stats['different_cluster']], [2.5, 97.5]),
        'std': np.percentile([x['std'] for x in bootstrap_stats['different_cluster']], [2.5, 97.5])
    }
}

# Print bootstrap results
print("\nBootstrap Analysis Results (95% Confidence Intervals):")
print("\nSame Cluster:")
print(f"Mean: {confidence_intervals['same_cluster']['mean'][0]:.4f} - {confidence_intervals['same_cluster']['mean'][1]:.4f}")
print(f"Median: {confidence_intervals['same_cluster']['median'][0]:.4f} - {confidence_intervals['same_cluster']['median'][1]:.4f}")
print(f"Std: {confidence_intervals['same_cluster']['std'][0]:.4f} - {confidence_intervals['same_cluster']['std'][1]:.4f}")

print("\nDifferent Clusters:")
print(f"Mean: {confidence_intervals['different_cluster']['mean'][0]:.4f} - {confidence_intervals['different_cluster']['mean'][1]:.4f}")
print(f"Median: {confidence_intervals['different_cluster']['median'][0]:.4f} - {confidence_intervals['different_cluster']['median'][1]:.4f}")
print(f"Std: {confidence_intervals['different_cluster']['std'][0]:.4f} - {confidence_intervals['different_cluster']['std'][1]:.4f}")

# Plot the distribution of cosine similarities
plt.figure(figsize=(12, 8))
sns.kdeplot(1 - np.array(same_cluster_similarities), label='Same Cluster', color='#470659', linewidth=10)
sns.kdeplot(1 - np.array(different_cluster_similarities), label='Different Cluster', color='#FF8C00', linewidth=10)
plt.xlabel('1- Pearson')
plt.ylabel('Density')
plt.legend(fontsize=22, title_fontsize='22', loc='upper right', frameon=True, shadow=True)
plt.tight_layout()
plt.savefig('cosine_similarity_distribution.png', dpi=600)
plt.close()

##PERMUTATION TEST
# Merge with cluster information

merged_rmsf_df_clusters = pd.merge(merged_rmsf_df, clustered_embeddings_df, on='PDB', how='inner')

# Pivot s2calc data: rows = PDBs, columns = residues
RMSF_pivot = merged_rmsf_df_clusters.pivot_table(
    index='PDB', 
    columns='resi', 
    values='delta_RMSF'  # Adjust column name if different
).fillna(0)

# Get cluster labels aligned with PDBs
pdb_list = RMSF_pivot.index.tolist()
labels = np.array([
    clustered_embeddings_df[clustered_embeddings_df['PDB'] == pdb]['Cluster'].values[0] 
    for pdb in pdb_list
])

D_condensed = pdist(s2calc_pivot.values, metric='eucledian')
D = squareform(D_condensed)

# Collect within- and between-cluster distances
within = []
between = []
for i, j in combinations(range(len(pdb_list)), 2):
    if labels[i] == labels[j]:
        within.append(D[i, j])
    else:
        between.append(D[i, j])

within = np.array(within, dtype=float)
between = np.array(between, dtype=float)

# ----- CLIFF'S DELTA (effect size) -----
def cliffs_delta(a, b):
    """
    Cliff's delta effect size: probability(a>b) - probability(a<b)
    Works well for non-normal, different-sized samples
    |δ| < 0.147 = negligible, < 0.33 = small, < 0.474 = medium, >= 0.474 = large
    """
    a = np.asarray(a)
    b = np.asarray(b)
    greater = 0
    less = 0
    for x in a:
        greater += np.sum(x > b)
        less += np.sum(x < b)
    n = a.size * b.size
    return (greater - less) / n if n > 0 else np.nan

# ----- MANN-WHITNEY U TEST -----
u_stat, mw_p = mannwhitneyu(between, within, alternative='greater')  # H1: between > within
delta = cliffs_delta(between, within)

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


print(f"Distance metric: Eucledian")
print(f"\nPairwise comparisons:")
print(f"  Within-cluster pairs: {len(within)}")
print(f"  Between-cluster pairs: {len(between)}")
print(f"\nDistance statistics:")
print(f"  Within-cluster mean: {within.mean():.3f} (SD: {within.std():.3f})")
print(f"  Between-cluster mean: {between.mean():.3f} (SD: {between.std():.3f})")
print(f"  Observed difference (between - within): {obs_diff:.3f}")
print(f"\nStatistical tests:")
print(f"  Mann-Whitney U test (between > within):")
print(f"    U-statistic: {u_stat:.0f}")
print(f"    p-value: {mw_p:.3g}")
print(f"  Permutation test (5000 iterations):")
print(f"    p-value: {perm_p:.3g}")


# ----- VISUALIZATION -----
# Plot distributions
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Histogram of distances
axes[0].hist(within, bins=30, alpha=0.6, label='Within-cluster', color='darkviolet', density=True)
axes[0].hist(between, bins=30, alpha=0.6, label='Between-cluster', color='darkorange', density=True)
axes[0].axvline(within.mean(), color='darkviolet', linestyle='--', linewidth=5, label=f'Within mean: {within.mean():.2f}')
axes[0].axvline(between.mean(), color='darkorange', linestyle='--', linewidth=5, label=f'Between mean: {between.mean():.2f}')
axes[0].set_xlabel(f'Eucledian Distance', fontweight='bold')
axes[0].set_ylabel('Density', fontweight='bold')
axes[0].legend(fontsize=11)
axes[0].grid(alpha=0.3)

# Permutation test null distribution
axes[1].hist(perm_diffs, bins=50, alpha=0.7, color='gray', edgecolor='black')
axes[1].axvline(obs_diff, color='red', linestyle='--', linewidth=3, label=f'Observed: {obs_diff:.3f}')
axes[1].set_xlabel('Difference (Between - Within)', fontweight='bold')
axes[1].set_ylabel('Frequency', fontweight='bold')
axes[1].legend(fontsize=14)
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig('RMSF_permutation_test_results.png', dpi=300, bbox_inches='tight')
plt.close()


# Pivot the delta RMSF DataFrame to have residues in the columns and delta RMSF values in the rows
delta_rmsf_pivot = merged_rmsf_df.pivot(index='PDB', columns='resi', values='delta_RMSF').fillna(0)

# Plot the clustermap with delta RMSF values
clustermap = sns.clustermap(delta_rmsf_pivot, cmap="magma", figsize=(14, 10))
plt.xlabel("Residue", fontsize=16)
plt.ylabel("PDB Structures", fontsize=16)
plt.tight_layout()
plt.savefig('delta_rmsf_clustermap.png', dpi=300)
plt.close()

# Get the order of rows and columns from the clustermap
row_order = clustermap.dendrogram_row.reordered_ind
col_order = clustermap.dendrogram_col.reordered_ind

# Reorder the DataFrame according to the clustermap order
ordered_delta_rmsf_pivot = delta_rmsf_pivot.iloc[row_order, col_order]

# Output the ordered DataFrame to a CSV file
ordered_delta_rmsf_pivot.to_csv('ordered_delta_rmsf_clustermap_data.csv')
# Calculate the average value for each column in ordered_delta_rmsf_pivot
average_values = ordered_delta_rmsf_pivot.median()

# Find residues where the average value is greater than 0.25
residues_increasing_flexibility = average_values[average_values < -0.25].index.tolist()

# Print the residues
print("Residues that always increase flexibility:", residues_increasing_flexibility)



# Determine correlation across residues
pivot_RMSF = merged_rmsf_df.pivot_table(index=["resi"], 
                       columns="PDB_ensemble", 
                       values="delta_RMSF")

# Compute correlation across residues
corr_matrix = pivot_RMSF.corrwith(pivot_RMSF, axis=1)  
corr_matrix = pivot_RMSF.T.corr()

# Plot heatmap
plt.figure(figsize=(14, 12))
sns.heatmap(corr_matrix, cmap="inferno", center=0, square=True, cbar_kws={"label": "Correlation"})
plt.title("Correlation of RMSF across residues", fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig(f'rmsf_correlation_all.png', dpi=300)
plt.close()

# Get the order of rows and columns from the clustermap
row_order = clustermap.dendrogram_row.reordered_ind
col_order = clustermap.dendrogram_col.reordered_ind

# Reorder the correlation matrix according to the clustermap order
ordered_corr_matrix = corr_matrix.iloc[row_order, col_order]

# Output the ordered correlation matrix to a CSV file
ordered_corr_matrix.to_csv('ordered_rmsf_correlation_matrix.csv')

# Calculate the correlation matrix for delta RMSF values
delta_rmsf_correlation = delta_rmsf_pivot.corr()

# Output the correlation matrix to a CSV file
delta_rmsf_correlation.to_csv('delta_rmsf_correlation_matrix.csv')

# Plot the heatmap of the correlation matrix
plt.figure(figsize=(14, 12))
sns.heatmap(delta_rmsf_correlation, cmap="cubehelix", center=0, square=True, cbar_kws={"label": "Correlation"})
plt.title("Correlation Matrix of Delta RMSF Values", fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Residue', fontsize=16)
plt.ylabel('Residue', fontsize=16)
plt.tight_layout()
plt.savefig('delta_rmsf_correlation_heatmap.png', dpi=300)
plt.close()
