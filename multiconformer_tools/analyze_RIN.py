# Load and analyze the uploaded RIN differential edges spreadsheet
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from math import log10

# Load
path = Path("RIN_diff_edges.csv")
xls = pd.ExcelFile(path)
sheet_names = xls.sheet_names
df = pd.read_excel(path, sheet_name=0)

# Standardize column names (strip spaces, lower)
df.columns = [c.strip().lower() for c in df.columns]

# Quick sanity: ensure expected columns exist
expected = {"cluster","node_u","node_v","freq_other","delta_freq","p_value_fisher_two_sided"}
missing = expected - set(df.columns)
summary_notes = {}
summary_notes['sheet_names'] = sheet_names
summary_notes['missing_columns'] = list(missing)

# Add helper columns
df['edge'] = df['node_u'].astype(str) + " — " + df['node_v'].astype(str)
df['direction'] = np.where(df['delta_freq']>=0, 'higher_in_cluster', 'lower_in_cluster')
df['abs_delta'] = df['delta_freq'].abs()
df['neglog10_p'] = -np.log10(df['p_value_fisher_two_sided'].replace(0, np.nextafter(0,1)))

# Basic descriptive stats
by_cluster_counts = df.groupby('cluster').size().rename('n_edges').reset_index()
by_cluster_sig = (df[df['p_value_fisher_two_sided'] < 0.05]
                  .groupby('cluster').size().rename('n_sig_edges').reset_index())
by_cluster = by_cluster_counts.merge(by_cluster_sig, on='cluster', how='left').fillna(0)
by_cluster['frac_sig'] = by_cluster['n_sig_edges']/by_cluster['n_edges']

# Strong-effect, significant edges (|delta|>=0.5 and p<0.05)
strong = df[(df['abs_delta']>=0.5) & (df['p_value_fisher_two_sided']<0.05)].copy()
strong_top10 = strong.sort_values(['abs_delta','p_value_fisher_two_sided'], ascending=[False,True]).head(10)

# Per-cluster “net shift” (sum of delta across significant edges), and counts by direction
sig = df[df['p_value_fisher_two_sided'] < 0.05].copy()
net_shift = sig.groupby('cluster')['delta_freq'].sum().rename('net_delta_sig').reset_index()
dir_counts = (sig.groupby(['cluster','direction']).size()
                .unstack(fill_value=0)
                .reset_index())
per_cluster_summary = by_cluster.merge(net_shift, on='cluster', how='left').merge(dir_counts, on='cluster', how='left').fillna(0)

# Volcano-like plot: delta_freq vs -log10(p)
plt.figure()
plt.scatter(df['delta_freq'], df['neglog10_p'], s=10, alpha=0.6)
plt.axvline(0, linestyle='--')
plt.xlabel("Delta frequency (cluster − others)")
plt.ylabel("−log10(p)")
plt.tight_layout()
plt.savefig("volcano_all_clusters.png", dpi=200)
plt.show()

# Bar plot: fraction of significant edges per cluster
plt.figure()
order = by_cluster.sort_values('frac_sig', ascending=False)['cluster']
plt.bar(by_cluster.set_index('cluster').loc[order].index.astype(str),
        by_cluster.set_index('cluster').loc[order]['frac_sig'])
plt.xlabel("Cluster")
plt.ylabel("Fraction of significant edges (p<0.05)")
plt.tight_layout()
plt.savefig(out_dir/"fraction_significant_by_cluster.png", dpi=200)
plt.show()

# Bar plot: net shift among significant edges per cluster
plt.figure()
order2 = per_cluster_summary.sort_values('net_delta_sig', ascending=False)['cluster']
plt.bar(per_cluster_summary.set_index('cluster').loc[order2].index.astype(str),
        per_cluster_summary.set_index('cluster').loc[order2]['net_delta_sig'])
plt.xlabel("Cluster")
plt.ylabel("Sum of delta frequency (significant edges only)")
plt.tight_layout()
plt.savefig("net_shift_by_cluster.png", dpi=200)
plt.show()

summary_notes
