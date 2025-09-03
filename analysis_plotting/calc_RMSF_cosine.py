import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt

def pairwise_similarity(
    df: pd.DataFrame,
    id_col
    x_col
    y_col
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


# Cosine similarity matrix across all PDBs
print('COSINE:')
cos_sim_df = pairwise_similarity(merged_rmsf_df, PDB_ensemble, resi, delta_RMSF)
log_cos_sim_df = np.log(cos_sim_df)
print(cos_sim_df.round(3))
cos_sim_df.to_csv("cosine_similarity_deltaRMSF.csv")

# Quick heatmap 
sns.clustermap(log_cos_sim_df, cmap="coolwarm", figsize=(6, 5))
plt.savefig('cosine_similarity_clustermap.png', dpi=300)
plt.close()
