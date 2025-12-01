#!/usr/bin/env python3

"""
H-bond graph analysis + mechanistic water replaceability + PATHWAY DIFFERENCES.
Outputs (new bits are marked [NEW]):
  - per_pdb_graph_metrics.csv (+ mechanistic columns when waters included)
  - per_pdb_node_centralities.csv
  - cluster_graph_summary.csv
  - top_central_residues_by_cluster.csv
  - per_water_replaceability_scores.csv
  - cluster_fraction_replaceable_waters.csv
  - [NEW] per_pdb_shortest_paths.csv               # all P–L shortest paths (<= max_len)
  - [NEW] per_pdb_edge_usage_from_paths.csv        # edges touched by any shortest path (binary)
  - [NEW] cluster_edge_freq_from_paths.csv         # edge frequency per cluster
  - [NEW] cluster_edge_differential.csv            # log2FC (+ FET q-values if SciPy)
  - [NEW] cluster_path_motif_counts.csv            # motif counts per cluster
  - [NEW] consensus_<Cluster>.graphml              # per-cluster consensus pathway graph
  - Plots incl. differential edge heatmap, motif distributions

Requirements: pandas, numpy, networkx, matplotlib, seaborn
Optional: scipy (for fisher_exact + BH-FDR)
"""

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

# Optional stats
try:
    from scipy.stats import fisher_exact
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False

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

# ----------------------------- Small helpers ---------------------------------

AA3 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "SEC","PYL"
}

def pdb_stem(name: str) -> str:
    return (name
            .replace("-pandda-model_refine_001.updated.pdb_hbonds", "")
            .replace("-pandda-model_refine_001.updated.pdb_summary", ""))

def load_clusters(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["PDB"] = df["PDB"].apply(lambda x: Path(x).stem.replace("_refitted.pdbredist_norm_lig.pdb",""))
    df["PDB"] = df["PDB"].str.replace('_qFit_refitted.pdbredist_norm_lig', '', regex=False)
    return df[["PDB","Cluster"]].drop_duplicates()

def load_subset_list(path: Path | None) -> set | None:
    if not path: return None
    keep = set()
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                keep.add(s)
    return keep

def load_subset_by_cluster(path: Path | None) -> dict | None:
    if not path: return None
    df = pd.read_csv(path)  # columns: Cluster,res (e.g., A:5)
    out = {}
    for c, sub in df.groupby("Cluster"):
        out[c] = set(map(str, sub["res"].tolist()))
    return out

def node_type_from_resname(resname: str, water_resnames: set[str]) -> str:
    r = (resname or "").upper()
    if r in water_resnames:
        return "WATER"
    if r in AA3:
        return "PROTEIN"
    return "LIGAND"

# ----------------------------- Graph building --------------------------------

def build_graph(hb: pd.DataFrame,
                subset: set | None,
                include_waters: bool,
                water_resnames: set[str]) -> tuple[nx.Graph, set]:
    """
    Build undirected weighted graph of residue nodes (chain:resi).
    If include_waters=False, rows with donor/acceptor resname in water_resnames are dropped.
    Returns (G, water_nodes_set).
    """
    hb = hb.copy()

    # Filter out waters unless requested
    if not include_waters:
        mask_d = ~hb["donor_resname"].astype(str).str.upper().isin(water_resnames)
        mask_a = ~hb["acceptor_resname"].astype(str).str.upper().isin(water_resnames)
        hb = hb[mask_d & mask_a]

    # Node labels
    hb["u"] = hb["donor_chain"].astype(str) + ":" + hb["donor_resi"].astype(str)
    hb["v"] = hb["acceptor_chain"].astype(str) + ":" + hb["acceptor_resi"].astype(str)
    hb = hb[hb["u"] != hb["v"]]

    # Track which nodes are waters (useful for annotation)
    water_nodes = set()
    water_nodes.update(
        hb.loc[hb["donor_resname"].astype(str).str.upper().isin(water_resnames), "u"].tolist()
    )
    water_nodes.update(
        hb.loc[hb["acceptor_resname"].astype(str).str.upper().isin(water_resnames), "v"].tolist()
    )

    # Collapse parallel edges -> undirected with weight = count
    pairs = np.sort(hb[["u","v"]].to_numpy(), axis=1)
    dfc = pd.DataFrame(pairs, columns=["u","v"]).value_counts().reset_index(name="w")

    G = nx.Graph()
    if subset is not None:
        G.add_nodes_from(subset)
        for _, r in dfc.iterrows():
            u, v, w = r["u"], r["v"], int(r["w"])
            if u in subset and v in subset:
                G.add_edge(u, v, weight=w)
    else:
        for _, r in dfc.iterrows():
            G.add_edge(r["u"], r["v"], weight=int(r["w"]))

    return G, water_nodes

def graph_metrics(G: nx.Graph) -> dict:
    n = G.number_of_nodes()
    if n == 0:
        return dict(n_nodes=0,n_edges=0,density=0,n_components=0,
                    giant_component_size=0,avg_clustering=np.nan,
                    avg_shortest_path_len_gc=np.nan,algebraic_connectivity_gc=np.nan)
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    GC = G.subgraph(comps[0]).copy()
    out = dict(
        n_nodes=n,
        n_edges=G.number_of_edges(),
        density=nx.density(G) if n>1 else 0.0,
        n_components=len(comps),
        giant_component_size=GC.number_of_nodes(),
        avg_clustering=nx.average_clustering(G) if n>1 else 0.0,
        avg_shortest_path_len_gc=np.nan,
        algebraic_connectivity_gc=np.nan
    )
    if GC.number_of_nodes()>1 and GC.number_of_edges()>0:
        try: out["avg_shortest_path_len_gc"] = nx.average_shortest_path_length(GC)
        except: pass
        try: out["algebraic_connectivity_gc"] = nx.algebraic_connectivity(GC)
        except: pass
    return out

def node_centrality(G: nx.Graph) -> pd.DataFrame:
    if G.number_of_nodes()==0:
        return pd.DataFrame(columns=["res","degree","weighted_degree","betweenness","closeness"])
    deg = dict(G.degree())
    wdeg = dict(G.degree(weight="weight"))
    btw = nx.betweenness_centrality(G, normalized=True)
    clo = nx.closeness_centrality(G)
    return pd.DataFrame({
        "res": list(G.nodes()),
        "degree": [deg[n] for n in G.nodes()],
        "weighted_degree": [wdeg[n] for n in G.nodes()],
        "betweenness": [btw[n] for n in G.nodes()],
        "closeness": [clo[n] for n in G.nodes()],
    })

# ---------------------- Mechanistic water analysis ----------------------------

def water_bridge_table(G: nx.Graph, node_types: dict, hb_geo: pd.DataFrame | None = None) -> pd.DataFrame:
    rows = []
    if G.number_of_nodes()==0:
        return pd.DataFrame(columns=[
            "water_node","degree","weighted_degree","betweenness","closeness",
            "clustering","is_articulation","nP","nL","nW","triangles",
            "bridge_like","redundancy","replaceability_score"
        ])
    deg = dict(G.degree())
    wdeg = dict(G.degree(weight="weight"))
    btw = nx.betweenness_centrality(G, normalized=True)
    clo = nx.closeness_centrality(G)
    try:
        clus = nx.clustering(G, weight=None)
    except Exception:
        clus = {n: np.nan for n in G.nodes()}
    try:
        arts = set(nx.articulation_points(G))
    except Exception:
        arts = set()

    # optional geometry
    geom_min_dist = {}
    geom_max_angle = {}
    if hb_geo is not None:
        cols = set([c.lower() for c in hb_geo.columns])
        dist_col = next((c for c in ["da_dist","distance","d_a_dist","donor_acceptor_dist"] if c in cols), None)
        ang_col  = next((c for c in ["angle","dha_angle","donor_hydrogen_acceptor_angle"] if c in cols), None)

        def node_series(side):
            if side == "donor":
                return hb_geo["donor_chain"].astype(str) + ":" + hb_geo["donor_resi"].astype(str)
            return hb_geo["acceptor_chain"].astype(str) + ":" + hb_geo["acceptor_resi"].astype(str)

        if dist_col:
            for u in node_series("donor"):
                if node_types.get(u)=="WATER":
                    geom_min_dist[u] = min(geom_min_dist.get(u, np.inf), float('inf'))
            for v in node_series("acceptor"):
                if node_types.get(v)=="WATER":
                    geom_min_dist[v] = min(geom_min_dist.get(v, np.inf), float('inf'))

    for n in G.nodes():
        if node_types.get(n) != "WATER":
            continue
        neigh = list(G.neighbors(n))
        nP = sum(node_types.get(x)=="PROTEIN" for x in neigh)
        nL = sum(node_types.get(x)=="LIGAND"  for x in neigh)
        nW = sum(node_types.get(x)=="WATER"   for x in neigh)

        triangles = 0
        for i in range(len(neigh)):
            for j in range(i+1, len(neigh)):
                if G.has_edge(neigh[i], neigh[j]):
                    triangles += 1

        bridge_like = int(dict(G.degree())[n]==2 and ((nL==1 and nP==1) or (nP==2) or (nL==2)))
        redundancy = clus.get(n, 0.0) if clus is not None else 0.0

        replaceability = (
            1.00*bridge_like +
            0.75*(1.0 if dict(G.degree())[n]==2 else 0.0) +
            0.50*nx.betweenness_centrality(G, normalized=True).get(n,0.0) -
            0.75*(1.0 if n in arts else 0.0) -
            0.50*min(1.0, redundancy)
        )

        rows.append(dict(
            water_node=n,
            degree=dict(G.degree())[n],
            weighted_degree=dict(G.degree(weight="weight"))[n],
            betweenness=nx.betweenness_centrality(G, normalized=True).get(n, np.nan),
            closeness=nx.closeness_centrality(G).get(n, np.nan),
            clustering=redundancy,
            is_articulation=int(n in arts),
            nP=nP, nL=nL, nW=nW, triangles=triangles,
            bridge_like=bridge_like,
            redundancy=redundancy,
            replaceability_score=float(replaceability)
        ))
    return pd.DataFrame(rows)

def count_water_motifs(G: nx.Graph, node_types: dict) -> dict:
    counts = dict(PWP=0, LWP=0, LWL=0, WWW_edges=0)
    for n in G.nodes():
        if node_types.get(n)!="WATER": 
            continue
        neigh = list(G.neighbors(n))
        for i in range(len(neigh)):
            for j in range(i+1,len(neigh)):
                ti, tj = node_types.get(neigh[i]), node_types.get(neigh[j])
                if (ti=="PROTEIN" and tj=="PROTEIN"):
                    counts["PWP"] += 1
                elif (ti=="LIGAND" and tj=="LIGAND"):
                    counts["LWL"] += 1
                elif {"LIGAND","PROTEIN"} == {ti,tj}:
                    counts["LWP"] += 1
        for m in neigh:
            if node_types.get(m)=="WATER":
                counts["WWW_edges"] += 1
    counts["WWW_edges"] //= 2
    return counts

def simulate_water_removal_effects(G: nx.Graph, node_types: dict) -> dict:
    prots = [n for n,t in node_types.items() if t=="PROTEIN" and n in G]
    ligs  = [n for n,t in node_types.items() if t=="LIGAND"  and n in G]
    if len(prots)==0 or len(ligs)==0 or G.number_of_edges()==0:
        return dict(n_PL_paths_with_water=0, n_PL_paths_no_water=0,
                    avg_PL_pathlen_with_water=np.nan, avg_PL_pathlen_no_water=np.nan)

    n_paths_w, lens_w = 0, []
    for p in prots:
        for l in ligs:
            if nx.has_path(G, p, l):
                n_paths_w += 1
                try:
                    lens_w.append(nx.shortest_path_length(G, p, l))
                except Exception:
                    pass

    nodes_noW = [n for n in G.nodes() if node_types.get(n)!="WATER"]
    G_noW = G.subgraph(nodes_noW).copy()
    n_paths_now, lens_now = 0, []
    if G_noW.number_of_nodes()>0:
        for p in prots:
            if p not in G_noW: 
                continue
            for l in ligs:
                if l not in G_noW:
                    continue
                if nx.has_path(G_noW, p, l):
                    n_paths_now += 1
                    try:
                        lens_now.append(nx.shortest_path_length(G_noW, p, l))
                    except Exception:
                        pass

    return dict(
        n_PL_paths_with_water=n_paths_w,
        n_PL_paths_no_water=n_paths_now,
        avg_PL_pathlen_with_water=(np.mean(lens_w) if lens_w else np.nan),
        avg_PL_pathlen_no_water=(np.mean(lens_now) if lens_now else np.nan)
    )

# ----------------------- PATHWAY extraction & diffs ---------------------------

def path_type_sequence(path, node_types):
    """Compress a node path to a type motif, e.g. ['P','W','P','L']."""
    tmap = {'PROTEIN':'P','LIGAND':'L','WATER':'W'}
    return "-".join([tmap.get(node_types.get(n,'PROTEIN'),'P') for n in path])

def record_edge_set_for_path(path):
    """Return canonical undirected edges for that path as tuple-sorted pairs."""
    edges = []
    for i in range(len(path)-1):
        u, v = path[i], path[i+1]
        edges.append(tuple(sorted((u,v))))
    return edges

def shortest_PL_paths(G: nx.Graph, node_types: dict, max_len: int = 6):
    """
    For each ligand node, get shortest paths to all proteins within cutoff.
    Returns list of dicts: { 'src_lig', 'dst_prot', 'length', 'path', 'type_motif' }
    """
    ligs  = [n for n,t in node_types.items() if t=="LIGAND"  and n in G]
    prots = [n for n,t in node_types.items() if t=="PROTEIN" and n in G]
    rows = []
    if not ligs or not prots:
        return rows
    for l in ligs:
        # single-source shortest paths up to max_len
        ssp = nx.single_source_shortest_path(G, l, cutoff=max_len)
        for p in prots:
            if p not in ssp:
                continue
            path = ssp[p]
            if len(path)-1 <= max_len:
                rows.append(dict(
                    src_lig=l,
                    dst_prot=p,
                    length=len(path)-1,
                    path=path,
                    type_motif=path_type_sequence(path, node_types)
                ))
    return rows

# ----------------------------- Plotting helpers ------------------------------

def save_bar(df: pd.DataFrame, xcol: str, ycol: str, title: str, path: Path):
    d = df.sort_values(ycol, ascending=False)
    plt.figure(figsize=(10,5))
    plt.bar(d[xcol], d[ycol])
    plt.xlabel(xcol); plt.ylabel(ycol.replace("_", " ")); plt.title(title)
    plt.xticks(rotation=90); plt.tight_layout()
    plt.savefig(path, dpi=300); plt.close()

def save_scatter(x, y, xlabel, ylabel, title, labels, sizes, path: Path):
    plt.figure(figsize=(6,5))
    plt.scatter(x, y, s=sizes)
    for xi, yi, lab in zip(x, y, labels):
        plt.annotate(lab, (xi, yi), fontsize=8, xytext=(3,3), textcoords="offset points")
    plt.xlabel(xlabel); plt.ylabel(ylabel); plt.title(title)
    plt.tight_layout(); plt.savefig(path, dpi=300); plt.close()

def save_heatmap(matrix: np.ndarray, yticks, xticks, title: str, path: Path, cbar_label: str):
    df = pd.DataFrame(matrix, index=yticks, columns=xticks)
    has_nan = df.isna().any().any() or np.isinf(df.values).any()
    if has_nan:
        df_plot = df.fillna(0)
        plt.figure(figsize=(max(8, 0.5*len(xticks)), max(4, 0.4*len(yticks))))
        sns.heatmap(df_plot, cmap="magma_r", center=0, cbar_kws={'label': cbar_label},
                    xticklabels=True, yticklabels=True, annot=False, fmt='.2f')
        plt.xticks(rotation=90); plt.title(title)
        plt.tight_layout(); plt.savefig(path, dpi=300, bbox_inches='tight'); plt.close()
    else:
        g = sns.clustermap(df, cmap="magma_r", center=0,
                           figsize=(max(8, 0.5*len(xticks)), max(4, 0.4*len(yticks))),
                           cbar_kws={'label': cbar_label},
                           xticklabels=True, yticklabels=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.suptitle(title)
        plt.savefig(path, dpi=300, bbox_inches='tight'); plt.close()

# ---------------------------------- Main -------------------------------------

def main():
    ap = argparse.ArgumentParser(description="H-bond graph analysis + mechanistic waters + pathway diffs.")
    ap.add_argument("--hbond_dir", required=True, help="Folder with *_hbonds.csv")
    ap.add_argument("--clusters_csv", required=True, help="HBDScan_clusters_3ormore_normpdb.csv")
    ap.add_argument("--outdir", default="hbond_graph_results_paths")
    ap.add_argument("--subset_residues", default=None, help="Txt file with residues like 'A:5' (optional)")
    ap.add_argument("--subset_residues_by_cluster", default=None, help="CSV Cluster,res (optional)")
    # Waters & paths
    ap.add_argument("--include_waters", action="store_true", help="Include water molecules as nodes")
    ap.add_argument("--max_path_len", type=int, default=6, help="Max hops for shortest P–L paths")
    # RMSF (optional)
    ap.add_argument("--rmsf_csv", default="avg_delta_rmsf.csv",
                    help="Optional CSV with columns PDB_ensemble or PDB, delta_RMSF")
    args = ap.parse_args()

    hbond_dir = Path(args.hbond_dir)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    clusters = load_clusters(Path(args.clusters_csv))
    pdb2cl = dict(zip(clusters.PDB, clusters.Cluster))
    print(f"Loaded {len(pdb2cl)} PDB entries from cluster CSV")

    global_subset = load_subset_list(Path(args.subset_residues)) if args.subset_residues else None
    subset_by_cluster = load_subset_by_cluster(Path(args.subset_residues_by_cluster)) if args.subset_residues_by_cluster else None

    graph_rows, node_rows = [], []
    per_water_all = []
    per_pdb_paths_rows = []
    per_pdb_edge_usage_rows = []

    # --- build graphs + metrics ---
    hbond_files = sorted(hbond_dir.glob("*_hbonds.csv"))
    print(f"Found {len(hbond_files)} hbond files")

    for f in hbond_files:
        pdb = str(f).split('/')[-1][:5]
        if pdb not in pdb2cl:
            continue
        cl = pdb2cl[pdb]
        hb = pd.read_csv(f)
        # Remove any rows where donor_chain or acceptor_chain is 'B'
        hb = hb[~((hb["donor_chain"] == "B") | (hb["acceptor_chain"] == "B"))].reset_index(drop=True)

        # normalize resnames for water identification
        hb["donor_resname"] = hb["donor_resname"].astype(str).str.upper()
        hb["acceptor_resname"] = hb["acceptor_resname"].astype(str).str.upper()

        subset = subset_by_cluster.get(cl) if subset_by_cluster else global_subset
        G, water_nodes = build_graph(hb, subset,
                                     include_waters=args.include_waters,
                                     water_resnames=['HOH'])

        # node typing (protein/ligand/water)
        node_types = {}
        u_nodes = hb["donor_chain"].astype(str)+":"+hb["donor_resi"].astype(str)
        for resname, node in zip(hb["donor_resname"].astype(str), u_nodes):
            t = node_type_from_resname(resname, 'HOH')
            if node not in node_types or t == "WATER":
                node_types[node] = t
        v_nodes = hb["acceptor_chain"].astype(str)+":"+hb["acceptor_resi"].astype(str)
        for resname, node in zip(hb["acceptor_resname"].astype(str), v_nodes):
            t = node_type_from_resname(resname, 'HOH')
            if node not in node_types or t == "WATER":
                node_types[node] = t
        node_types_G = {n: node_types.get(n, "PROTEIN") for n in G.nodes()}

        # core metrics
        gm = graph_metrics(G); gm.update(PDB=pdb, Cluster=cl)

        # mechanistic analysis (waters)
        motifs = {}
        mechanistic = {}
        water_df = pd.DataFrame()
        if args.include_waters:
            water_df = water_bridge_table(G, node_types_G, hb_geo=hb)
            if not water_df.empty:
                water_df.insert(0, "PDB", pdb)
                water_df.insert(1, "Cluster", cl)
                water_df.sort_values("replaceability_score", ascending=False, inplace=True)
                per_water_all.append(water_df)
                water_df.to_csv(outdir/f"{pdb}_top_replaceable_waters.csv", index=False)

            motifs = count_water_motifs(G, node_types_G)
            mechanistic = simulate_water_removal_effects(G, node_types_G)

        gm.update({f"motif_{k}": v for k,v in motifs.items()})
        gm.update(mechanistic)
        graph_rows.append(gm)

        # node centralities
        ndf = node_centrality(G)
        if not ndf.empty:
            ndf.insert(0, "PDB", pdb)
            ndf.insert(1, "Cluster", cl)
            ndf["is_water"] = ndf["res"].isin(water_nodes).astype(int)
            node_rows.append(ndf)

        # ---------------- PATHWAYS (shortest P–L) ----------------
        path_rows = shortest_PL_paths(G, node_types_G, max_len=args.max_path_len)
        if path_rows:
            # expand to DataFrame
            tmp = pd.DataFrame([{
                "PDB": pdb,
                "Cluster": cl,
                "src_lig": r["src_lig"],
                "dst_prot": r["dst_prot"],
                "length": r["length"],
                "type_motif": r["type_motif"],
                "path_nodes": "->".join(r["path"])
            } for r in path_rows])
            per_pdb_paths_rows.append(tmp)

            # mark which edges were used by any shortest path (binary per PDB)
            used_edges = set()
            for r in path_rows:
                for e in record_edge_set_for_path(r["path"]):
                    used_edges.add(e)
            for (u,v) in used_edges:
                per_pdb_edge_usage_rows.append(dict(PDB=pdb, Cluster=cl, u=u, v=v, used_in_any_PL_shortest_path=1))

    # --- write CSVs (core) ---
    graph_df = (pd.DataFrame(graph_rows).sort_values(["Cluster","PDB"])
                if graph_rows else
                pd.DataFrame(columns=["PDB","Cluster","n_nodes","n_edges","density","n_components",
                                      "giant_component_size","avg_clustering","avg_shortest_path_len_gc",
                                      "algebraic_connectivity_gc"]))
    print(graph_df.head())
    graph_df.to_csv(outdir/"per_pdb_graph_metrics.csv", index=False)

    node_df = (pd.concat(node_rows, ignore_index=True)
               if node_rows else
               pd.DataFrame(columns=["PDB","Cluster","res","degree","weighted_degree","betweenness","closeness","is_water"]))
    node_df.to_csv(outdir/"per_pdb_node_centralities.csv", index=False)

    # cluster summary (means)
    if not graph_df.empty:
        clust_sum = (graph_df.groupby("Cluster")[[
            "n_nodes","n_edges","density","n_components",
            "giant_component_size","avg_clustering",
            "avg_shortest_path_len_gc","algebraic_connectivity_gc"
        ]].mean().reset_index())
        clust_sum.to_csv(outdir/"cluster_graph_summary.csv", index=False)
    else:
        clust_sum = pd.DataFrame()

    # -- per-PDB PATH lists
    paths_df = pd.concat(per_pdb_paths_rows, ignore_index=True) if per_pdb_paths_rows else pd.DataFrame(
        columns=["PDB","Cluster","src_lig","dst_prot","length","type_motif","path_nodes"])
    paths_df.to_csv(outdir/"per_pdb_shortest_paths.csv", index=False)

    edge_usage_df = pd.DataFrame(per_pdb_edge_usage_rows) if per_pdb_edge_usage_rows else pd.DataFrame(
        columns=["PDB","Cluster","u","v","used_in_any_PL_shortest_path"])
    edge_usage_df.to_csv(outdir/"per_pdb_edge_usage_from_paths.csv", index=False)

    # ---------------------- Cluster DIFFERENCES from pathways ------------------
    if not edge_usage_df.empty:
        # Frequency of each edge per cluster (normalize by # PDBs in cluster)
        pdbs_per_cluster = (edge_usage_df[["PDB","Cluster"]].drop_duplicates()
                            .groupby("Cluster")["PDB"].nunique().rename("n_pdbs")).reset_index()
        edge_counts = (edge_usage_df.groupby(["Cluster","u","v"])["used_in_any_PL_shortest_path"]
                       .sum().rename("n_pdbs_with_edge").reset_index())
        edge_freq = edge_counts.merge(pdbs_per_cluster, on="Cluster", how="left")
        edge_freq["freq"] = edge_freq["n_pdbs_with_edge"] / edge_freq["n_pdbs"].clip(lower=1)
        edge_freq.to_csv(outdir/"cluster_edge_freq_from_paths.csv", index=False)

        # Differential edge signal: compare each cluster vs the rest (one-vs-all)
        all_clusters = sorted(edge_freq["Cluster"].unique().tolist())
        diffs = []
        for cl in all_clusters:
            cl_pdbs = pdbs_per_cluster.loc[pdbs_per_cluster["Cluster"]==cl, "n_pdbs"].values[0]
            rest_pdbs = pdbs_per_cluster["n_pdbs"].sum() - cl_pdbs
            sub = edge_freq[edge_freq["Cluster"]==cl][["u","v","n_pdbs_with_edge","n_pdbs"]].copy()
            sub = sub.rename(columns={"n_pdbs_with_edge":"k_cl", "n_pdbs":"n_cl"})
            # get rest counts for same edges
            rest = (edge_freq[edge_freq["Cluster"]!=cl]
                    .groupby(["u","v"])["n_pdbs_with_edge"].sum().rename("k_rest").reset_index())
            merged = sub.merge(rest, on=["u","v"], how="outer").fillna(0)
            merged["n_rest"] = rest_pdbs
            # compute log2FC & Fisher exact if available
            eps = 1e-9
            merged["freq_cl"] = merged["k_cl"] / merged["n_cl"].clip(lower=1)
            merged["freq_rest"] = merged["k_rest"] / merged["n_rest"].clip(lower=1)
            merged["log2FC_cl_vs_rest"] = np.log2( (merged["freq_cl"]+eps) / (merged["freq_rest"]+eps) )
            if SCIPY_OK:
                pvals = []
                for _, r in merged.iterrows():
                    table = [[int(r["k_cl"]), int(r["n_cl"]-r["k_cl"])],
                             [int(r["k_rest"]), int(r["n_rest"]-r["k_rest"])]]
                    _, p = fisher_exact(table, alternative="two-sided")
                    pvals.append(p)
                merged["pval"] = pvals
            merged["Cluster"] = cl
            diffs.append(merged)

        diff_df = pd.concat(diffs, ignore_index=True) if diffs else pd.DataFrame()
        if SCIPY_OK and not diff_df.empty:
            # BH-FDR within each cluster
            qvals = []
            for cl, sub in diff_df.groupby("Cluster"):
                m = len(sub)
                order = np.argsort(sub["pval"].values)
                q = np.empty(m, dtype=float); q.fill(np.nan)
                ranked = np.arange(1, m+1)
                p_sorted = sub["pval"].values[order]
                q_sorted = p_sorted * m / ranked
                # monotone
                for i in range(m-2, -1, -1):
                    q_sorted[i] = min(q_sorted[i], q_sorted[i+1])
                q[order] = np.minimum(1.0, q_sorted)
                qvals.append(pd.Series(q, index=sub.index))
            diff_df["qval"] = pd.concat(qvals).sort_index()
        diff_df.to_csv(outdir/"cluster_edge_differential.csv", index=False)

        # Plot: heatmap of top differential edges per cluster by |log2FC|
        if not diff_df.empty:
            topN = 30
            top_rows = []
            for cl, sub in diff_df.groupby("Cluster"):
                sub2 = sub.copy()
                if "qval" in sub2.columns:
                    sub2 = sub2.sort_values(["qval","log2FC_cl_vs_rest"], ascending=[True, False])
                else:
                    sub2 = sub2.sort_values("log2FC_cl_vs_rest", ascending=False)
                top_rows.append(sub2.head(topN))
            tops = pd.concat(top_rows, ignore_index=True) if top_rows else pd.DataFrame()
            if not tops.empty:
                tops["edge"] = tops.apply(lambda r: f"{r['u']}–{r['v']}", axis=1)
                piv = tops.pivot_table(index="Cluster", columns="edge", values="log2FC_cl_vs_rest", aggfunc="first")
                save_heatmap(piv.values, yticks=piv.index.tolist(), xticks=piv.columns.tolist(),
                             title="Top differential edges per cluster (from P–L shortest paths)",
                             path=outdir/"differential_edges_heatmap.png",
                             cbar_label="log2FC (cluster vs rest)")

        # Export per-cluster consensus graphs (GraphML) weighted by edge frequency
        for cl, sub in edge_freq.groupby("Cluster"):
            Gc = nx.Graph()
            for _, r in sub.iterrows():
                u, v = r["u"], r["v"]
                Gc.add_edge(u, v, weight=float(r["freq"]))
            nx.write_graphml(Gc, outdir/f"consensus_{cl}.graphml")

    # Path motif distributions
    if not paths_df.empty:
        motif_counts = (paths_df.groupby(["Cluster","type_motif"])["PDB"]
                        .count().rename("n_paths").reset_index())
        motif_counts.to_csv(outdir/"cluster_path_motif_counts.csv", index=False)

        # Stacked bar of most common motifs
        # pick top K motifs overall
        K = 10
        top_motifs = (motif_counts.groupby("type_motif")["n_paths"].sum()
                      .sort_values(ascending=False).head(K).index.tolist())
        sub = motif_counts[motif_counts["type_motif"].isin(top_motifs)].copy()
        piv = sub.pivot_table(index="Cluster", columns="type_motif", values="n_paths", aggfunc="sum").fillna(0)
        piv = piv.div(piv.sum(axis=1).replace(0,np.nan), axis=0)  # fraction
        piv.plot(kind="bar", stacked=True, figsize=(14,7))
        plt.ylabel("Fraction of shortest P–L paths")
        plt.title("Pathway motif distribution by cluster")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(outdir/"pathway_motif_stackedbar.png", dpi=300, bbox_inches="tight")
        plt.close()

    # ------------------------------ PLOTS (core) ------------------------------
    print(clust_sum.head())
    if not clust_sum.empty:
        clust_sum["Cluster"] = clust_sum["Cluster"].astype(str)

        plt.figure(figsize=(12, 8))
        cluster_medians = graph_df.groupby('Cluster')['n_edges'].median().sort_values()
        cluster_order = cluster_medians.index.tolist()
        ax = sns.boxplot(
            data=graph_df, 
            x='Cluster', 
            y='n_edges', 
            palette='magma_r', 
            order=cluster_order
        )
        plt.xlabel('Cluster')
        plt.ylabel('Number of Hydrogen Bonds')
        plt.xticks(rotation=45, ha='right')
        plt.yticks()
        # Add horizontal grid lines for better readability
        ax.yaxis.grid(True, which='major', linestyle='-', linewidth=1, color='gray', alpha=0.5)
        plt.tight_layout()
        plt.savefig(outdir/'hbond_edges_by_cluster_boxplot.png', dpi=300, bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(12, 8))
        cluster_medians = graph_df.groupby('Cluster')['density'].median().sort_values()
        cluster_order = cluster_medians.index.tolist()
        ax = sns.boxplot(
            data=graph_df, 
            x='Cluster', 
            y='density', 
            palette='magma_r', 
            order=cluster_order
        )
        plt.xlabel('Cluster')
        plt.ylabel('Hydrogen Bond Density')
        plt.xticks(rotation=45, ha='right')
        plt.yticks()
        # Add horizontal grid lines for better readability
        ax.yaxis.grid(True, which='major', linestyle='-', linewidth=1, color='gray', alpha=0.5)
        plt.tight_layout()
        plt.savefig(outdir/'hbond_density_by_cluster_boxplot.png', dpi=300, bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(12, 8))
        cluster_medians = graph_df.groupby('Cluster')['giant_component_size'].median().sort_values()
        cluster_order = cluster_medians.index.tolist()
        ax = sns.boxplot(
            data=graph_df, 
            x='Cluster', 
            y='giant_component_size', 
            palette='magma_r', 
            order=cluster_order
        )
        plt.xlabel('Cluster')
        plt.ylabel('Large Component')
        plt.xticks(rotation=45, ha='right')
        plt.yticks()
        # Add horizontal grid lines for better readability
        ax.yaxis.grid(True, which='major', linestyle='-', linewidth=1, color='gray', alpha=0.5)
        plt.tight_layout()
        plt.savefig(outdir/'hbond_giant_component_size_by_cluster_boxplot.png', dpi=300, bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(12, 8))
        cluster_medians = graph_df.groupby('Cluster')['avg_shortest_path_len_gc'].median().sort_values()
        cluster_order = cluster_medians.index.tolist()
        ax = sns.boxplot(
            data=graph_df, 
            x='Cluster', 
            y='avg_shortest_path_len_gc', 
            palette='magma_r', 
            order=cluster_order
        )
        plt.xlabel('Cluster')
        plt.ylabel('Average Shortest Length')
        plt.xticks(rotation=45, ha='right')
        plt.yticks()
        # Add horizontal grid lines for better readability
        ax.yaxis.grid(True, which='major', linestyle='-', linewidth=1, color='gray', alpha=0.5)
        plt.tight_layout()
        plt.savefig(outdir/'hbond_shortest_length_by_cluster_boxplot.png', dpi=300, bbox_inches='tight')
        plt.close()


        sizes = 20 + 2.0 * clust_sum["giant_component_size"]
        save_scatter(clust_sum["density"], clust_sum["avg_shortest_path_len_gc"],
                     "Graph density", "Avg shortest path length (GC)",
                     "Density vs Path Length (size ~ GC size)",
                     clust_sum["Cluster"], sizes,
                     outdir/"density_vs_pathlen.png")

        save_scatter(clust_sum["algebraic_connectivity_gc"], clust_sum["avg_shortest_path_len_gc"],
                     "Algebraic connectivity (GC)", "Avg shortest path length (GC)",
                     "Connectivity vs Path Length",
                     clust_sum["Cluster"], np.full(len(clust_sum), 30),
                     outdir/"conn_vs_pathlen.png")

        d = clust_sum.sort_values("n_components", ascending=True)
        save_bar(d, "Cluster", "n_components",
                 "# Connected Components by Cluster",
                 outdir/"components_by_cluster.png")

        metrics = ["n_nodes","n_edges","density","n_components","giant_component_size",
                   "avg_clustering","avg_shortest_path_len_gc","algebraic_connectivity_gc"]
        Z = clust_sum.copy()
        for m in metrics:
            vals = pd.to_numeric(Z[m], errors='coerce')
            mu, sd = vals.mean(), vals.std(ddof=1)
            Z[m] = (vals - mu) / sd if sd and not np.isnan(sd) and sd != 0 else 0.0
        dfz = Z.set_index("Cluster")[metrics]
        has_nan = dfz.isna().any().any() or np.isinf(dfz.values).any()
        if has_nan:
            plt.figure(figsize=(12, 8))
            sns.heatmap(dfz.fillna(0), cmap="magma_r", center=0,
                        yticklabels=True, xticklabels=True,
                        cbar_kws={"label": "z-score"})
            plt.xticks(rotation=45, ha='right')
            plt.title("Cluster Metric Z-scores")
            plt.tight_layout()
            plt.savefig(outdir/"cluster_metric_heatmap.png", bbox_inches="tight", dpi=300)
            plt.close()
        else:
            g = sns.clustermap(dfz, cmap="magma_r", center=0, figsize=(12,8),
                               yticklabels=True, xticklabels=True,
                               cbar_kws={"label": "z-score"})
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
            g.fig.suptitle("Cluster Metric Z-scores", y=1.01)
            plt.savefig(outdir/"cluster_metric_heatmap.png", bbox_inches="tight", dpi=300)
            plt.close()

    # ------------------------------ RMSF option --------------------------------
    rmsf_path = Path(args.rmsf_csv)
    if rmsf_path.exists() and not graph_df.empty:
        try:
            rmsf_df = pd.read_csv(rmsf_path)
            if "PDB" not in rmsf_df.columns and "PDB_ensemble" in rmsf_df.columns:
                rmsf_df["PDB"] = rmsf_df["PDB_ensemble"]
            if "delta_RMSF" not in rmsf_df.columns:
                raise ValueError("RMSF CSV must contain 'delta_RMSF'")

            cluster_rmsf = (graph_df.merge(rmsf_df[["PDB","delta_RMSF"]], on="PDB", how="left")
                            .groupby("Cluster")["delta_RMSF"].mean()
                            .rename("avg_delta_RMSF"))
            cluster_metrics = graph_df.groupby("Cluster").agg({
                "n_nodes": "mean",
                "n_edges": "mean",
                "density": "mean",
                "n_components": "mean",
                "giant_component_size": "mean",
                "avg_clustering": "mean",
                "avg_shortest_path_len_gc": "mean",
                "algebraic_connectivity_gc": "mean"
            })
            corr_df = cluster_metrics.merge(cluster_rmsf, left_index=True, right_index=True)

            correlations = {}
            for metric in cluster_metrics.columns:
                mask = ~(corr_df[metric].isna() | corr_df["avg_delta_RMSF"].isna())
                if mask.sum() > 1:
                    r = np.corrcoef(corr_df.loc[mask, metric], corr_df.loc[mask, "avg_delta_RMSF"])[0,1]
                    correlations[metric] = r
                    plt.figure(figsize=(6,5))
                    plt.scatter(corr_df.loc[mask, metric], corr_df.loc[mask, "avg_delta_RMSF"])
                    plt.xlabel(metric); plt.ylabel("Average ΔRMSF")
                    plt.title(f"Correlation: r = {r:.3f}")
                    plt.tight_layout()
                    plt.savefig(outdir/f"rmsf_corr_{metric}.png", dpi=300, bbox_inches="tight")
                    plt.close()

            pd.Series(correlations).to_csv(outdir/"rmsf_metric_correlations.csv")
        except Exception as e:
            print(f"[RMSF] Skipped due to error: {e}")

    print("Done.")

if __name__ == "__main__":
    main()
