#!/usr/bin/env python3
"""
SASA comparison with per-PDB + per-cluster summaries.

What it adds vs previous:
- NEW: Per-PDB summary table with cluster attached (which PDBs become most exposed/buried).

Inputs/Assumptions
------------------
- FreeSASA outputs in PDB format (per-atom SASA in B-factor).
- Cluster CSV columns: PDB,ID,Affinity,Cluster  (ID like 'x3430')
- Reference structure appears in filename (default: 7KQO).

Outputs
-------
- residue_sasa_changes_overall.csv         (per-residue vs ALL)
- residue_sasa_changes_by_cluster.csv      (per-residue within cluster)
- cluster_sasa_summary.csv                 (cluster-level exposure ranking)
- pdb_sasa_summary.csv                     (NEW: per-PDB exposure summary + cluster)

Usage
-----
python sasa_compare_with_cluster_summary.py \
  --sasa_dir results_freesasa \
  --clusters_csv HBDScan_clusters_3ormore_normpdb.csv \
  --ref 7KQO \
  --outdir sasa_out \
  --topk 20
"""
import argparse
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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

# ---------- tiny helpers ----------
def is_atom(line: str) -> bool:
    return line.startswith("ATOM  ") or line.startswith("HETATM")

def ffloat(s: str) -> float:
    try:
        return float(s.strip())
    except Exception:
        return 0.0

def elem_from_fields(atom_name: str, element_field: str) -> str:
    e = element_field.strip().upper()
    if e:
        return e
    m = re.match(r"\s*([A-Za-z])", atom_name)
    return m.group(1).upper() if m else "X"

def is_polar(elem: str) -> bool:
    # Simple rule: C = non-polar; N/O/S/P = polar; others -> polar.
    return False if elem == "C" else True

def res_key(chain: str, resi: int, icode: str, resn: str):
    return (chain or "_", int(resi), (icode or " ").strip() or " ", resn.strip())

def res_label(key) -> str:
    chain, resi, icode, resn = key
    icode = "" if icode == " " else icode
    return f"{resn}{resi}{icode}:{chain}"

def extract_id_from_name(name: str):
    """
    Prefer your 'x####' IDs; else accept classic 4-char PDB codes.
    """
    m = re.search(r'(^|[^A-Za-z0-9])(x\d{4})([^A-Za-z0-9]|$)', name, flags=re.I)
    if m: return m.group(2).lower()
    m = re.search(r'(^|[^A-Za-z0-9])([0-9][A-Za-z0-9]{3})([^A-Za-z0-9]|$)', name)
    if m: return m.group(2).upper()
    return None

# ---------- parsing ----------
def parse_residue_sasa_from_pdb(pdb_path: Path) -> dict:
    """
    Sum per-atom SASA (from B-factor) to residue totals, also polar/nonpolar.
    Returns { (chain, resi, icode, resn): {total, polar, nonpolar} }
    """
    out = defaultdict(lambda: {"total":0.0, "polar":0.0, "nonpolar":0.0})
    with pdb_path.open("r", errors="ignore") as fh:
        for line in fh:
            if not is_atom(line): continue
            atom  = line[12:16]
            resn  = line[17:20]
            chain = line[21].strip() or "_"
            resi  = ffloat(line[22:26])
            icode = (line[26].strip() or " ")
            bfac  = ffloat(line[60:66])
            elem  = elem_from_fields(atom, line[76:78] if len(line) >= 78 else "")
            k = res_key(chain, int(resi), icode, resn)
            out[k]["total"] += bfac
            if is_polar(elem): out[k]["polar"] += bfac
            else:              out[k]["nonpolar"] += bfac
    return out

# ---------- core ----------
def load_cluster_map(csv_path: Path) -> dict:
    """
    Expect columns: PDB,ID,Affinity,Cluster
    Map: {ID -> Cluster}  (ID like 'x3430')
    """
    df = pd.read_csv(csv_path)
    df["ID"] = df["ID"].astype(str).str.strip().str.lower()
    return {row.ID: str(row.Cluster) for _, row in df.iterrows()}

def collect_files(sasa_dir: Path) -> list:
    return sorted([p for p in sasa_dir.iterdir()
                   if p.is_file() and p.suffix.lower() in {".pdb",".ent",".txt"}])

def find_ref_file(files: list, ref_id: str) -> Path:
    rid = ref_id.lower()
    for f in files:
        if rid in f.name.lower():
            return f
    raise FileNotFoundError(f"Reference '{ref_id}' not found in: {[f.name for f in files][:10]} ...")

def compare_vs_ref(files: list, ref_file: Path, id2cluster: dict):
    """
    Returns:
      overall_df   : per-residue mean/|mean| deltas across ALL others
      byclu_df     : per-residue mean/|mean| deltas within each cluster
      pdb_summary  : NEW per-PDB summary (one row per PDB) + cluster
    """
    ref = parse_residue_sasa_from_pdb(ref_file)

    # per-residue aggregations across all PDBs and by cluster
    overall = defaultdict(lambda: {"d_tot":[], "d_pol":[], "d_npl":[]})
    byclu   = defaultdict(lambda: defaultdict(lambda: {"d_tot":[], "d_pol":[], "d_npl":[]}))

    # NEW: per-PDB aggregation structures
    # sums across residues for each PDB (to characterize exposure/burial of that PDB vs ref)
    pdb_accum = {}  # pid -> dict of running sums/lists

    for f in files:
        if f == ref_file: continue
        pid = extract_id_from_name(f.name)
        if pid is None: pid = f.stem
        pid_key = pid.lower()

        cur = parse_residue_sasa_from_pdb(f)

        # Initialize container for this PDB
        entry = pdb_accum.setdefault(pid, {
            "pdb_id": pid,
            "cluster": id2cluster.get(pid_key, None),
            "residue_count": 0,
            # totals (signed)
            "sum_delta_total": 0.0,
            "sum_delta_polar": 0.0,
            "sum_delta_nonpolar": 0.0,
            # absolute/positive/negative components
            "sum_abs_delta_total": 0.0,
            "sum_positive_total": 0.0,
            "sum_negative_total": 0.0,  # (negative number or magnitude? We'll store magnitude as positive)
            "sum_positive_polar": 0.0,
            "sum_positive_nonpolar": 0.0,
        })

        for k, refv in ref.items():
            if k not in cur: continue
            d_tot = cur[k]["total"]    - refv["total"]
            d_pol = cur[k]["polar"]    - refv["polar"]
            d_npl = cur[k]["nonpolar"] - refv["nonpolar"]

            # Per-residue overall/cluster collections
            overall[k]["d_tot"].append(d_tot)
            overall[k]["d_pol"].append(d_pol)
            overall[k]["d_npl"].append(d_npl)
            clu = id2cluster.get(pid_key)
            if clu is not None:
                byclu[clu][k]["d_tot"].append(d_tot)
                byclu[clu][k]["d_pol"].append(d_pol)
                byclu[clu][k]["d_npl"].append(d_npl)

            # ---- Per-PDB aggregations (across residues) ----
            entry["residue_count"] += 1
            entry["sum_delta_total"]    += d_tot
            entry["sum_delta_polar"]    += d_pol
            entry["sum_delta_nonpolar"] += d_npl
            entry["sum_abs_delta_total"] += abs(d_tot)
            if d_tot > 0:
                entry["sum_positive_total"] += d_tot
            elif d_tot < 0:
                entry["sum_negative_total"] += -d_tot  # store magnitude
            if d_pol > 0:
                entry["sum_positive_polar"] += d_pol
            if d_npl > 0:
                entry["sum_positive_nonpolar"] += d_npl

    # ---- Build per-residue overall DF ----
    rows = []
    for k, d in overall.items():
        tots = d["d_tot"]; pols = d["d_pol"]; npls = d["d_npl"]
        rows.append({
            "residue": res_label(k),
            "chain": k[0], "resi": k[1], "icode": "" if k[2]==" " else k[2], "resname": k[3],
            "mean_delta_total": sum(tots)/len(tots) if tots else 0.0,
            "mean_abs_delta_total": sum(abs(x) for x in tots)/len(tots) if tots else 0.0,
            "mean_delta_polar": sum(pols)/len(pols) if pols else 0.0,
            "mean_delta_nonpolar": sum(npls)/len(npls) if npls else 0.0,
            "n": len(tots),
        })
    overall_df = pd.DataFrame(rows).sort_values("mean_abs_delta_total", ascending=False)

    # ---- Build per-residue by-cluster DF ----
    rows = []
    for clu, resmap in byclu.items():
        for k, d in resmap.items():
            tots = d["d_tot"]; pols = d["d_pol"]; npls = d["d_npl"]
            rows.append({
                "cluster": str(clu),
                "residue": res_label(k),
                "chain": k[0], "resi": k[1], "icode": "" if k[2]==" " else k[2], "resname": k[3],
                "mean_delta_total": sum(tots)/len(tots) if tots else 0.0,
                "mean_abs_delta_total": sum(abs(x) for x in tots)/len(tots) if tots else 0.0,
                "mean_delta_polar": sum(pols)/len(pols) if pols else 0.0,
                "mean_delta_nonpolar": sum(npls)/len(npls) if npls else 0.0,
                "n_in_cluster": len(tots),
            })
    byclu_df = pd.DataFrame(rows).sort_values(["cluster","mean_abs_delta_total"], ascending=[True, False])

    # ---- Build NEW per-PDB summary DF ----
    pdb_rows = []
    for pid, e in pdb_accum.items():
        n = max(e["residue_count"], 1)
        pdb_rows.append({
            "pdb_id": e["pdb_id"],
            "cluster": e["cluster"],
            "residues_compared": e["residue_count"],
            # signed exposure/burial across residues
            "net_delta_total": e["sum_delta_total"],
            "net_delta_polar": e["sum_delta_polar"],
            "net_delta_nonpolar": e["sum_delta_nonpolar"],
            # positive/negative magnitudes (exposure vs burial components)
            "net_positive_total": e["sum_positive_total"],
            "net_negative_total": e["sum_negative_total"],
            "net_positive_polar": e["sum_positive_polar"],
            "net_positive_nonpolar": e["sum_positive_nonpolar"],
            # means (per-residue averages for that PDB)
            "mean_delta_total": e["sum_delta_total"]/n,
            "mean_abs_delta_total": e["sum_abs_delta_total"]/n,
        })
    pdb_summary = pd.DataFrame(pdb_rows)
    if not pdb_summary.empty:
        # Rank by exposure (net_positive_total) descending
        pdb_summary = pdb_summary.sort_values(["net_positive_total","mean_delta_total"], ascending=[False, False])

    return overall_df, byclu_df, pdb_summary

def summarize_clusters(byclu_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate per-residue cluster deltas to cluster-level exposure scores.
    Positive deltas -> increased solvent exposure upon ligand binding.
    """
    if byclu_df.empty:
        return pd.DataFrame(columns=[
            "cluster","n_residues",
            "mean_delta_total","net_delta_total",
            "mean_positive_total","net_positive_total",
            "mean_delta_polar","net_delta_polar",
            "mean_positive_polar","net_positive_polar",
            "mean_delta_nonpolar","net_delta_nonpolar",
            "mean_positive_nonpolar","net_positive_nonpolar",
        ])
    def pos(x): return x.clip(lower=0.0)
    g = byclu_df.groupby("cluster", as_index=False)
    summary = g.apply(lambda df: pd.Series({
        "n_residues": len(df),
        # total
        "mean_delta_total": df["mean_delta_total"].mean(),
        "net_delta_total": df["mean_delta_total"].sum(),
        "mean_positive_total": pos(df["mean_delta_total"]).mean(),
        "net_positive_total": pos(df["mean_delta_total"]).sum(),
        # polar
        "mean_delta_polar": df["mean_delta_polar"].mean(),
        "net_delta_polar": df["mean_delta_polar"].sum(),
        "mean_positive_polar": pos(df["mean_delta_polar"]).mean(),
        "net_positive_polar": pos(df["mean_delta_polar"]).sum(),
        # nonpolar
        "mean_delta_nonpolar": df["mean_delta_nonpolar"].mean(),
        "net_delta_nonpolar": df["mean_delta_nonpolar"].sum(),
        "mean_positive_nonpolar": pos(df["mean_delta_nonpolar"]).mean(),
        "net_positive_nonpolar": pos(df["mean_delta_nonpolar"]).sum(),
    })).reset_index(drop=True)
    return summary.sort_values("net_positive_total", ascending=False)

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sasa_dir", required=True, help="Folder with FreeSASA PDB-formatted outputs")
    ap.add_argument("--clusters_csv", required=True, help="CSV with columns: PDB,ID,Affinity,Cluster")
    ap.add_argument("--ref", default="7KQO", help="Reference ID present in a filename (e.g., 7KQO or x3430)")
    ap.add_argument("--outdir", default="sasa_out")
    ap.add_argument("--topk", type=int, default=20)
    args = ap.parse_args()

    sasa_dir = Path(args.sasa_dir); outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    files = collect_files(sasa_dir)
    if not files: raise SystemExit(f"No PDB-like files in {sasa_dir}")
    ref_file = find_ref_file(files, args.ref)
    id2cluster = load_cluster_map(Path(args.clusters_csv))

    print(f"[info] Reference file: {ref_file.name}")
    print(f"[info] Loaded {len(id2cluster)} cluster assignments from {args.clusters_csv}")

    overall_df, byclu_df, pdb_summary = compare_vs_ref(files, ref_file, id2cluster)

    # Save per-residue tables
    overall_csv = outdir/"residue_sasa_changes_overall.csv"
    byclu_csv   = outdir/"residue_sasa_changes_by_cluster.csv"
    overall_df.to_csv(overall_csv, index=False)
    byclu_df.to_csv(byclu_csv, index=False)

    # Cluster and PDB summaries
    cluster_summary = summarize_clusters(byclu_df)
    cluster_csv = outdir/"cluster_sasa_summary.csv"
    cluster_summary.to_csv(cluster_csv, index=False)

    pdb_csv = outdir/"pdb_sasa_summary.csv"
    pdb_summary.to_csv(pdb_csv, index=False)

    # ---- quick console summaries ----
    k = min(args.topk, len(overall_df))
    print("\n=== Top residues by mean |ΔSASA_total| (overall) ===")
    print(overall_df[["residue","resname","chain","resi","mean_delta_total","mean_abs_delta_total","mean_delta_polar","mean_delta_nonpolar","n"]].head(k).to_string(index=False))

    if not byclu_df.empty:
        print("\n=== Top residues by mean |ΔSASA_total| per cluster ===")
        for clu, sub in byclu_df.groupby("cluster"):
            kk = min(args.topk, len(sub))
            print(f"\n[Cluster {clu}]")
            print(sub[["residue","resname","chain","resi","mean_delta_total","mean_abs_delta_total","mean_delta_polar","mean_delta_nonpolar","n_in_cluster"]].head(kk).to_string(index=False))

    if not cluster_summary.empty:
        print("\n=== Cluster ranking by increased solvent exposure (net_positive_total) ===")
        print(cluster_summary[["cluster","n_residues","net_positive_total","net_positive_polar","net_positive_nonpolar","mean_delta_total"]].to_string(index=False))
        print(f"[info] Saved cluster summary to: {cluster_csv}")

    if not pdb_summary.empty:
        print("\n=== PDBs ranked by increased solvent exposure (net_positive_total) ===")
        cols = ["pdb_id","cluster","residues_compared","net_positive_total","net_positive_polar","net_positive_nonpolar","net_delta_total","mean_delta_total"]
        print(pdb_summary[cols].head(args.topk).to_string(index=False))
        print(f"[info] Saved per-PDB summary to: {pdb_csv}")
    else:
        print("\n[warn] No per-PDB summary produced.")

    print(f"\n[done] Wrote:\n - {overall_csv}\n - {byclu_csv}\n - {cluster_csv}\n - {pdb_csv}")

    if "pdb_summary" in locals() and not pdb_summary.empty:
        # Compute cluster order by mean ΔSASA
        cluster_means = pdb_summary.groupby("cluster")["mean_delta_total"].mean().sort_values()
        ordered_clusters = cluster_means.index.tolist()


        # --- Create boxplot ---
        plt.figure(figsize=(12, 6))
        ax = sns.boxplot(
            data=pdb_summary,
            x="cluster",
            y="mean_delta_total",
            palette="magma",
            order=ordered_clusters,
            width=0.6,
            fliersize=2
        )
        ax.yaxis.grid(True, linestyle='-', which='major', color='gray', alpha=0.5)

        plt.xlabel("Cluster")
        plt.ylabel("Mean ΔSASA per PDB (Å²)")
        plt.xticks(rotation=45)
        plt.tight_layout()

        out_path = outdir / "cluster_mean_SASA_delta_total_boxplot.png"
        plt.savefig(out_path, dpi=300)



if __name__ == "__main__":
    main()
