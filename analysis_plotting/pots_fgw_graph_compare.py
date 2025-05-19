#!/usr/bin/env python3

import os, glob, math, argparse, warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, NeighborSearch, is_aa
from scipy.spatial.distance import cdist
import ot


# ───────────────────────── 1. Pocket parsing ──────────────────────────────
ELEMENTS = ["C", "N", "O", "S", "P",
            "F", "CL", "BR", "I",          # halogens
            "ZN", "MG", "CA", "MN", "FE",  # common metals
            "OTHER"]

def _ligand_atoms(struct):
    """Return heavy atoms of the largest non-protein hetero residue."""
    het = [r for r in struct.get_residues()
           if not is_aa(r, standard=True) and r.get_resname() not in ("HOH", "WAT")]
    if not het:
        raise ValueError("no ligand found")
    target = max(het, key=lambda r: len([a for a in r if a.element != "H"]))
    return [a for a in target if a.element != "H"]

def one_hot_element(elem: str):
    e = elem.upper()
    if e in ELEMENTS[:-1]:
        idx = ELEMENTS.index(e)
    else:
        idx = len(ELEMENTS) - 1 
    vec = np.zeros(len(ELEMENTS), dtype=float)
    vec[idx] = 1.0
    return vec

def pocket_from_pdb(pdb_file, shell=7.0, dcut=5.0):
    """
    Returns:
        C  : (n × n) geometry matrix  (Å)
        F  : (n × d) feature matrix   (one-hot element + ligand flag)
        w  : (n,)   uniform mass      (probability)
    """
    struct = PDBParser(QUIET=True).get_structure("p", pdb_file)
    lig    = _ligand_atoms(struct)

    atoms  = [a for a in struct.get_atoms() if a.element != "H"]
    ns     = NeighborSearch(atoms)

    pocket = {a for lat in lig for a in ns.search(lat.coord, shell, "A")}
    pocket = sorted(pocket, key=lambda a: a.serial_number)
    n      = len(pocket)
    if n == 0:
        raise ValueError(f"{pdb_file}: empty pocket")

    # geometry
    xyz = np.stack([a.coord for a in pocket])
    C   = cdist(xyz, xyz)                                     # Å

    # features
    F   = []
    lig_serials = {a.serial_number for a in lig}
    for a in pocket:
        # element one-hot
        feat = one_hot_element(a.element)
        # ligand flag
        feat = np.concatenate([feat, [float(a.serial_number in lig_serials)]])
        F.append(feat)
    F = np.vstack(F)

    w = np.ones(n) / n
    return C, F, w


# ───────────────────────── 2. FGW distance ────────────────────────────────
def fgw_distance2(Ca, Cb, Fa, Fb, pa, pb,
                  alpha=0.7, epsilon=5e-3, iters=100, tol=1e-9):
    """
    Squared Fused GW distance between two pockets.
    Uses POT > 0.9 fused_gromov_wasserstein; returns scalar.
    """
    M = ot.dist(Fa, Fb, metric="sqeuclidean") 
    T, log = ot.gromov.fused_gromov_wasserstein(
                M, Ca, Cb, pa, pb,
                loss_fun='square_loss',
                alpha=alpha,
                epsilon=epsilon,
                max_iter=iters,
                tol=tol,
                verbose=False, log=True)
    return log['fgw_dist']  # squared distance


# ───────────────────────── 3. CLI + pairwise computation ─────────────────
def main(args):
    pdbs = sorted(glob.glob(os.path.join(args.pdb_dir, "*.pdb")))
    if not pdbs:
        raise SystemExit(f"no PDBs found in {args.pdb_dir}")

    pockets = {}
    for p in pdbs:
        try:
            pockets[p] = pocket_from_pdb(p, args.shell, args.dcut)
        except ValueError as e:
            warnings.warn(str(e))

    n = len(pockets)
    names = list(pockets.keys())
    D2   = np.zeros((n, n), dtype=float)

    for i in range(n):
        Ci, Fi, pi = pockets[names[i]]
        for j in range(i+1, n):
            Cj, Fj, pj = pockets[names[j]]
            d2 = fgw_distance2(Ci, Cj, Fi, Fj, pi, pj,
                               alpha=args.alpha,
                               epsilon=args.epsilon,
                               iters=args.max_iter,
                               tol=args.tol)
            D2[i, j] = D2[j, i] = d2

    np.save(args.npy_out, D2)
    pd.DataFrame(D2, index=[os.path.basename(n) for n in names],
                 columns=[os.path.basename(n) for n in names]).to_csv(args.csv_out)


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=UserWarning)
    ap = argparse.ArgumentParser(description="FGW distances between binding pockets")
    ap.add_argument("pdb_dir", help="folder with .pdb files")
    ap.add_argument("--shell",   type=float, default=7.0,
                    help="pocket radius Å around ligand (default 7)")
    ap.add_argument("--dcut",    type=float, default=5.0,
                    help="edge cutoff Å for geometry matrix (default 5)")
    ap.add_argument("--alpha",   type=float, default=0.7,
                    help="FGW balance: 1=geometry only, 0=chemistry only (default 0.7)")
    ap.add_argument("--epsilon", type=float, default=5e-3,
                    help="entropic reg ε (default 5e-3)")
    ap.add_argument("--max_iter", type=int,   default=100,
                    help="max FGW iterations (default 100)")
    ap.add_argument("--tol",     type=float, default=1e-9,
                    help="convergence tol (default 1e-9)")
    ap.add_argument("--npy_out", default="fgw_dist.npy",
                    help="output .npy filename")
    ap.add_argument("--csv_out", default="fgw_dist.csv",
                    help="output .csv filename")
    args = ap.parse_args()
    main(args)

