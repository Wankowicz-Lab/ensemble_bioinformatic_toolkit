#!/usr/bin/env python3
"""
Usage:
  python residue_packing_biotite.py --pdb x1234.pdb --chain A --cutoff 5.0 \
      --out packing_x1234_A.csv

Notes:
- Hydrogens (and deuterium) are ignored.
- Waters (HOH/WAT) are ignored.
"""

import argparse
from pathlib import Path
import csv
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree

import biotite.structure as struc
import biotite.structure.io as strucio

AA3 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL"
}
WATER_NAMES = {"HOH"}


def parse_args():
    ap = argparse.ArgumentParser(description="Simple residue packing from a PDB (Biotite version).")
    ap.add_argument("--pdb", required=True, help="Input PDB file")
    ap.add_argument("--chain", default=None,
                    help="Chain ID to restrict to (default: all protein chains)")
    ap.add_argument("--cutoff", type=float, default=5.0,
                    help="Distance cutoff in Å (default: 5.0)")
    ap.add_argument("--model", type=int, default=0,
                    help="Model index to use if multiple models are present (0-based, default: 0)")
    ap.add_argument("--out", type=str, default=None,
                    help="Output CSV path (default: <pdb_stem>_packing.csv)")
    ap.add_argument("--plot", action="store_true",
                    help="Also save a quick PNG plot of packing vs residue index")
    return ap.parse_args()


def atom_is_hydrogen(element: str, atom_name: str) -> bool:
    """
    Decide if an atom is hydrogen/deuterium by element or by name.
    """
    e = (element or "").upper()
    if e in {"H", "D"}:
        return True
    name = (atom_name or "").strip().upper()
    return name.startswith("H") or name.startswith("D")


def main():
    args = parse_args()
    pdb_path = Path(args.pdb)

    if args.out:
        out_csv = Path(args.out)
    else:
        out_csv = pdb_path.with_suffix("").with_name(f"{pdb_path.stem}_packing.csv")

    # Load structure (AtomArray or AtomArrayStack)
    struct = strucio.load_structure(str(pdb_path))

    # Select model if stack
    if isinstance(struct, struc.AtomArrayStack):
        if args.model < 0 or args.model >= struct.stack_depth():
            raise IndexError(
                f"Model index {args.model} out of range (0..{struct.stack_depth()-1})"
            )
        array = struct[args.model]
    else:
        if args.model != 0:
            print("Warning: --model ignored; PDB contains only a single model.")
        array = struct

    # Convenience aliases
    chain_id = array.chain_id
    res_id   = array.res_id
    ins_code = array.ins_code
    res_name = array.res_name
    atom_name = array.atom_name
    element = array.element
    b_factor = array.b_factor
    occupancy = array.occupancy
    coords = array.coord

    n_atoms_total = array.array_length()

    # Masks for heavy atoms and non-water atoms (for neighbor search)
    is_water = np.isin(res_name, list(WATER_NAMES))
    is_heavy = np.array(
        [not atom_is_hydrogen(e, an) for e, an in zip(element, atom_name)],
        dtype=bool
    )
    ns_mask = (~is_water) & is_heavy  # atoms included in neighbor search

    if not np.any(ns_mask):
        raise RuntimeError("No heavy non-water atoms found in structure.")

    ns_indices = np.nonzero(ns_mask)[0]
    ns_coords = coords[ns_indices]

    # Build KDTree over all heavy atoms
    tree = cKDTree(ns_coords)

    # Residue boundaries
    res_starts = struc.get_residue_starts(array)
    res_ends = struc.get_residue_ends(array)

    rows = []

    for s, e in zip(res_starts, res_ends):
        # Representative atom for residue-level metadata
        r_chain = chain_id[s]
        r_resid = int(res_id[s])
        r_icode = (ins_code[s] or "").strip()
        r_resn = (res_name[s] or "").upper()

        # Chain filter
        if args.chain is not None and r_chain != args.chain:
            continue

        # Only consider standard amino-acid residues (not water)
        if r_resn in WATER_NAMES:
            continue
        if r_resn not in AA3:
            continue

        # Heavy atoms of this residue (all altlocs)
        res_indices = np.arange(s, e)
        res_heavy_mask = is_heavy[res_indices]
        if not np.any(res_heavy_mask):
            continue
        res_heavy_indices = res_indices[res_heavy_mask]

        # Map these residue atoms to the KDTree index space
        # Build a quick index map: original_index -> tree_index
        # (do this once outside the loop in a more optimized version)
        idx_map = {orig: i for i, orig in enumerate(ns_indices)}
        tree_indices = [idx_map[i] for i in res_heavy_indices if i in idx_map]
        if not tree_indices:
            continue

        # Collect neighboring residue keys (distinct) within cutoff of ANY atom
        neighbor_res_keys = set()
        for ti in tree_indices:
            center = ns_coords[ti]
            # Neighbor atoms within cutoff
            neighbor_k = tree.query_ball_point(center, r=args.cutoff)
            for nk in neighbor_k:
                orig_idx = ns_indices[nk]
                # Skip self
                if orig_idx in res_heavy_indices:
                    continue

                # Info about neighbor atom
                n_resn = (res_name[orig_idx] or "").upper()
                if n_resn in WATER_NAMES:
                    continue
                if n_resn not in AA3:
                    continue

                n_chain = chain_id[orig_idx]
                if args.chain is not None and n_chain != args.chain:
                    continue

                n_resid = int(res_id[orig_idx])
                n_icode = (ins_code[orig_idx] or "").strip()

                neighbor_res_keys.add((n_chain, n_resid, n_icode))

        # Simple stats for this residue: use all heavy atoms (all altlocs) in this residue
        res_heavy_elements = element[res_heavy_indices]
        res_heavy_atom_names = atom_name[res_heavy_indices]
        res_heavy_b = b_factor[res_heavy_indices]
        res_heavy_occ = occupancy[res_heavy_indices]

        num_atoms = len(res_heavy_indices)
        num_neighbors = len(neighbor_res_keys)
        contact_density = num_neighbors / max(1, num_atoms)

        avg_b = float(np.nanmean(res_heavy_b)) if num_atoms > 0 else np.nan
        # Default occupancy to 1.0 where missing
        occ_clean = np.where(np.isnan(res_heavy_occ), 1.0, res_heavy_occ)
        avg_occ = float(np.mean(occ_clean)) if num_atoms > 0 else np.nan

        rows.append({
            "chain": r_chain,
            "resi": r_resid,
            "icode": r_icode,
            "resn": r_resn,
            "n_atoms": num_atoms,
            "n_neighbor_residues": num_neighbors,
            "contact_density": round(contact_density, 4),
            "avg_b": round(avg_b, 2) if not np.isnan(avg_b) else "",
            "avg_occupancy": round(avg_occ, 2) if not np.isnan(avg_occ) else "",
        })

    if not rows:
        raise RuntimeError("No protein residues found with given filters.")



if __name__ == "__main__":
    main()
