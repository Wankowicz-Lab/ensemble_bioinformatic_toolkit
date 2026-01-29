#!/usr/bin/env python3
"""
Protein interaction analysis script.

Calculates various types of residue-residue interactions from a PDB file
and outputs results to CSV files.

"""
from __future__ import annotations
import argparse
import math
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any
from collections import namedtuple
from scipy.spatial import cKDTree

import numpy as np
import pandas as pd
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile


## INTERACTION DEFINITIONS

ACIDIC_RESIDUES = {'ASP', 'GLU'}
BASIC_RESIDUES = {'LYS', 'ARG', 'HIS'}
AROMATIC_RESIDUES = {'PHE', 'TYR', 'TRP'}
CATIONIC_RESIDUES = {'LYS', 'ARG'}  

SALT_BRIDGE_ATOMS = {
    'ASP': ['OD1', 'OD2'],
    'GLU': ['OE1', 'OE2'],
    'LYS': ['NZ'],
    'ARG': ['NH1', 'NH2', 'NE'],
    'HIS': ['ND1', 'NE2'],
}

AROMATIC_RING_ATOMS = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
}

VDW_RADII = {
    'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'H': 1.20,
}

# Hydrogen bond parameters
DA_MAX = 3.5           # Å donor–acceptor distance
H_A_MAX = 2.6          # Å H–acceptor distance  
ANGLE_MIN = 120.0      # degrees minimum angle


DonorSite = namedtuple(
    "DonorSite",
    "res_key resname chain resi atom_name altloc "
    "base_atom base_altloc coord_base coord_donor "
    "has_explicit_H H_name H_altloc H_coord "
    "is_backbone"
)

AcceptorSite = namedtuple(
    "AcceptorSite",
    "res_key resname chain resi atom_name altloc coord is_backbone",
)

def get_residue_atoms(array: struc.AtomArray, chain: str, resi: int, atom_names: list) -> np.ndarray:
    """Get coordinates of specific atoms in a residue."""
    mask = (array.chain_id == chain) & (array.res_id == resi) & np.isin(array.atom_name, atom_names)
    return array.coord[mask]


def get_ring_center(array: struc.AtomArray, chain: str, resi: int, resname: str) -> np.ndarray | None:
    """Get the center of an aromatic ring."""
    if resname not in AROMATIC_RING_ATOMS:
        return None
    coords = get_residue_atoms(array, chain, resi, AROMATIC_RING_ATOMS[resname])
    if len(coords) < 3:
        return None
    return coords.mean(axis=0)


def get_ring_normal(array: struc.AtomArray, chain: str, resi: int, resname: str) -> np.ndarray | None:
    """Get the normal vector of an aromatic ring plane."""
    if resname not in AROMATIC_RING_ATOMS:
        return None
    coords = get_residue_atoms(array, chain, resi, AROMATIC_RING_ATOMS[resname])
    if len(coords) < 3:
        return None
    v1 = coords[1] - coords[0]
    v2 = coords[2] - coords[0]
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)
    if norm < 1e-6:
        return None
    return normal / norm

def is_heavy(atom_name: str) -> bool:
    """Return True if atom is a heavy atom (not hydrogen)."""
    return not atom_name.strip().startswith("H")


def residue_key(chain: str, resi: int) -> str:
    """Create a unique key for a residue."""
    return f"{chain}:{resi}"


def norm_alt(val: Any) -> str:
    """Normalize altloc value to string."""
    if val is None:
        return ""
    s = str(val).strip()
    return "" if s in (".", "?", " ") else s


def load_structure(path: str | Path) -> struc.AtomArray:
    """Load a protein structure from a PDB file."""
    path = Path(path)
    if path.suffix.lower() != ".pdb":
        raise ValueError(f"Only .pdb files are supported: {path}")

    pdb = PDBFile.read(str(path))
    arr = pdb.get_structure(model=1, extra_fields=["b_factor", "occupancy"])
    
    # Ensure altloc annotation exists
    if "altloc" not in arr.get_annotation_categories():
        if "altloc_id" in arr.get_annotation_categories():
            altloc_vals = np.array([
                str(a).strip() if a is not None else '' 
                for a in arr.altloc_id
            ])
            arr.set_annotation("altloc", altloc_vals)
        else:
            arr.set_annotation("altloc", np.array([''] * arr.array_length()))
    
    return arr

def _get_donor_acceptor_atoms(resname: str) -> Tuple[List, List]:
    """Get donor and acceptor atom definitions for a residue."""
    donors = []
    acceptors = []
    
    # Backbone
    donors.append(("N", "CA", None))
    acceptors.append("O")
    
    # Sidechains
    if resname == "SER":
        donors.append(("OG", "CB", None))
        acceptors.append("OG")
    elif resname == "THR":
        donors.append(("OG1", "CB", None))
        acceptors.append("OG1")
    elif resname == "TYR":
        donors.append(("OH", "CZ", None))
        acceptors.append("OH")
    elif resname == "ASN":
        acceptors.append("OD1")
        donors.append(("ND2", "CG", None))
    elif resname == "GLN":
        acceptors.append("OE1")
        donors.append(("NE2", "CD", None))
    elif resname == "ASP":
        acceptors.extend(["OD1", "OD2"])
    elif resname == "GLU":
        acceptors.extend(["OE1", "OE2"])
    elif resname == "LYS":
        donors.append(("NZ", "CE", None))
    elif resname == "ARG":
        donors.append(("NE", "CZ", None))
        donors.append(("NH1", "CZ", None))
        donors.append(("NH2", "CZ", None))
    elif resname == "HIS":
        donors.append(("ND1", "CG", None))
        donors.append(("NE2", "CD2", None))
        acceptors.extend(["ND1", "NE2"])
    elif resname == "TRP":
        donors.append(("NE1", "CD1", None))
    elif resname == "CYS":
        # SH can be both donor and acceptor
        donors.append(("SG", "CB", None))
        acceptors.append("SG")
    
    return donors, acceptors


def build_sites(array: struc.AtomArray) -> Tuple[List[DonorSite], List[AcceptorSite]]:
    """Build donor and acceptor sites from structure."""
    donors = []
    acceptors = []
    
    res_starts = struc.get_residue_starts(array)
    res_starts = np.append(res_starts, array.array_length())
    
    for i in range(len(res_starts) - 1):
        # Get residue atoms via contiguous slice
        start_idx = res_starts[i]
        end_idx = res_starts[i + 1]
        res_atoms = array[start_idx:end_idx]
        
        chain = res_atoms.chain_id[0]
        resi = res_atoms.res_id[0]
        resname = res_atoms.res_name[0]
        res_key = f"{chain}:{resi}:{resname}"
        
        # Build atom lookup
        atom_lookup = {}
        for i in range(res_atoms.array_length()):
            name = res_atoms.atom_name[i]
            alt = norm_alt(res_atoms.altloc[i])
            atom_lookup[(name, alt)] = res_atoms.coord[i]
            if alt != "":
                atom_lookup[(name, "")] = res_atoms.coord[i]
        
        donor_defs, acceptor_defs = _get_donor_acceptor_atoms(resname)
        
        # Process donors
        for donor_def in donor_defs:
            donor_name, base_name, h_name = donor_def
            altlocs = sorted({
                alt for (name, alt) in atom_lookup.keys() if name == donor_name
            })
            
            for alt in altlocs:
                if base_name and (base_name, alt) not in atom_lookup and (base_name, "") not in atom_lookup:
                    continue
                
                donor_coord = atom_lookup[(donor_name, alt)]
                base_coord = atom_lookup.get((base_name, alt))
                if base_coord is None:
                    base_coord = atom_lookup.get((base_name, ""))
                
                is_bb = donor_name == "N"
                
                donors.append(DonorSite(
                    res_key=res_key, resname=resname, chain=chain, resi=resi,
                    atom_name=donor_name, altloc=alt,
                    base_atom=base_name, base_altloc=alt,
                    coord_base=base_coord, coord_donor=donor_coord,
                    has_explicit_H=False, H_name=None, H_altloc=None, H_coord=None,
                    is_backbone=is_bb
                ))
        
        # Process acceptors
        for acc_name in acceptor_defs:
            altlocs = sorted({
                alt for (name, alt) in atom_lookup.keys() if name == acc_name
            })
            for alt in altlocs:
                is_bb = acc_name == "O"
                
                acceptors.append(AcceptorSite(
                    res_key=res_key, resname=resname, chain=chain, resi=resi,
                    atom_name=acc_name, altloc=alt,
                    coord=atom_lookup[(acc_name, alt)],
                    is_backbone=is_bb
                ))
    
    return donors, acceptors


def detect_hbonds(donors: List[DonorSite], acceptors: List[AcceptorSite]) -> List[Dict]:
    """Detect hydrogen bonds between donors and acceptors."""
    hbonds = []
    
    for d in donors:
        if d.coord_donor is None or d.coord_base is None:
            continue
            
        for a in acceptors:
            if a.coord is None:
                continue
            
            # Skip same residue
            if d.res_key == a.res_key:
                continue
            
            # Check D-A distance
            da_vec = a.coord - d.coord_donor
            da_dist = np.linalg.norm(da_vec)
            
            if da_dist > DA_MAX:
                continue
            
            # Check angle (base-donor-acceptor)
            bd_vec = d.coord_donor - d.coord_base
            bd_norm = np.linalg.norm(bd_vec)
            if bd_norm < 1e-6:
                continue
            
            cos_angle = np.dot(bd_vec, da_vec) / (bd_norm * da_dist)
            angle = math.degrees(math.acos(np.clip(cos_angle, -1, 1)))
            
            if angle < ANGLE_MIN:
                continue
            
            # Determine category
            if d.is_backbone and a.is_backbone:
                category = "backbone-backbone"
            elif d.is_backbone or a.is_backbone:
                category = "backbone-sidechain"
            else:
                category = "sidechain-sidechain"
            
            hbonds.append({
                'donor_chain': d.chain,
                'donor_resi': d.resi,
                'donor_resname': d.resname,
                'donor_atom': d.atom_name,
                'acceptor_chain': a.chain,
                'acceptor_resi': a.resi,
                'acceptor_resname': a.resname,
                'acceptor_atom': a.atom_name,
                'distance': da_dist,
                'angle': angle,
                'category': category,
            })
    
    return hbonds


# Other interaction calculations

def calculate_salt_bridges(array: struc.AtomArray, cutoff: float = 4.0) -> pd.DataFrame:
    """Identify salt bridge interactions."""
    res_starts = struc.get_residue_starts(array)
    chains = array.chain_id[res_starts]
    res_ids = array.res_id[res_starts]
    resnames = array.res_name[res_starts]
    
    results = []
    
    acidic_indices = [i for i, rn in enumerate(resnames) if rn in ACIDIC_RESIDUES]
    basic_indices = [i for i, rn in enumerate(resnames) if rn in BASIC_RESIDUES]
    
    cutoff2 = cutoff * cutoff
    
    for acid_idx in acidic_indices:
        acid_chain, acid_resi, acid_resn = chains[acid_idx], res_ids[acid_idx], resnames[acid_idx]
        acid_atoms = get_residue_atoms(array, acid_chain, acid_resi, SALT_BRIDGE_ATOMS[acid_resn])
        if len(acid_atoms) == 0:
            continue
            
        for base_idx in basic_indices:
            base_chain, base_resi, base_resn = chains[base_idx], res_ids[base_idx], resnames[base_idx]
            base_atoms = get_residue_atoms(array, base_chain, base_resi, SALT_BRIDGE_ATOMS[base_resn])
            if len(base_atoms) == 0:
                continue
            
            diff = acid_atoms[:, None, :] - base_atoms[None, :, :]
            d2 = np.einsum("ijk,ijk->ij", diff, diff)
            
            if d2.min() <= cutoff2:
                results.append({
                    'chain': acid_chain, 'resi': int(acid_resi), 'resn': acid_resn,
                    'has_salt_bridge': 'Y',
                    'partner_chain': base_chain, 'partner_resi': int(base_resi), 'partner_resn': base_resn
                })
                results.append({
                    'chain': base_chain, 'resi': int(base_resi), 'resn': base_resn,
                    'has_salt_bridge': 'Y',
                    'partner_chain': acid_chain, 'partner_resi': int(acid_resi), 'partner_resn': acid_resn
                })
    
    if not results:
        return pd.DataFrame(columns=['chain', 'resi', 'resn', 'has_salt_bridge', 
                                      'partner_chain', 'partner_resi', 'partner_resn'])
    
    return pd.DataFrame(results)


def calculate_disulfide_bonds(array: struc.AtomArray, cutoff: float = 2.5) -> pd.DataFrame:
    """Identify disulfide bonds."""
    res_starts = struc.get_residue_starts(array)
    chains = array.chain_id[res_starts]
    res_ids = array.res_id[res_starts]
    resnames = array.res_name[res_starts]
    
    results = []
    
    cys_indices = [i for i, rn in enumerate(resnames) if rn == 'CYS']
    
    cutoff2 = cutoff * cutoff
    
    for i, cys1_idx in enumerate(cys_indices):
        cys1_chain, cys1_resi = chains[cys1_idx], res_ids[cys1_idx]
        cys1_sg = get_residue_atoms(array, cys1_chain, cys1_resi, ['SG'])
        if len(cys1_sg) == 0:
            continue
            
        for cys2_idx in cys_indices[i+1:]:
            cys2_chain, cys2_resi = chains[cys2_idx], res_ids[cys2_idx]
            cys2_sg = get_residue_atoms(array, cys2_chain, cys2_resi, ['SG'])
            if len(cys2_sg) == 0:
                continue
            
            d2 = np.sum((cys1_sg[0] - cys2_sg[0])**2)
            
            if d2 <= cutoff2:
                results.append({
                    'chain': cys1_chain, 'resi': int(cys1_resi), 'resn': 'CYS',
                    'has_disulfide': 'Y',
                    'partner_chain': cys2_chain, 'partner_resi': int(cys2_resi), 'partner_resn': 'CYS'
                })
                results.append({
                    'chain': cys2_chain, 'resi': int(cys2_resi), 'resn': 'CYS',
                    'has_disulfide': 'Y',
                    'partner_chain': cys1_chain, 'partner_resi': int(cys1_resi), 'partner_resn': 'CYS'
                })
    
    if not results:
        return pd.DataFrame(columns=['chain', 'resi', 'resn', 'has_disulfide',
                                      'partner_chain', 'partner_resi', 'partner_resn'])
    
    return pd.DataFrame(results)


def calculate_pi_stacking(array: struc.AtomArray, distance_cutoff: float = 5.5, 
                          angle_cutoff: float = 30.0) -> pd.DataFrame:
    """Identify pi-stacking interactions."""
    res_starts = struc.get_residue_starts(array)
    chains = array.chain_id[res_starts]
    res_ids = array.res_id[res_starts]
    resnames = array.res_name[res_starts]
    
    results = []
    
    aromatic_data = []
    for i, (ch, ri, rn) in enumerate(zip(chains, res_ids, resnames)):
        if rn in AROMATIC_RESIDUES:
            center = get_ring_center(array, ch, ri, rn)
            normal = get_ring_normal(array, ch, ri, rn)
            if center is not None and normal is not None:
                aromatic_data.append((i, ch, ri, rn, center, normal))
    
    cutoff2 = distance_cutoff * distance_cutoff
    angle_rad = np.radians(angle_cutoff)
    
    for i, (idx1, ch1, ri1, rn1, center1, normal1) in enumerate(aromatic_data):
        for idx2, ch2, ri2, rn2, center2, normal2 in aromatic_data[i+1:]:
            d2 = np.sum((center1 - center2)**2)
            if d2 > cutoff2:
                continue
            
            dot = abs(np.dot(normal1, normal2))
            angle = np.arccos(np.clip(dot, 0, 1))
            
            is_parallel = angle < angle_rad
            is_perpendicular = abs(angle - np.pi/2) < angle_rad
            
            if is_parallel or is_perpendicular:
                geometry = 'parallel' if is_parallel else 't-shaped'
                results.append({
                    'chain': ch1, 'resi': int(ri1), 'resn': rn1,
                    'has_pi_stacking': 'Y',
                    'partner_chain': ch2, 'partner_resi': int(ri2), 'partner_resn': rn2,
                    'geometry': geometry
                })
                results.append({
                    'chain': ch2, 'resi': int(ri2), 'resn': rn2,
                    'has_pi_stacking': 'Y',
                    'partner_chain': ch1, 'partner_resi': int(ri1), 'partner_resn': rn1,
                    'geometry': geometry
                })
    
    if not results:
        return pd.DataFrame(columns=['chain', 'resi', 'resn', 'has_pi_stacking',
                                      'partner_chain', 'partner_resi', 'partner_resn', 'geometry'])
    
    return pd.DataFrame(results)


def calculate_cation_pi(array: struc.AtomArray, cutoff: float = 6.0) -> pd.DataFrame:
    """Identify cation-pi interactions."""
    res_starts = struc.get_residue_starts(array)
    chains = array.chain_id[res_starts]
    res_ids = array.res_id[res_starts]
    resnames = array.res_name[res_starts]
    
    results = []
    
    cation_atoms = {'LYS': ['NZ'], 'ARG': ['CZ']}
    cation_data = []
    for i, (ch, ri, rn) in enumerate(zip(chains, res_ids, resnames)):
        if rn in CATIONIC_RESIDUES:
            atoms = get_residue_atoms(array, ch, ri, cation_atoms.get(rn, []))
            if len(atoms) > 0:
                cation_data.append((i, ch, ri, rn, atoms[0]))
    
    aromatic_data = []
    for i, (ch, ri, rn) in enumerate(zip(chains, res_ids, resnames)):
        if rn in AROMATIC_RESIDUES:
            center = get_ring_center(array, ch, ri, rn)
            if center is not None:
                aromatic_data.append((i, ch, ri, rn, center))
    
    cutoff2 = cutoff * cutoff
    
    for cat_idx, cat_ch, cat_ri, cat_rn, cat_coord in cation_data:
        for aro_idx, aro_ch, aro_ri, aro_rn, aro_center in aromatic_data:
            d2 = np.sum((cat_coord - aro_center)**2)
            
            if d2 <= cutoff2:
                results.append({
                    'chain': cat_ch, 'resi': int(cat_ri), 'resn': cat_rn,
                    'has_cation_pi': 'Y',
                    'partner_chain': aro_ch, 'partner_resi': int(aro_ri), 'partner_resn': aro_rn,
                    'role': 'cation'
                })
                results.append({
                    'chain': aro_ch, 'resi': int(aro_ri), 'resn': aro_rn,
                    'has_cation_pi': 'Y',
                    'partner_chain': cat_ch, 'partner_resi': int(cat_ri), 'partner_resn': cat_rn,
                    'role': 'aromatic'
                })
    
    if not results:
        return pd.DataFrame(columns=['chain', 'resi', 'resn', 'has_cation_pi',
                                      'partner_chain', 'partner_resi', 'partner_resn', 'role'])
    
    return pd.DataFrame(results)


def calculate_hydrogen_bonds(array: struc.AtomArray) -> pd.DataFrame:
    """Identify hydrogen bond interactions."""
    donors, acceptors = build_sites(array)
    hbonds = detect_hbonds(donors, acceptors)
    
    results = []
    
    for h in hbonds:
        results.append({
            'chain': h['donor_chain'], 
            'resi': int(h['donor_resi']), 
            'resn': h['donor_resname'],
            'has_hbond': 'Y',
            'partner_chain': h['acceptor_chain'],
            'partner_resi': int(h['acceptor_resi']),
            'partner_resn': h['acceptor_resname'],
            'role': 'donor',
            'category': h['category']
        })
        results.append({
            'chain': h['acceptor_chain'],
            'resi': int(h['acceptor_resi']),
            'resn': h['acceptor_resname'],
            'has_hbond': 'Y',
            'partner_chain': h['donor_chain'],
            'partner_resi': int(h['donor_resi']),
            'partner_resn': h['donor_resname'],
            'role': 'acceptor',
            'category': h['category']
        })
    
    if not results:
        return pd.DataFrame(columns=['chain', 'resi', 'resn', 'has_hbond',
                                      'partner_chain', 'partner_resi', 'partner_resn',
                                      'role', 'category'])
    
    return pd.DataFrame(results)


def calculate_vdw_contacts(array: struc.AtomArray, cutoff_factor: float = 1.0) -> pd.DataFrame:
    """Identify van der Waals contacts."""
    
    
    res_starts = struc.get_residue_starts(array)
    chains = array.chain_id[res_starts]
    res_ids = array.res_id[res_starts]
    resnames = array.res_name[res_starts]
    
    heavy_mask = np.array([is_heavy(n) for n in array.atom_name], dtype=bool)
    heavy_array = array[heavy_mask]
    
    if heavy_array.array_length() == 0:
        return pd.DataFrame(columns=['chain', 'resi', 'resn', 'has_vdw_contact',
                                      'partner_chain', 'partner_resi', 'partner_resn'])
    
    atom_chains = heavy_array.chain_id
    atom_res_ids = heavy_array.res_id
    atom_res_names = heavy_array.res_name
    atom_elements = np.array([n[0] for n in heavy_array.atom_name])
    coords = heavy_array.coord
    
    radii = np.array([VDW_RADII.get(e, 1.70) for e in atom_elements])
    
    tree = cKDTree(coords)
    max_radius = max(VDW_RADII.values())
    max_cutoff = 2 * max_radius + cutoff_factor
    
    pairs = tree.query_pairs(max_cutoff)
    
    seen_pairs = set()
    results = []
    
    for i, j in pairs:
        key_i = f"{atom_chains[i]}:{atom_res_ids[i]}"
        key_j = f"{atom_chains[j]}:{atom_res_ids[j]}"
        if key_i == key_j:
            continue
        
        pair_key = tuple(sorted([key_i, key_j]))
        if pair_key in seen_pairs:
            continue
        
        dist = np.linalg.norm(coords[i] - coords[j])
        vdw_sum = radii[i] + radii[j] + cutoff_factor
        
        if dist <= vdw_sum:
            seen_pairs.add(pair_key)
            results.append({
                'chain': atom_chains[i], 'resi': int(atom_res_ids[i]), 'resn': atom_res_names[i],
                'has_vdw_contact': 'Y',
                'partner_chain': atom_chains[j], 'partner_resi': int(atom_res_ids[j]), 'partner_resn': atom_res_names[j]
            })
            results.append({
                'chain': atom_chains[j], 'resi': int(atom_res_ids[j]), 'resn': atom_res_names[j],
                'has_vdw_contact': 'Y',
                'partner_chain': atom_chains[i], 'partner_resi': int(atom_res_ids[i]), 'partner_resn': atom_res_names[i]
            })
    
    if not results:
        return pd.DataFrame(columns=['chain', 'resi', 'resn', 'has_vdw_contact',
                                      'partner_chain', 'partner_resi', 'partner_resn'])
    
    return pd.DataFrame(results)


def analyze_interactions(
    structure_path: str | Path,
    output_dir: str | Path,
    prefix: Optional[str] = None,
) -> dict[str, pd.DataFrame]:
    """
    Analyze all interactions in a protein structure and save to CSV files.

    Parameters
    ----------
    structure_path : str or Path
        Path to PDB or mmCIF file.
    output_dir : str or Path
        Directory to save output CSV files.
    prefix : str, optional
        Prefix for output filenames. Defaults to structure filename stem.

    Returns
    -------
    dict[str, pd.DataFrame]
        Dictionary mapping interaction type to DataFrame.
    """
    structure_path = Path(structure_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if prefix is None:
        prefix = structure_path.stem
    
    arr = load_structure(structure_path)
    
    aa = arr[struc.filter_amino_acids(arr)]
    
    results = {}
    
    results['salt_bridges'] = calculate_salt_bridges(aa)
    results['disulfide_bonds'] = calculate_disulfide_bonds(aa)
    results['pi_stacking'] = calculate_pi_stacking(aa)
    results['cation_pi'] = calculate_cation_pi(aa)
    results['hydrogen_bonds'] = calculate_hydrogen_bonds(arr)
    results['vdw_contacts'] = calculate_vdw_contacts(aa)
    
    pdb_name = structure_path.name
    base_cols = ["resn", "resi", "chain", "pdb"]
    for name, df in results.items():
        df = df.copy()
        df["pdb"] = pdb_name
        rest_cols = [col for col in df.columns if col not in base_cols]
        df = df[base_cols + rest_cols]
        output_path = output_dir / f"{prefix}_{name}.csv"
        df.to_csv(output_path, index=False)
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Analyze protein interactions from a PDB file."
    )
    parser.add_argument(
        "structure",
        type=str,
        help="Path to the structure file"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default=".",
        help="Output directory for CSV files. Default is current directory."
    )
    parser.add_argument(
        "--prefix", "-p",
        type=str,
        default=None,
        help="Prefix for output filenames. Defaults to structure filename."
    )
    
    args = parser.parse_args()
    
    analyze_interactions(
        structure_path=args.structure,
        output_dir=args.output_dir,
        prefix=args.prefix,
    )


if __name__ == "__main__":
    main()



