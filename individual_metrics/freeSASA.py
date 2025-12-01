#!/usr/bin/env python3
"""
Parse all FreeSASA JSON outputs into CSV format with per-residue statistics.

Usage:
    python parse_freesasa_json.py
"""

import json
import csv
import os
from glob import glob

# === CONFIG ===
JSON_DIR = "/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/freesasa_per_atom"
CSV_DIR = "/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/freesasa_per_atom_csv"


def parse_freesasa_json(json_file, output_csv):
    """
    Parse FreeSASA JSON output and create CSV with per-residue data.
    """
    # Load JSON
    with open(json_file, "r") as f:
        data = json.load(f)

    # Extract PDB ID from filename
    pdb_id = os.path.basename(json_file).replace("_sasa_atom.json", "")

    rows = []

    # Navigate through the JSON structure
    for result in data.get("results", []):
        for structure in result.get("structure", []):
            for chain_data in structure.get("chains", []):
                chain = chain_data.get("label", "")

                for residue in chain_data.get("residues", []):
                    resn = residue.get("name", "")
                    resi = residue.get("number", "")

                    area = residue.get("area", {})
                    total_sasa = area.get("total", 0.0)
                    polar_sasa = area.get("polar", 0.0)
                    apolar_sasa = area.get("apolar", 0.0)

                    atoms = residue.get("atoms", [])
                    n_atoms_total = len(atoms)
                    n_atoms_polar = sum(1 for atom in atoms if atom.get("is-polar", False))
                    n_atoms_apolar = n_atoms_total - n_atoms_polar

                    total_sasa_normalized = (
                        total_sasa / n_atoms_total if n_atoms_total > 0 else 0.0
                    )
                    apolar_sasa_normalized = (
                        apolar_sasa / n_atoms_apolar if n_atoms_apolar > 0 else 0.0
                    )
                    polar_sasa_normalized = (
                        polar_sasa / n_atoms_polar if n_atoms_polar > 0 else 0.0
                    )

                    is_exposed = 1 if total_sasa > 0 else 0

                    rows.append(
                        {
                            "pdb_id": pdb_id,
                            "chain": chain,
                            "resi": resi,
                            "resn": resn,
                            "total_sasa": total_sasa,
                            "apolar_sasa": apolar_sasa,
                            "polar_sasa": polar_sasa,
                            "n_atoms_total": n_atoms_total,
                            "n_atoms_apolar": n_atoms_apolar,
                            "n_atoms_polar": n_atoms_polar,
                            "total_sasa_normalized": total_sasa_normalized,
                            "apolar_sasa_normalized": apolar_sasa_normalized,
                            "polar_sasa_normalized": polar_sasa_normalized,
                            "is_exposed": is_exposed,
                        }
                    )

    if not rows:
        return False, 0

    fieldnames = [
        "pdb_id",
        "chain",
        "resi",
        "resn",
        "total_sasa",
        "apolar_sasa",
        "polar_sasa",
        "n_atoms_total",
        "n_atoms_apolar",
        "n_atoms_polar",
        "total_sasa_normalized",
        "apolar_sasa_normalized",
        "polar_sasa_normalized",
        "is_exposed",
    ]

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return True, len(rows)


def main():
    """Process all JSON files in the input directory."""
    os.makedirs(CSV_DIR, exist_ok=True)

    json_files = sorted(glob(os.path.join(JSON_DIR, "*_sasa_atom.json")))
    if not json_files:
        print(f"No JSON files found in {JSON_DIR}")
        raise SystemExit(1)

    success_count = 0
    fail_count = 0
    skip_count = 0

    for json_file in json_files:
        pdb_id = os.path.basename(json_file).replace("_sasa_atom.json", "")
        output_csv = os.path.join(CSV_DIR, f"{pdb_id}_freesasa.csv")

        # Skip if output already exists
        if os.path.exists(output_csv):
            skip_count += 1
            continue
        ok, nrows = parse_freesasa_json(json_file, output_csv)
   


if __name__ == "__main__":
    main()
