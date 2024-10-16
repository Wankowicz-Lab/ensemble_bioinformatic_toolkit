import argparse
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

def compute_sasa(pdb_file):
    # Initialize PDB parser and structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB_structure", pdb_file)

    # Initialize Shrake-Rupley for SASA calculation
    sr = ShrakeRupley()
    sr.compute(structure, level="A")  # Compute SASA at the atom level

    # Prepare data for CSV output
    sasa_data = []

    # Iterate through all models, chains, and residues
    for model in structure:
        for chain in model:
            for residue in chain:
                # Iterate over all altlocs of the residue
                altlocs = {atom.get_altloc() for atom in residue}  # Get unique altlocs
                
                # Iterate through each atom in the residue
                for atom in residue:
                    # Skip hydrogen atoms
                    if atom.get_name().startswith("H"):
                        continue

                    base_sasa = round(atom.sasa, 2)

                    # Record the SASA for the first altloc
                    sasa_data.append({
                        'Residue': residue.get_resname(),
                        'Chain': chain.id,
                        'Residue ID': residue.get_id()[1],
                        'Atom': atom.get_name(),
                        'Alt Loc': atom.get_altloc(),
                        'SASA': base_sasa
                    })

                    # Check and record SASA for alternate conformations
                    for altloc in altlocs:
                        if altloc and atom.get_altloc() == altloc:
                            # Get SASA for the atom with the current altloc
                            sasa_value = round(atom.sasa, 2)
                            sasa_data.append({
                                'Residue': residue.get_resname(),
                                'Chain': chain.id,
                                'Residue ID': residue.get_id()[1],
                                'Atom': atom.get_name(),
                                'Alt Loc': altloc,
                                'SASA': sasa_value
                            })

    # Convert to DataFrame and save to CSV
    df = pd.DataFrame(sasa_data)
    df = df.drop_duplicates()  # Remove duplicate rows
    output_csv = pdb_file.replace('.pdb', '_sasa.csv')
    df.to_csv(output_csv, index=False)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Calculate SASA for all atoms in a PDB file, excluding hydrogen atoms.')
    parser.add_argument('pdb_file', type=str, help='Path to the PDB file')
    args = parser.parse_args()

    compute_sasa(args.pdb_file)

if __name__ == "__main__":
    main()
