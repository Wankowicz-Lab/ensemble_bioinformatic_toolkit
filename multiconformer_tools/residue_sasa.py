import argparse
import os
import pandas as pd
from Bio.PDB import PDBParser, ShrakeRupley, is_aa

# Maximum SASA values for RSA calculation (Kabsch & Sander, 1983)
max_sasa = {
    "A": 113, "R": 241, "N": 158, "D": 151, "C": 140, "E": 183, "Q": 189, "G": 85,
    "H": 194, "I": 182, "L": 180, "K": 211, "M": 204, "F": 218, "P": 143, "S": 122,
    "T": 146, "W": 259, "Y": 229, "V": 160
}

def compute_rsa(pdb_file):
    """Computes Relative Solvent Accessibility (RSA) using Shrake-Rupley algorithm."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]  # Use first model

    # Compute SASA
    sr = ShrakeRupley()
    sr.compute(model, level="R")  # Compute at residue level

    rsa_data = []
    
    for chain in model:
        for residue in chain:
            if is_aa(residue):  # Ignore non-amino acid residues
                res_name = residue.get_resname()
                res_id = residue.get_id()[1]  # Residue ID number
                sasa = residue.sasa  # SASA from Shrake-Rupley

                one_letter_code = three_to_one_dict.get(res_name, "X")  # Convert to one-letter code
                
                if one_letter_code in max_sasa:
                    rsa = sasa / max_sasa[one_letter_code]  # Normalize SASA
                    rsa_data.append({
                        "Residue": one_letter_code + str(res_id),
                        "Residue ID": res_id,
                        "SASA": round(sasa, 2),
                        "RSA": round(rsa, 3),
                        "Exposure": "Solvent-Exposed" if rsa > 0.25 else "Buried"
                    })

    # Convert to DataFrame
    df = pd.DataFrame(rsa_data)

    # Save results
    output_file = os.path.splitext(pdb_file)[0] + "_rsa.csv"
    df.to_csv(output_file, index=False)

    print(f"Saved RSA results to {output_file}")

# Convert three-letter to one-letter amino acid code
three_to_one_dict = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RSA from PDB using Shrake-Rupley SASA method")
    parser.add_argument("pdb_file", help="Input PDB file")
    args = parser.parse_args()
    
    compute_rsa(args.pdb_file)
