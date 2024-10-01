import Bio.PDB
import pandas as pd
from Bio.PDB import *
import argparse

# Example usage from command line:
# python OP_to_bfactor.py input.pdb data.csv output.pdb --column_name s2calc

def process_pdb(pdb_filename, df_full_filename, output_filename, column_name='OP_Diff'):
    structure = Bio.PDB.PDBParser().get_structure('structure', pdb_filename)
    df_full = pd.read_csv(df_full_filename, index_col=None, sep=',', header=0)

    for model in structure:
        for chain in model:
            df = df_full[df_full['chain'] == chain.get_id()]
            df['resi'] = df['resi'].astype(int)
            for residue in chain:
                if residue.get_full_id()[3][0] != " ":
                    b = 0
                    
                if column_name.isin(['s2calc', 'OP_Diff']):
                    if residue.get_resname() in ['PRO', 'GLY']:
                        b = 0
                else:
                    b = df[df['resi'] == int(residue.get_full_id()[3][1])][column_name].values
                    b = b[0] if len(b) > 0 else 0

                for atom in residue.get_unpacked_list():
                    atom.set_bfactor(b)

    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(output_filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB file and set B-factors based on CSV data.")
    parser.add_argument("pdb_filename", help="Input PDB file")
    parser.add_argument("df_full_filename", help="Input CSV file with data")
    parser.add_argument("output_filename", help="Output PDB file")
    parser.add_argument("--column_name", default='OP_Diff',
                        help="Column name to use for B-factor (default: OP_Diff)")

    args = parser.parse_args()

    process_pdb(args.pdb_filename, args.df_full_filename, args.output_filename, args.column_name)

