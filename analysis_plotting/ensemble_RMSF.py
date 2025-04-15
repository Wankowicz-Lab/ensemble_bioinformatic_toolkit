import os
import glob
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, is_aa
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_rmsf(residue_coords):
    mean_coords = np.mean(residue_coords, axis=0)
    rmsf = np.sqrt(np.mean(np.sum((residue_coords - mean_coords) ** 2, axis=1)))
    return rmsf

pdb_files = glob.glob('/Users/stephaniewanko/Downloads/vanderbilt/mac1/oP/ensemble/*.pdb')
parser = PDBParser(QUIET=True)


all_rmsf_dfs = []  # List to store all RMSF DataFrames

for pdb_file in pdb_files:
    residue_rmsf = {}
    structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)
    pdb_base_name = os.path.splitext(os.path.basename(pdb_file))[0]  # Extract base name without extension
    print(pdb_base_name)
    for model in structure:
        for chain in model:
            if chain.id != 'A':
                continue
            for residue in chain:
                if is_aa(residue):
                    res_id = (chain.id, residue.id[1])
                    if res_id not in residue_rmsf:
                        residue_rmsf[res_id] = []
                    atom_coords = [atom.coord for atom in residue]
                    residue_rmsf[res_id].append(np.mean(atom_coords, axis=0))

    rmsf_values = {res_id: calculate_rmsf(np.array(coords)) for res_id, coords in residue_rmsf.items()}

    rmsf_df = pd.DataFrame(list(rmsf_values.items()), columns=['Residue', 'RMSF'])
    rmsf_df[['Chain', 'Resi']] = pd.DataFrame(rmsf_df['Residue'].tolist(), index=rmsf_df.index)
    rmsf_df = rmsf_df.drop(columns=['Residue'])
    rmsf_df['PDB'] = pdb_base_name  # Add a column with pdb_base_name
    all_rmsf_dfs.append(rmsf_df)  # Append the DataFrame to the list
    rmsf_df.to_csv(f'rmsf_output_{pdb_base_name}.csv', index=False)

    residues = [str(res) for res in rmsf_df['Resi']]
    rmsf = rmsf_df['RMSF']

    plt.figure(figsize=(10, 6))
    plt.plot(residues, rmsf, marker='o')
    plt.xlabel('Residue')
    plt.ylabel('RMSF')
    plt.title('RMSF per Residue')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f'rmsf_plot_{pdb_base_name}.png')

# Concatenate all RMSF DataFrames into a single DataFrame
concatenated_rmsf_df = pd.concat(all_rmsf_dfs, ignore_index=True)
concatenated_rmsf_df.to_csv('concatenated_rmsf_output.csv', index=False)

# Plot all RMSF on top of each other with different colors per PDB
plt.figure(figsize=(12, 8))
for pdb in concatenated_rmsf_df['PDB'].unique():
    pdb_rmsf_df = concatenated_rmsf_df[concatenated_rmsf_df['PDB'] == pdb]
    plt.plot(pdb_rmsf_df['resi'], pdb_rmsf_df['RMSF'], label=pdb)

plt.xlabel('Residue')
plt.ylabel('RMSF')
plt.title('RMSF per Residue for All PDBs')
plt.legend(title='PDB Structures', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('all_rmsf_overlay.png')
plt.close()
