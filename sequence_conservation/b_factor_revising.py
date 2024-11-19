import Bio.PDB  
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import warnings
import argparse
from Bio import AlignIO, BiopythonWarning
import os


# Suppress PDB warnings
warnings.simplefilter('ignore', BiopythonWarning)

# Mapping three-letter amino acid codes to one-letter codes
aa_map = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'CSD': 'C', 'KCX': 'K', 'ALY': 'K', 'OCS': 'C',
    'TPO': 'T',
}

# Step 1: Read metadata and Resi_initial
def read_metadata_with_initial(metadata_file):
    sequences = {}
    initial_indices = {}
    with open(metadata_file, "r") as f:
        next(f)
        lines = f.readlines()
        for line in lines[1:]:
            parts = line.strip().split()
            pdb_id = parts[0]
            try:
                chain_index = parts.index('A') + 1
                chain_a_seq = parts[chain_index]
                resi_index = chain_index + 1
                resi_initial = int(parts[resi_index])
            except (ValueError, IndexError):
                print(f"Warning: Invalid metadata for {pdb_id}")
                continue
            one_letter_seq = ""
            for i in range(0, len(chain_a_seq), 3):
                three_letter_code = chain_a_seq[i:i+3]
                if three_letter_code in aa_map:
                    one_letter_seq += aa_map[three_letter_code]
                else:
                    print(f"Warning: Unrecognized amino acid code '{three_letter_code}' in {pdb_id}")
                    one_letter_seq += '?'
            sequences[pdb_id] = one_letter_seq
            initial_indices[pdb_id] = resi_initial
    return sequences, initial_indices


def create_fasta(sequences, fasta_file):
    with open(fasta_file, 'w') as f:
        for pdb_id, seq in sequences.items():
            f.write(f'>{pdb_id}\n{seq}\n')

def run_mafft(fasta_file="sequences.fasta"):
    aligned_file = "aligned.fasta"
    subprocess.run(['wsl', 'mafft', '--auto', fasta_file], stdout=open(aligned_file, 'w'))
    print(f"Alignment completed. Output saved to {aligned_file}.")
    return aligned_file

def calculate_shannon_entropy(aligned_sequences, output_file="residue_entropy_list.txt"):
    entropy_values = []
    residue_entropy_list = [] 
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY-')
    for idx, column in enumerate(zip(*aligned_sequences), start=1):
        counts = np.array([column.count(aa) for aa in amino_acids])
        frequencies = counts / np.sum(counts) if np.sum(counts) > 0 else np.zeros_like(counts)
        entropy = -np.sum(frequencies[frequencies > 0] * np.log2(frequencies[frequencies > 0]))
        max_entropy = np.log2(len(frequencies[frequencies > 0])) if len(frequencies[frequencies > 0]) > 0 else 0
        normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0
        entropy_values.append(normalized_entropy)
        residue_entropy_list.append((idx, normalized_entropy))
    with open(output_file, 'w') as f:
        f.write("Residue Number\tEntropy\n")
        for resi, entropy in residue_entropy_list:
            f.write(f"{resi}\t{entropy:.4f}\n")
    return entropy_values

# Step 2: Process PDBs with aligned_pos and actual_resi
def process_all_pdbs(folder_path, entropy_values, aligned_sequences, initial_indices, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    for filename in os.listdir(folder_path):
        if filename.endswith('.pdb'):
            pdb_filename = os.path.join(folder_path, filename)
            output_filename = os.path.join(output_folder, filename)

            structure = Bio.PDB.PDBParser().get_structure('structure', pdb_filename)
            pdb_id = filename.split('_')[0]  
            resi_initial = initial_indices[pdb_id]
            aligned_sequence = aligned_sequences[0]  
            
            aligned_pos = 0
            residue_pos = 0

            for model in structure:
                for chain in model:
                    if chain.id == 'A':  
                        for residue in chain:
                            while aligned_pos < len(aligned_sequence) and aligned_sequence[aligned_pos] == '-':
                                aligned_pos += 1

                            if aligned_pos < len(aligned_sequence):
                                current_char = aligned_sequence[aligned_pos]

                                if current_char in list('ACDEFGHIKLMNPQRSTVWY'):
                                    residue_pos += 1
                                    actual_resi = residue_pos + resi_initial - 1
                                    entropy_value = entropy_values[aligned_pos]

                                    for atom in residue:
                                        if residue.id[1] == actual_resi:
                                            atom.bfactor = entropy_value

                                aligned_pos += 1

            io = Bio.PDB.PDBIO()
            io.set_structure(structure)
            io.save(output_filename)

def plot_entropy(entropy_values):
    resi_numbers = list(range(1, len(entropy_values) + 1))
    plt.plot(resi_numbers, entropy_values, marker='o')
    plt.xlabel('Residue Number')
    plt.ylabel('Shannon Entropy')
    plt.title('Residue Entropy in Aligned Sequences')
    plt.grid()
    plt.savefig('residue_entropy_plot.png')
    plt.show()

# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB files and calculate Shannon entropy.")
    parser.add_argument("metadata_file", help="Input metadata file")
    parser.add_argument("output_folder", help="Output folder for revised PDB files")
    args = parser.parse_args()
    # Read metadata and initial residue index
    sequences, initial_indices = read_metadata_with_initial(args.metadata_file)
    # Create FASTA file for alignment
    fasta_file = 'sequences.fasta'
    create_fasta(sequences, fasta_file)
    # Run MAFFT alignment
    aligned_file = run_mafft(fasta_file)
    alignment = AlignIO.read(aligned_file, "fasta")
    # Calculate Shannon entropy for the aligned sequences
    entropy_values = calculate_shannon_entropy([str(record.seq) for record in alignment])
    # Plot the entropy values
    plot_entropy(entropy_values)
    # Process all PDB files for B-factor assignment
    process_all_pdbs("CDK_PDBs", entropy_values, alignment, initial_indices, args.output_folder)
