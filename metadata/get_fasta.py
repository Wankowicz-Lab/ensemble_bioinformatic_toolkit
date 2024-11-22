import warnings
import argparse
from Bio import BiopythonWarning

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

# Step 1: Read sequences from PDB sequence file
def read_pdb_sequences(seq_file):
    sequences = {}
    with open(seq_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) == 3:
                pdb_id = parts[0]
                seq = parts[1]
                one_letter_seq = ""
                for i in range(0, len(seq), 3):
                    three_letter_code = seq[i:i+3]
                    if three_letter_code in aa_map:
                        one_letter_seq += aa_map[three_letter_code]
                    else:
                        print(f"Warning: Unrecognized amino acid code '{three_letter_code}' in {pdb_id}")
                        one_letter_seq += '?'
                sequences[pdb_id] = one_letter_seq
    return sequences

# Step 2: Create FASTA file
def create_fasta(sequences, fasta_file):
    with open(fasta_file, 'w') as f:
        for pdb_id, seq in sequences.items():
            f.write(f'>{pdb_id}\n{seq}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert PDB sequences to FASTA format.")
    parser.add_argument("seq_file", help="Input sequence file from get_seq.py")
    parser.add_argument("fasta_file", help="Output FASTA file")
    args = parser.parse_args()
    
    # Read sequences from input file
    sequences = read_pdb_sequences(args.seq_file)
    
    # Create FASTA file
    create_fasta(sequences, args.fasta_file)
