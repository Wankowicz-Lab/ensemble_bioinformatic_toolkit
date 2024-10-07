import pandas as pd
import argparse
import os

# Set up argument parsing
parser = argparse.ArgumentParser(description='Convert qFit info text file to CSV.')
parser.add_argument('input_file', type=str, help='Path to the input qFit info text file.')
parser.add_argument('ligand_name', type=str, help='Name of the ligand.')
parser.add_argument('ligand_residue_number', type=int, help='Residue number of the ligand.')
args = parser.parse_args()

# Extract the PDB name from the input file
pdb_name = os.path.basename(args.input_file)[:4]
print(pdb_name)

# Initialize a dictionary to hold the data
data = {
    "PDB": [],
    "Ligand Name": [],
    "Ligand Residue Number": [],
}

# Read the file and parse the contents
with open(args.input_file, 'r') as file:
    lines = file.readlines()
    for line in lines:
        # Strip leading/trailing whitespace and split by ':'
        parts = line.strip().split(':')
        if len(parts) == 2:
            metric = parts[0].strip().replace(' ', '_')  # Replace spaces with underscores
            value = parts[1].strip()
            data[metric] = [value]  # Create a new column for each metric with the value

# Add the PDB, ligand name, and residue number
data["PDB"].append(pdb_name)
data["Ligand Name"].append(args.ligand_name)
data["Ligand Residue Number"].append(args.ligand_residue_number)

# Create a DataFrame from the dictionary
df = pd.DataFrame(data)

# Define the output CSV file path
output_file_path = '{}_qFit_info.csv'.format(pdb_name)

# Save the DataFrame to a CSV file
df.to_csv(output_file_path, index=False)

