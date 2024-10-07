import pandas as pd
import argparse
import os

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Convert qFit info text file to CSV.')
    parser.add_argument('input_file', type=str, help='Path to the input qFit info text file.')
    parser.add_argument('ligand_name', type=str, help='Name of the ligand.')
    parser.add_argument('ligand_residue_number', type=int, help='Residue number of the ligand.')
    args = parser.parse_args()

    # Extract the PDB name from the input file
    pdb_name = os.path.basename(args.input_file)[:4]  

    # Initialize a dictionary to hold the data
    data = {
        "Metric": [],
        "Value": [],
        "PDB": []
    }

    # Read the file and parse the contents
    with open(args.input_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # Strip leading/trailing whitespace and split by ':'
            parts = line.strip().split(':')
            if len(parts) == 2:
                metric = parts[0].strip()
                value = parts[1].strip()
                data["Metric"].append(metric)
                data["Value"].append(value)
                data["PDB"].append(pdb_name)

    # Create a DataFrame from the dictionary
    df = pd.DataFrame(data)

    # Add ligand name and residue number as additional columns
    df['Ligand Name'] = args.ligand_name
    df['Ligand Residue Number'] = args.ligand_residue_number

    # Define the output CSV file path
    output_file_path = f'{pdb_name}_qFit_info.csv'

    # Save the DataFrame to a CSV file
    df.to_csv(output_file_path, index=False)


if __name__ == '__main__':
    main()
