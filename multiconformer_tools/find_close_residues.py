import argparse
import csv
from pymol import cmd

## THIS SCRIPT REQUIRES PYMOL. SOURCE IT USING: source /sb/sbgrid/programs/sbgrid.shrc

def get_min_distance(sel1, sel2):
    coords1 = cmd.get_coords(sel1)
    coords2 = cmd.get_coords(sel2)
    min_distance = None
    for c1 in coords1:
        for c2 in coords2:
            d = ((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2) ** 0.5
            if min_distance is None or d < min_distance:
                min_distance = d
    return min_distance

def find_close_residues(pdb_file, residue_number, residue_name, chain, distance):
    # Load the PDB file
    print(f"Loading PDB file: {pdb_file}")
    cmd.load(pdb_file)

    # Get the PDB ID from the filename (without extension)
    pdb_id = pdb_file.rsplit('.', 1)[0]

    # Define the selection for the target residue using resi and chain
    target_selection = f'chain {chain} and resi {residue_number}'

    # Select the target residue
    cmd.select('target_residue', target_selection)

    # Select all residues within the specified distance from the target residue
    cmd.select('close_residues', f'byres all within {distance} of target_residue')

    # Prepare to write the output to a CSV file
    output_file = f"{pdb_id}_{residue_name}_{distance}_closeres.csv"

    with open(output_file, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['resi', 'chain', 'PDB', 'distance'])

        # Get the list of unique residues in close_residues
        residues = set()
        model = cmd.get_model('close_residues')
        for atom in model.atom:
            resi = atom.resi
            res_chain = atom.chain
            residues.add((resi, res_chain))

        # For each residue, calculate the minimum distance to the target residue
        for resi, res_chain in residues:
            # Skip the target residue itself
            if resi == residue_number and res_chain == chain:
                continue

            # Define the selection for this residue
            residue_selection = f'chain {res_chain} and resi {resi}'

            # Calculate the minimum distance between any atom of this residue and the target residue
            distance_value = get_min_distance('target_residue', residue_selection)

            csv_writer.writerow([resi, res_chain, pdb_id, distance_value])

    # Clean up selections
    cmd.delete('target_residue')
    cmd.delete('close_residues')



parser = argparse.ArgumentParser(description='Find residues within a certain distance of a target residue in a PDB file.')
parser.add_argument('pdb_file', type=str, help='Path to the PDB file')
parser.add_argument('residue_number', type=str, help='Residue number of the target residue')
parser.add_argument('residue_name', type=str, help='Residue name of the target residue')
parser.add_argument('chain', type=str, help='Chain identifier of the target residue')
parser.add_argument('distance', type=float, help='Distance threshold in angstroms')

# Parse arguments
args = parser.parse_args()

find_close_residues(args.pdb_file, args.residue_number, args.residue_name, args.chain, args.distance)

