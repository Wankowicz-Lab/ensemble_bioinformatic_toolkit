import argparse
import freesasa

def main():
    parser = argparse.ArgumentParser(description='Process PDB file and output PDB with SASA.')
    parser.add_argument('input_pdb', type=str, help='Input PDB file')
    args = parser.parse_args()

    input_pdb = args.input_pdb
    output_pdb = input_pdb.replace('.pdb', '_sasa.pdb')

    # Load the structure from the PDB file
    structure = freesasa.Structure(input_pdb)
    result = freesasa.calc(structure)
    
    # Print total SASA
    print("Total : %.2f A2" % result.totalArea())

    # Dictionary to hold SASA for each residue, keyed by residue identifier
    residue_sasa = {}

    # Iterate through residues to calculate SASA per residue
    for i, res in enumerate(structure):
        res_id = res.id  # Unique identifier for residue (e.g., (chain_id, res_seq, icode))
        
        # Get SASA for the current residue
        sasa_value = result.residueArea(i)
        
        # If residue has alternate locations, adjust how we store the results
        if res_id not in residue_sasa:
            residue_sasa[res_id] = sasa_value
        else:
            residue_sasa[res_id] += sasa_value  # Accumulate SASA for alternate locations

    # Print residue-by-residue SASA
    for res_id, sasa in residue_sasa.items():
        print(f"Residue {res_id}: %.2f A2" % sasa)

    # Write the result to a new PDB file
    result.write_pdb(output_pdb)

if __name__ == "__main__":
    main()
