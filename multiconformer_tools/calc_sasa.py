import argparse
import freesasa

##THIS IS NOT WEIGHTED BY MULTICONF

def main():
    parser = argparse.ArgumentParser(description='Process PDB file and output PDB with SASA.')
    parser.add_argument('input_pdb', type=str, help='Input PDB file')
    args = parser.parse_args()

    input_pdb = args.input_pdb
    output_pdb = input_pdb.replace('.pdb', '_sasa.pdb')

    structure = freesasa.Structure(input_pdb)
    result = freesasa.calc(structure)
    result.write_pdb(output_pdb)

if __name__ == "__main__":
    main()
