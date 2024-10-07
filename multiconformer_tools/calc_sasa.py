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
    area_classes = freesasa.classifyResults(result, structure)
    print("Total : %.2f A2" % result.totalArea())
    for key in area_classes:
        print(key, ": %.2f A2" % area_classes[key])

    selections = freesasa.selectArea(('alanine, resn ala', 'r1_10, resi 1-10', 'B_alt, altloc B'),
                                 structure, result)
    for key in selections:
      print(key, ": %.2f A2" % selections[key])

    result = freesasa.calc(structure)
    result.write_pdb(output_pdb)

if __name__ == "__main__":
    main()
