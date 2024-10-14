import sys

def filter_altloc(input_pdb, output_pdb, altloc):
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                altloc_char = line[16]
                if altloc_char == ' ' or altloc_char == altloc:
                    outfile.write(line)
            else:
                outfile.write(line)

def main():
    if len(sys.argv) != 2:
        sys.exit(1)

    input_pdb = sys.argv[1]
    existing_altlocs = set()

    # First, scan the input file to find existing alternate locations
    with open(input_pdb, 'r') as infile:
        for line in infile:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                altloc_char = line[16]
                if altloc_char != ' ':
                    existing_altlocs.add(altloc_char)

    # Then, generate files only for existing alternate locations
    for altloc in existing_altlocs:
        output_pdb = f"{input_pdb.split('.')[0]}_altloc_{altloc}.pdb"
        filter_altloc(input_pdb, output_pdb, altloc)

if __name__ == "__main__":
    main()
