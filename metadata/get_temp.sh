#!/bin/bash

# Output file
PDB_file=/dors/wankowicz_lab/stephanie/RT/pdbs.txt
base_dir='/dors/wankowicz_lab/stephanie/RT/'
output_file="/dors/wankowicz_lab/stephanie/RT/RT_temperature.csv"

# Create the output file and add headers
echo "PDB,temp" > "$output_file"

cd ${base_dir}
for pdb_id in PDB_file; do
    cd ${pdb_id}
    temp=$(grep "REMARK 200  TEMPERATURE           (KELVIN) :" "$pdb_file" | awk '{print $NF}')        
    if [[ ! -z "$temp" ]]; then
        # Append the PDB ID and temperature to the output file
        echo "$pdb_id,$temp" >> "$output_file"
    fi
done