#!/usr/bin/env bash
set -o errexit          # abort on un-handled errors
set -o nounset          # treat unset vars as errors
set -o pipefail         # catch errors inside pipelines

source /programs/sbgrid.shrc
export PHENIX_OVERWRITE_ALL=true

PDB_file='/dors/wankowicz_lab/stephanie/macrodomain/set1/all_pdbs2.txt'
PDB_dir='/dors/wankowicz_lab/stephanie/macrodomain/set1'
output_dir='/dors/wankowicz_lab/stephanie/macrodomain/set1/all_output'

mkdir -p "$output_dir"
while IFS= read -r PDB || [[ -n $PDB ]]; do
    pdb_dir="${PDB_dir}/${PDB}"
    pdb="${pdb_dir}/${PDB}.pdb"
    cd $pdb_dir
    # Skip cleanly if the directory or qFit model is absent
    [[ -d $pdb_dir ]]        || { echo "⨯ dir $pdb_dir not found – skipping"; continue; }
    [[ -f $pdb_qfit ]]       || { echo "⨯ $pdb_qfit not found – skipping";  continue; }

    phenix.pdbtools "$pdb" remove_alt_confs=True remove="water"

    {
        echo "DIFFMODE OFF"
        echo "OUTPUT ATOM"
    } > area_input.txt

    if [[ -f $lig_file ]]; then
        lig=$(<"$lig_file")
        grep "$lig" "${PDB}_modified.pdb" > ligand.pdb      || true
    fi

    /programs/x86_64-linux/ccp4/8.0/bin.capsules/areaimol \
        XYZIN "${PDB}_modified.pdb" \
        XYZOUT "${PDB}_asa_nowater.pdb"      < area_input.txt    || true

done < "$PDB_file"
