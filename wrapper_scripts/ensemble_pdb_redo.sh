#!/bin/bash
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --time=60:00:00
#SBATCH --array=1-32%100  # Adjust 1-64 to match number of entries in abl_ids.txt
#SBATCH --output=extreme_ensemble_stdout_%A_%a.out

#__________________SOURCE PHENIX/QFIT________________________________________________#
source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit

source /dors/wankowicz_lab/phenix-installer-1.21.2-5419-intel-linux-2.6-x86_64-centos6/phenix-1.21.2-5419/phenix_env.sh
export PHENIX_OVERWRITE_ALL=true

# Set directories and file paths
#parent_dir="/dors/wankowicz_lab/PDBRedo/pdb-redo3/"
#full_id_file="/dors/wankowicz_lab/stephanie/script/all_kinase.txt"

parent_dir="/dors/wankowicz_lab/PDBRedo/pdb-redo3/"
#full_id_file=/dors/wankowicz_lab/all_pdb/1_10000/pdb_1_10000.txt
full_id_file=/dors/wankowicz_lab/serine_protease/Chymotrypsin/pdbs.txt

# Get the PDB ID for this job array task
pdb=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$full_id_file")
pdb_lower=$(echo "$pdb" | tr '[:upper:]' '[:lower:]')
two_digit_folder="${pdb_lower:1:2}"  # get 2nd and 3rd letters
pdb_path="$parent_dir/$two_digit_folder/$pdb_lower"

# Run if path exists
if [[ -d "$pdb_path" ]]; then
    echo "Running ensemble refinement for: $pdb_lower"
    cd "$pdb_path"
    if [[ -f "${pdb_lower}_final_ensemble.pdb.gz" ]]; then
       continue
    else
      phenix.elbow "${pdb_lower}_final.pdb" --do_all
      if [[ -f "elbow.${pdb_lower}_final_pdb.all.001.cif" ]]; then 
        phenix.ensemble_refinement "${pdb_lower}_final.pdb" "${pdb_lower}_final.mtz" "elbow.${pdb_lower}_final_pdb.all.001.cif" ptls=0.8 wxray_coupled_tbath_offset=5 tx=0.8
      else
        phenix.ensemble_refinement "${pdb_lower}_final.pdb" "${pdb_lower}_final.mtz" ptls=0.8 wxray_coupled_tbath_offset=5 tx=0.8
      fi
    fi
else
    echo "Missing folder: $pdb_path"
fi

