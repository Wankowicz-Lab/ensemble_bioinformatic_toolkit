#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=1-69%100
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --time=0-65:15:00     
#SBATCH --output=1_refine._stdout
#SBATCH --job-name=1_refine
#
#__________________SOURCE PHENIX/QFIT________________________________________________#
source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit2

source /dors/wankowicz_lab/Phenix-2.0/phenix_env.sh 
export PHENIX_OVERWRITE_ALL=true

# Set directories and file paths
parent_dir="/dors/wankowicz_lab/PDBRedo/pdb-redo3/"
full_id_file=/dors/wankowicz_lab/stephanie/script/pdbs_to_get_done.txt
#full_id_file=/dors/wankowicz_lab/serine_protease2/Trypsin/pdbs.txt
#full_id_file=/dors/wankowicz_lab/all_pdb/70001_80000/pdb_70001_80000.txt

# Get the PDB ID for this job array task
pdb=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$full_id_file")
pdb_lower=$(echo "$pdb" | tr '[:upper:]' '[:lower:]')
two_digit_folder="${pdb_lower:1:2}"  # get 2nd and 3rd letters
pdb_path="$parent_dir/$two_digit_folder/$pdb_lower"

# Run if path exists
if [[ -d "$pdb_path" ]]; then
    echo "Running qFit for: $pdb_lower"
    cd "$pdb_path"
    if [[ ! -f "${PDB}_qFit.pdb" ]]; then
      if [[ ! -f "composite_omit_map.mtz" ]]; then
        phenix.composite_omit_map ${pdb_lower}_final.mtz ${pdb_lower}_final.pdb omit-type=refine nproc=2 r_free_flags.generate=True
      fi
      if [[ ! -f "multiconformer_model2.pdb" ]]; then
        qfit_protein composite_omit_map.mtz -l 2FOFCWT,PH2FOFCWT ${pdb_lower}_final.pdb
      fi
      qfit_final_refine_xray.sh ${pdb_lower}_final.mtz
    fi
else
    echo "Missing folder: $pdb_path"
fi

