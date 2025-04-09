#!/bin/bash
#SBATCH --job-name=ensemble_refine          # Give your job a name
#SBATCH --array=1-380               # Task array (1-32)
#SBATCH --time=48:00:00            # Runtime limit (48 hours)
#SBATCH --mem=8G                   # Memory per job (8GB)
#SBATCH --cpus-per-task=1          # Number of CPUs per task
#SBATCH --export=ALL               # Export all environment variables
#SBATCH --requeue                  # Ensure job is requeued if it fails


#________________________________________________SET PATHS________________________________________________#
source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit2
source /dors/wankowicz_lab/phenix-installer-dev-5366-intel-linux-2.6-x86_64-centos6/phenix-dev-5366/phenix_env.sh
export PHENIX_OVERWRITE_ALL=true


PDB_file=/dors/wankowicz_lab/stephanie/macrodomain/set1/all_pdbs2.txt
base_dir='/dors/wankowicz_lab/stephanie/macrodomain/set1/'

PDB=$(cat $PDB_file | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

echo ${PDB}

cd ${base_dir}/${PDB}


if [ -f "elbow.LIG.${PDB}.cif" ]; then
phenix.ensemble_refinement ${PDB}.pdb ${PDB}.mtz elbow.LIG.${PDB}.cif ptls=0.8 wxray_coupled_tbath_offset=5 tx=0.8
else
phenix.ensemble_refinement ${PDB}.pdb ${PDB}.mtz ptls=0.8 wxray_coupled_tbath_offset=5 tx=0.8
fi
