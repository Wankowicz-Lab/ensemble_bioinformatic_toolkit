#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=1-1300%200
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --time=0-10:15:00     
#SBATCH --partition=batch
#SBATCH --output=dSASA._stdout
#SBATCH --job-name=dSASA

source /sb/sbgrid/programs/sbgrid.shrc

PDB_file=/dors/wankowicz_lab/stephanie/macrodomain/set1/all_pdbs2.txt
base_dir='/dors/wankowicz_lab/stephanie/macrodomain/set1/'

PDB=$(cat $PDB_file | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
cd ${base_dir}/${PDB}
echo ${PDB}

phenix.pdbtools ${PDB}_qFit_norm_lig.pdb remove_alt_confs=True remove="water"
echo "DIFFMODE OFF" > area_input.txt
echo "OUTPUT ATOM" >> area_input.txt
/programs/x86_64-linux/ccp4/8.0/bin.capsules/areaimol XYZIN ${PDB}_qFit_norm_lig_modified.pdb XYZOUT ${PDB}_qFit_norm_lig_sasa.pdb < area_input.txt
