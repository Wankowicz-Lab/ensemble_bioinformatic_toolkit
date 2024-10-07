#!/bin/bash
#SBATCH --mem=4G  # Request 4GB of memory
#SBATCH --cpus-per-task=1  # Request 8 CPUs per task
#SBATCH --time=02:00:00  # Request a walltime of 2 hours
#SBATCH --cpus-per-task=1
#SBATCH --array=1-20%10  # Set up a job array with the number of PDBs you have 
#SBATCH --output=bioinformatics._stdout   # Set up with the output you want. 



#__________________SOURCE PHENIX/QFIT________________________________________________#
source /dors/wankowicz_lab/phenix-installer-dev-5366-intel-linux-2.6-x86_64-centos6/phenix-dev-5366/phenix_env.sh
export PHENIX_OVERWRITE_ALL=true

source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit


#________________PDB INFO__________________________________
PDB_file=/dors/wankowicz_lab/stephanie/Kojetin_lab/qfit_pdbs.txt
PDB_dir='/dors/wankowicz_lab/stephanie/Kojetin_lab/'
output_dir='/dors/wankowicz_lab/stephanie/Kojetin_lab/output'

category=''
apo_pdb=''
cd ${output_dir}

PDB=$(cat $PDB_file | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
echo $PDB

find_largest_ligand.py ${PDB_dir}/${PDB}/${PDB}_qFit.pdb ${PDB}
lig=$(cat "${output_dir}/${PDB}_ligand_name.txt")
echo $lig

get_lig_chain_res.py ${PDB}.pdb $lig
chain_resi=$(cat "${lig}_chain_resi.txt")
#IFS=',' read -r chain resi <<< "$chain_resi"


##get information on pocket
conda activate PE
fpocket -f ${PDB_dir}/${PDB}/${PDB}_qFit.pdb -r ${resi}:${lig}:${chain} -x

python 


#_______________________GET DISTANCE OF RESIDUES FROM LIGAND OF INTEREST___________________
source /sb/sbgrid/programs/sbgrid.shrc

pymol -c  /dors/wankowicz_lab/stephanie/script/bioinformatics/find_close_residues.py -- ${PDB_dir}/${PDB}/${PDB}_qFit.pdb ${resi} ${lig} ${chain} 5.0
pymol -c  /dors/wankowicz_lab/stephanie/script/bioinformatics/find_close_residues.py -- ${PDB_dir}/${PDB}/${PDB}_qFit.pdb ${resi} ${lig} ${chain} 10.0


python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/subset_output_apo.py ${PDB} ${apo_PDB} -dist 5.0 -qFit=N -lig ${lig}
python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/subset_output_apo.py ${PDB} ${apo_PDB} -dist 10.0 -qFit=N -lig ${lig}
