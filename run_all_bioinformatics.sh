#!/bin/bash
#SBATCH --mem=4G  # Request 16GB of memory
#SBATCH --cpus-per-task=1  # Request 8 CPUs per task
#SBATCH --time=01:00:00  # Request a walltime of 24 hours
#SBATCH --cpus-per-task=1
#SBATCH --array=1-20%10  # Set up a job array
#SBATCH --output=qFit_testset._stdout

#__________________SOURCE PHENIX/QFIT________________________________________________#
source /dors/wankowicz_lab/phenix-installer-dev-5366-intel-linux-2.6-x86_64-centos6/phenix-dev-5366/phenix_env.sh
export PHENIX_OVERWRITE_ALL=true

source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit


#________________PDB INFO__________________________________

PDB_file=/dors/wankowicz_lab/stephanie/Kojetin_lab/qfit_pdbs.txt ## FIX AS NEEDED
PDB_dir='/dors/wankowicz_lab/stephanie/Kojetin_lab/'  ## FIX AS NEEDED
output_dir='/dors/wankowicz_lab/stephanie/Kojetin_lab/output' ## FIX AS NEEDED

cd ${output_dir}

PDB=$(cat $PDB_file | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
echo $PDB

#cp ${PDB_dir}/${PDB}/${category}/${PDB}.updated_refine_001_qFit.pdb ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb
#cp ${PDB_dir}/${PDB}/${category}/${PDB}.updated_refine_001_qFit.log ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.log

#______________________PDB STATS_______________________

##NUM ALTLOC
b_fac=$(b_factor.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb --pdb=${PDB})

##ROTAMER
phenix.rotalyze model=${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb outliers_only=False > ${output_dir}/${PDB}_rotamer_output.txt

## RFREE
python /dors/wankowicz_lab/stephanie/script/single_parse_log.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.log  ${PDB}

##VALIDATION SCRIPTS
mmtbx.validation_summary ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb > ${PDB}_validation.txt


##HBOND
python /dors/wankowicz_lab/stephanie/script/bioinformatics/calc_hbond_alt.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb
python /dors/wankowicz_lab/stephanie/script/bioinformatics/plot_hbond_contactmap.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb ${PDB}_qFit_hbonds.csv


##DIHEDRAL
python /dors/wankowicz_lab/stephanie/bioinformatics/calc_dihedral.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb

#phenix.reduce -NOFLIP ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb > ${PDB}_qFit_H.pdb
make_methyl_df.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb --pdb ${PDB}
calc_OP.py ${output_dir}/${PDB}.dat ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb ${output_dir}/${PDB}_OP.out -r 1.5 -b $b_fac
python /dors/wankowicz_lab/stephanie/script/OP_to_bfactor.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb ${output_dir}/${PDB}_OP.out ${output_dir}/${PDB}_OP.pdb --column_name s2calc



#__________________ACTIVATE ANOTHER ENVIORNMENT________________________________________________#
conda activate PE
##SASA
python /dors/wankowicz_lab/stephanie/script/bioinformatics/calc_sasa.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb
python /dors/wankowicz_lab/stephanie/script/bioinformatics/plot_sasa.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb



