#!/bin/bash
#SBATCH --mem=4G  # Request 4GB of memory
#SBATCH --cpus-per-task=1  # Request 1 CPUs per task
#SBATCH --time=02:00:00  # Request a walltime of 2 hours
#SBATCH --cpus-per-task=1
#SBATCH --array=1-20%10  # Set up a job array with the number of PDBs you have
#SBATCH --output=bioinformatics._stdout # Set up an output file with 


#__________________SOURCE PHENIX/QFIT________________________________________________#
source /dors/wankowicz_lab/phenix-installer-dev-5366-intel-linux-2.6-x86_64-centos6/phenix-dev-5366/phenix_env.sh
export PHENIX_OVERWRITE_ALL=true

source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit


#________________PDB INFO__________________________________

PDB_file=/dors/wankowicz_lab/stephanie/apo_holo_dataset/pdbs.txt ## FIX AS NEEDED
PDB_dir='/dors/wankowicz_lab/stephanie/apo_holo_dataset/'  ## FIX AS NEEDED
output_dir='/dors/wankowicz_lab/stephanie/apo_holo_dataset/output' ## FIX AS NEEDED


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
python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/multiconformer_tools/single_parse_log.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.log  ${PDB}

##VALIDATION SCRIPTS
mmtbx.validation_summary ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb > ${PDB}_validation.txt

##HBOND
python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/multiconformer_tools/weighted_hydrogen_bond.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb

##DIHEDRAL
python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/multiconformer_tools/calc_dihedral.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb

#GET RESOLUTION
mtzmetadata=`phenix.mtz.dump "${pdb_name}.mtz"`
resrange=`grep "Resolution range:" <<< "${mtzmetadata}"`

echo "${resrange}"

res=`echo "${resrange}" | cut -d " " -f 4 | cut -c 1-5`
echo ${res}

#RUN OP
phenix.reduce -NOFLIP ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb > ${PDB}_qFit_H.pdb
make_methyl_df.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit_H.pdb --pdb ${PDB}
calc_OP.py ${output_dir}/${PDB}.dat ${PDB_dir}/${PDB}/${category}/${PDB}_qFit_H.pdb ${output_dir}/${PDB}_OP.out -r ${res} -b ${b_fac}

#__________________ACTIVATE ANOTHER ENVIORNMENT________________________________________________#
conda activate PE
##SASA
python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/multiconformer_tools/calc_sasa.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb




