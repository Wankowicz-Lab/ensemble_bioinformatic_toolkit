#!/bin/bash
#SBATCH --mem=10G  # Request 4GB of memory
#SBATCH --cpus-per-task=1  # Request 1 CPUs per task
#SBATCH --time=24:00:00  # Request a walltime of 2 hours
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1%1  # Set up a job array with the number of PDBs you have
#SBATCH --output=bioinformatics._stdout # Set up an output file with 


#__________________SOURCE PHENIX/QFIT________________________________________________#
source /dors/wankowicz_lab/phenix-installer-dev-5366-intel-linux-2.6-x86_64-centos6/phenix-dev-5366/phenix_env.sh
export PHENIX_OVERWRITE_ALL=true

source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit


#________________PDB INFO__________________________________
PDB_file=/dors/wankowicz_lab/all_pdb/10001_20000/pdb_10001_20000.txt ## FIX AS NEEDED
PDB_dir='/dors/wankowicz_lab/all_pdb/10001_20000/'  ## FIX AS NEEDED
output_dir='/dors/wankowicz_lab/all_pdb/10001_20000/inforamtics_output' ## FIX AS NEEDED

cd ${output_dir}
while IFS= read -r PDB; do
    cd ${PDB_dir}/${PDB}
    echo ${PDB}
    if [[ ! -f "${PDB}_qFit.pdb" ]]; then
        echo "${PDB}_qFit.pdb not found."
        continue  # Skip to the next PDB if the file is not found
    fi
cd ${output_dir}
#PDB=$(cat $PDB_file | head -n 1 | tail -n 1)
#echo $PDB

#cp ${PDB_dir}/${PDB}/${category}/${PDB}.updated_refine_001_qFit.pdb ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb
#cp ${PDB_dir}/${PDB}/${category}/${PDB}.updated_refine_001_qFit.log ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.log

#source /sb/sbgrid/programs/sbgrid.shrc
#python /dors/wankowicz_lab/stephanie/script/align_structures.py ${PDB_dir}/${PDB}/${PDB}_qFit.pdb ${PDB_dir}/6EAU/6EAU_qFit.pdb  

#______________________PDB STATS_______________________
source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit
##NUM ALTLOC
b_fac=$(b_factor.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb --pdb=${PDB})

##ROTAMER
phenix.rotalyze model=${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb outliers_only=False > ${output_dir}/${PDB}_rotamer_output.txt

#RUN PDBTOOLS TO REMOVE HYDROGENS
phenix.pdbtools remove="name H*" ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb

#RUN CALC CHI ON MODIFIED PDB
calc_chi ${PDB}_qFit_modified.pdb ${PDB}_qFit

## RFREE
python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/multiconformer_tools/single_parse_log.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.log  ${PDB}

##VALIDATION SCRIPTS
mmtbx.validation_summary ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb > ${PDB}_validation.txt

##HBOND
python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/multiconformer_tools/weighted_hydrogen_bond.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb

##DIHEDRAL
#python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/multiconformer_tools/calc_dihedral.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb

#phenix.reduce -NOFLIP ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb > ${PDB}_qFit_H.pdb
make_methyl_df.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb --pdb ${PDB}
calc_OP.py ${output_dir}/${PDB}.dat ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb ${output_dir}/${PDB}_OP.out -r 1.5 -b $b_fac
python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/analysis_plotting/rename_b_factor.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb ${output_dir}/${PDB}_OP.out ${output_dir}/${PDB}_OP.pdb --column_name s2calc


#__________________ACTIVATE ANOTHER ENVIORNMENT________________________________________________#
conda activate PE
##SASA
python /dors/wankowicz_lab/ensemble_bioinformatic_toolkit/multiconformer_tools/calc_sasa.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb

##RUN FPOCKET


## RUN PYMOL CLOSE SCRIPT


## SUBSET PYTHON SCRIPT

done < "$PDB_file"

