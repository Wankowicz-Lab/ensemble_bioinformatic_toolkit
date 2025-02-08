#!/bin/bash
#SBATCH --mem=8G  # Request 4GB of memory
#SBATCH --cpus-per-task=1  # Request 1 CPUs per task
#SBATCH --time=60:00:00  # Request a walltime of 2 hours
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1%1  # Set up a job array with the number of PDBs you have
#SBATCH --output=CATH._stdout # Set up an output file with 


#__________________SOURCE PHENIX/QFIT________________________________________________#
source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit


# Set the parent directory where the two-digit folders are located
parent_dir="/dors/wankowicz_lab/PDBRedo/pdb-redo3/"
output_dir="/dors/wankowicz_lab/PDBRedo/pdb-redo3/CATH"

# Loop through each two-digit folder
for two_digit_folder in "$parent_dir"/*/; do
    # Ensure we are in the correct folder
    if [[ -d "$two_digit_folder" ]]; then
        echo $two_digit_folder
        cd $two_digit_folder
        # Loop through each folder inside the two-digit folder
        for folder in "$two_digit_folder"*/; do
            conda activate qfit
            if [[ -d "$folder" ]]; then
                echo $folder
                cd $folder
                pdb=$(basename "$folder")
                cd ${output_dir}
                wget https://www.ebi.ac.uk/pdbe/api/mappings/cath/${pdb}
            fi
        done
     fi
done
