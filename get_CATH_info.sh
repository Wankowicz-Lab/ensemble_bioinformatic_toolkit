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
parent_dir="/dors/wankowicz_lab/all_pdb/60001_70000"
output_dir="/dors/wankowicz_lab/all_pdb/CATH/"
uniprot_dir="/dors/wankowicz_lab/all_pdb/uniprot/"

# Loop through each two-digit folder
for folder in "$parent_dir"/*; do
            if [[ -d "$folder" ]]; then
                echo $folder
                cd $folder
                pdb=$(basename "$folder")
                cd ${output_dir}
                wget https://www.ebi.ac.uk/pdbe/api/mappings/cath/${pdb}
                cd ${uniprot_dir}
                wget https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/${pdb}
            fi
done
