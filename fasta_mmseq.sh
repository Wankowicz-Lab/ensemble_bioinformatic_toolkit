#!/bin/bash
#SBATCH --mem=10G  # Request 10GB of memory
#SBATCH --cpus-per-task=1  # Request 1 CPUs per task
#SBATCH --time=24:00:00  # Request a walltime of 24 hours
#SBATCH --array=1-1%1  # Set up a job array with 1 job
#SBATCH --output=bioinformatics._stdout

#__________________SOURCE PHENIX/QFIT________________________________________________#
source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit

#________________PDB INFO__________________________________
PDB_file="/dors/wankowicz_lab/all_pdb/1_10000/pdb_1_10000.txt"
PDB_dir="/dors/wankowicz_lab/all_pdb/1_10000/" 
OUTPUT_FASTA="PDB_seq.fasta"
MMSEQ_RESULT_dir="/dors/wankowicz_lab/mmseq_output"  # Update 

> $OUTPUT_FASTA  

# Step 1: Extract sequences for each PDB
while IFS= read -r PDB_ID; do
    PDB_FILE="$PDB_DIR/${PDB_ID}.pdb"
    if [ -f "$PDB_FILE" ]; then
        python get_seq.py "$PDB_FILE" > "${PDB_ID}_seq.txt"
    else
        echo "Warning: PDB file $PDB_FILE not found."
        continue
    fi

done < "$PDB_TXT"

# Step 2: Generate FASTA file from sequence files
for PDB_ID in $(cat $PDB_TXT); do
    if [ -f "${PDB_ID}_seq.txt" ]; then
        python get_fasta.py "${PDB_ID}_seq.txt" temp.fasta
        cat temp.fasta >> $OUTPUT_FASTA
        rm temp.fasta
    fi
done

# Step 3: Run MMseqs
mmseqs easy-search $OUTPUT_FASTA $OUTPUT_FASTA $MMSEQ_RESULT_DIR --threads 4 --format-output "query,target,pident,alnlen"

echo "Processing complete. Sequences are saved in $OUTPUT_FASTA and MMseqs output in $MMSEQ_RESULT_DIR."
