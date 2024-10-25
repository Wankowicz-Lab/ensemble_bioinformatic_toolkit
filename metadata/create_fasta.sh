#!/bin/bash

#________________PDB INFO__________________________________
PDB_file=/dors/wankowicz_lab/stephanie/CDK/pdbs.txt
base_dir='/dors/wankowicz_lab/stephanie/CDK/'
file_name='CDK'

for i in {1..237}; do
  chain_info=$(python /dors/wankowicz_lab/stephanie/script/extact_PDB_seq.py ${base_dir}/${PDB}/${PDB}.pdb)
  echo $chain_info >> ${base_dir}/${file_name}_sequence.fa
done < $PDB_file
