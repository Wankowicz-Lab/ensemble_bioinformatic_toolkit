#!/bin/bash
#Stephanie Wankowicz
#04/25/2019

source /dors/wankowicz_lab/phenix-installer-dev-5366-intel-linux-2.6-x86_64-centos6/phenix-dev-5366/phenix_env.sh

#________________________________________________INPUTS________________________________________________
pdb_filelist=/dors/wankowicz_lab/stephanie/serine_protease/trypsin/Renumbered_unaligned_pdbs/pdbs.txt
base_folder='/dors/wankowicz_lab/stephanie/serine_protease/trypsin/Renumbered_unaligned_pdbs/'
while read -r line; do
  PDB=$line
  cd $base_folder
  if [ -d "/$PDB" ]; then
    echo "Folder exists."
  else
    mkdir $PDB
  fi
  cd $PDB
  phenix.fetch_pdb ${PDB}
  phenix.fetch_pdb ${PDB} -x
  mv ../${PDB}.pdb .
