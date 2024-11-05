#!/bin/bash


#__________________SOURCE PHENIX/QFIT________________________________________________#
source /dors/wankowicz_lab/phenix-installer-dev-5366-intel-linux-2.6-x86_64-centos6/phenix-dev-5366/phenix_env.sh
export PHENIX_OVERWRITE_ALL=true

#________________PDB INFO__________________________________
PDB_file=/dors/wankowicz_lab/stephanie/kemp_elim/pdbs.txt
base_dir='/dors/wankowicz_lab/stephanie/kemp_elim/'
file_name='kemp'

cd ${base_dir}
echo 'PDB' 'Resolution' 'SpaceGroup' 'UnitCell_L1' 'UnitCell_L2' 'UnitCell_L3' 'UnitCell_A1' 'UnitCell_A2' 'UnitCell_A3' 'Lig' 'qFit_Rwork' 'qFit_Rfree' 'Chain_A_Seq' 'Chain_A_Res' 'Chain_B_Seq' 'Chain_A_Res'  >> ${base_dir}/${file_name}_space_unit_reso.txt


for i in {1..47}; do
  PDB=$(cat $PDB_file | head -n $i | tail -n 1)
  echo ${PDB}
  phenix.cif_as_mtz ${PDB}-sf.cif --merge
  mv ${PDB}-sf.mtz ${PDB}.mtz
  phenix.mtz.dump ${base_dir}/${PDB}/${PDB}.mtz > ${base_dir}/${PDB}/${PDB}.dump
  SPACE1=$(grep "^Space group symbol from file:" ${base_dir}/${PDB}/${PDB}.dump | awk '{print $6,$7}')
  UNIT1=$(grep "Unit cell:" ${base_dir}/${PDB}/${PDB}.dump | tail -n 1 | sed "s/[(),]//g" | awk '{print $3,$4,$5,$6,$7,$8}')
  RESO1=$(grep "^Resolution" ${base_dir}/${PDB}/${PDB}.dump | head -n 1 | awk '{print $4}')
  python /dors/wankowicz_lab/stephanie/script/single_parse_log.py ${base_dir}/${PDB}/${PDB}_qFit.log  ${PDB}
  log_file="${base_dir}/${PDB}_rvalues.csv"
  read -r log_info < <(tail -n +2 "$log_file" | head -n 1)
  echo $log_info
  IFS=',' read -r pdb_id pdb_rwork pdb_rfree <<< "$log_info"  # Split the CSV into variables  

  if [ -z "$pdb_rwork" ]; then
    pdb_rwork="no_rwork"
  fi

  if [ -z "$pdb_rfree" ]; then
    pdb_rwork="no_rfree"
  fi
  find_largest_ligand.py ${base_dir}/${PDB}/${PDB}.pdb ${PDB}
  lig=$(cat "${base_dir}/${PDB}_ligand_name.txt")
  if [ -z "$lig" ]; then
    lig="no_lig"
  fi
  chain_info=$(python /dors/wankowicz_lab/stephanie/script/extact_PDB_seq.py ${base_dir}/${PDB}/${PDB}.pdb)
  echo $PDB $RESO1 $SPACE1 $UNIT1 $pdb_rwork $pdb_rfree $lig $chain_info >> ${base_dir}/${file_name}_space_unit_reso.txt
done
