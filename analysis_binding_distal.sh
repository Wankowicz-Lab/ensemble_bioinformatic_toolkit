

#find_largest_ligand.py ${PDB_dir}/${PDB}/${PDB}_qFit.pdb ${PDB}
#lig=$(cat "${output_dir}/${PDB}_ligand_name.txt")
#echo $lig
#find_close_residues.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb ${PDB} ${lig} 5.0
#find_close_residues.py ${PDB_dir}/${PDB}/${category}/${PDB}_qFit.pdb ${PDB} ${lig} 10.0
#python /wynton/group/fraser/swankowicz/script/subset_output_apo.py ${PDB} 7KQO -dist 5.0 -qFit=N -lig ${lig}
#python /wynton/group/fraser/swankowicz/script/subset_output_apo.py ${PDB} 7KQO -dist 10.0 -qFit=N -lig ${lig}
