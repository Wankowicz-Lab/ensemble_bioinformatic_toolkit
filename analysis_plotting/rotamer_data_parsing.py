#packages
import pandas as pd
import numpy as np
import os
import glob
import sys



# ## LOAD IN ROTAMERS
all_files = glob.glob("./rotamer/*_rotamer_output.txt")
li = []
for filename in all_files:
    try:
        df = pd.read_csv(filename, index_col=None, header=0, sep=':')
        df['PDB'] = filename[10:15]
        li.append(df)
    except pd.errors.EmptyDataError:
        continue

rotamer = pd.concat(li, axis=0, ignore_index=True)
rotamer = rotamer[rotamer['residue']!= 'SUMMARY'].reset_index()
split = rotamer['residue'].str.split(" ")
for i in range(0, (len(rotamer.index)-1)):
    rotamer.loc[i,'chain'] = split[i][1]
    STUPID = str(rotamer.loc[i,'residue'])[3:6]
    tmp = []
    try:
        tmp = (int(''.join(i for i in STUPID if i.isdigit())))
    except:
        newstr = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in STUPID)
        tmp = [float(i) for i in newstr.split()]
    rotamer.loc[i, 'resi'] = tmp
