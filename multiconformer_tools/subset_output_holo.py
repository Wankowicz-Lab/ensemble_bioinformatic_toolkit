import numpy as np
import pandas as pd
import argparse
import os
import sys

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("holo", type=str, help="name of holo")
    p.add_argument("-dist", type=str,
                   help="distance of close residues")
    p.add_argument("-qFit", type=str)
    p.add_argument("-lig", type=str)
    args = p.parse_args()
    return args


def main():
    args = parse_args()
    if not args.qFit == None:
            qFit = '_qFit'
    else:
            qFit = ''
    rmsd = None
    B_factor = None
    sasa = None
    rotamer = None
    order = None

    
    try:
        B_factor = pd.read_csv(args.apo + qFit + "_B_factors.csv")
        B_factor['AA'] = B_factor.AA.str.replace('[','')
        B_factor['AA'] = B_factor.AA.str.replace(']','')
        B_factor['Chain'] = B_factor.Chain.str.replace(']','')
        B_factor['Chain'] = B_factor.Chain.str.replace('[','')
        B_factor['resseq'] = B_factor.resseq.str.replace('[','')
        B_factor['resseq'] = B_factor.resseq.str.replace(']','')
        B_factor['Chain'] = B_factor.Chain.str.replace("\'", '')
        B_factor['resseq'] = B_factor['resseq'].astype(int)

    except IOError:
        pass
    
    try:
    	sasa = pd.read_csv(args.apo + qFit + "_sasa.csv", index_col=0)
    except IOError:
       pass  
    
    try:
    	hbond = pd.read_csv(args.apo + "_qFit_hbonds.csv", index_col=0)
    except IOError:
       pass  
    

    try:
       rotamer = pd.read_csv(args.apo + qFit + "_rotamer_output.txt", sep = ':')
       split = rotamer['residue'].str.split(" ") 
       for i in range(0,len(rotamer.index)-1):
           rotamer.loc[i,'chain'] = split[i][1]
           STUPID = str(rotamer.loc[i,'residue'])[3:8]
           rotamer.loc[i,'resi'] = [int(s) for s in STUPID.split() if s.isdigit()]
    except IOError:
       pass
    
    try:
        order = pd.read_csv(args.holo + '_OP.out', sep=',')
    except IOError:
        pass
    close_res = pd.read_csv(args.holo + "_" + args.lig + "_" + args.dist + "_closeres.csv")
    close_res = close_res.drop_duplicates()
    

    li_b = []
    li_r = []
    li_sasa = []
    li_rotamer = []
    li_order = []
    li_hbond = []    

    for i in close_res.chain.unique():
        output = close_res[close_res['chain'] == i]
        residue = output.res_id.unique()
    
        if B_factor is not None:
          b_s = B_factor[(B_factor['Chain'] == i) & (B_factor['resseq'].isin(residue))]
          li_b.append(b_s)

        rot_s = rotamer[(rotamer['chain'] == i) & (rotamer['resi'].isin(residue))]
        li_rotamer.append(rot_s)
    
        sasa_s = sasa.loc[(sasa['Chain'] == i) & (sasa['Residue ID'].isin(residue))]
        li_sasa.append(sasa_s)

        order_s = order[(order['chain'] == i) & (order['resi'].isin(residue))]
        li_order.append(order_s)

        hbond_s = hbond[((hbond['donor_chain'] == i) & (hbond['donor_residue_number'].isin(residue))) |
                        ((hbond['acceptor_chain'] == i) & (hbond['acceptor_residue_number'].isin(residue)))]
        li_hbond.append(hbond_s)

    order_subset = pd.concat(li_order, axis=0, ignore_index=True)
    
    
    if B_factor is not None:
       b_factor_subset = pd.concat(li_b, axis=0, ignore_index=True)
       b_factor_subset.to_csv(args.holo + '_' + args.dist + '_bfactor_subset.csv', index=False)
    rotamer_subset = pd.concat(li_rotamer, axis=0, ignore_index=True)
    sasa_subset = pd.concat(li_sasa, axis=0, ignore_index=True)
    
    b_factor_subset.to_csv(args.holo + '_' + args.dist + '_bfactor_subset.csv', index=False)
    order_subset.to_csv(args.holo + '_' + args.dist + '_order_param_subset.csv', index=False)
    sasa_subset.to_csv(args.holo + '_' + args.dist + '_sasa_subset.csv', index=False)
    rotamer_subset.to_csv(args.holo + '_' + args.dist + '_rotamer_subset.csv', index=False)

if __name__ == "__main__":
   main()