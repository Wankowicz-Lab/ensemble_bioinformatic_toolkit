#!/usr/bin/env python

import pandas as pd
import argparse

def parse_log(qfit_log, pdb, qFit):
    rval = pd.DataFrame()
    log = open(qfit_log, 'r')
    rval.loc[1,'PDB'] = pdb
    for line in log:
        if line.startswith('Final R-work'):
            rval.loc[1,'Rwork'] = line.split('=')[1][1:6]
            rval.loc[1,'Rfree'] = line.split('=')[2][1:6]
    rval.to_csv(pdb + qFit +'_rvalues.csv', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('qFit_Log_File')
    parser.add_argument('PDB')
    parser.add_argument("--qFit", action='store_true')
    args = parser.parse_args()
    qFit = '_qFit' if args.qFit else ''
    parse_log(args.qFit_Log_File, args.PDB, qFit)

