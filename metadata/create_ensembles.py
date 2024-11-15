#!/usr/bin/env python3
import sys,os
from EnsemblePDB import search, build, refine, analyze
from EnsemblePDB.utils import file_management
import pandas as pd
import argparse

'''
Example: python ensemblepdb.py --dir=/dors/wankowicz_lab/stephanie/serine_protease --protease='elastase' --seq='VVGGTEAQRNSWPSQISLQYRSGSSWAHTCGGTLIRQNWVMTAAHCVDRELTFRVVVGEHNLNQNNGTEQYVGVQKIVVHPYWNTDDVAAGYDIALLRLAQSVTLNSYVQLGVLPRAGTILANNSPCYITGWGLTRTNGQLAQTLQQAYLPTVDYAICSSSSYWGSTVKNSMVCAGGDGVRSGCQGDSGGPLHCLVNGQYAVHGVTSFVSRLGCNVTRKPTVFTRVSAYISWINNVIASN' --ref_pdb='6EST' --ref_chain='A' --protein_name='elastate' --organism='Porcine'
'''

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process some input parameters.")
    parser.add_argument("--dir", required=True, help="Direct'ory path")
    parser.add_argument("--protease", required=True, help="Prorease parameter")
    parser.add_argument("--seq", required=True, help="Sequence parameter")
    parser.add_argument("--ref_pdb", required=True, help="Reference PDB file")
    parser.add_argument("--ref_chain", required=True, help="Reference chain")
    parser.add_argument("--protein_name", required=True, help="Name of the protein (lowercase)")
    parser.add_argument("--organism", required=True, help="Name of the organism")
    return parser.parse_args()

args = parse_arguments()
output_dir = args.dir
protease = args.protease
seq = args.seq
ref_pdb = args.ref_pdb
ref_chains = args.ref_chain


summary = search.query_pdb.query_by_sequence(sequences=[seq], label=protease,macromolecule_name=None, seq_id=0.95,xray_only=True, evalue=10,pairwise_align_options = {'alignment_type':"local",'gap_open_score':-0.5, 'gap_extend_score':-0.1},output_dir=output_dir,min_length=50)
data = search.get_data.get_pdb_data(f'{output_dir}/seq_search_output_{protease}.csv', label=protease)
print(args.organism)
organism_list = [args.organism]

data = refine.filter.filter_pdb_data(summary_report_csv=f'{output_dir}/summary_report_{protease}.csv'
, protein_name=args.protein_name, exclude_terms=['zymogen','radiation damage','proenzyme','chymotrypsinogen'], max_res=2.0, organism=organism_list, max_muts=3)

data = refine.check_chains.check_multimer_chains(ensemble_csv=f'{output_dir}/summary_report_{protease}_filtered.csv', allowance=10)

data = build.renumber.download_and_renumber(summary_report_csv=f'{output_dir}/summary_report_{protease}_filtered_checked.csv', reference_pdb = ref_pdb,reference_chains=ref_chains,multimer_save=True, combine_chains_by='all')
