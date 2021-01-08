#!/usr/bin/env python3

import pandas as pd
import scipy.stats
import argparse


def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--file", help ="file ex BRCA_fts_by_VALUE.tsv", required=True, type=str)
    parser.add_argument("-o", "--outfile", help ="outfile like data/heatmap/BRCA_vals_meth.tsv", required=True, type=str)
    parser.add_argument("-d", "--datatype", help ="data type prefix as seen in tarball ex N:METH", required=True, type=str)
    return parser.parse_args()
args = get_arguments()


###### hardcoded
# file = 'data/heatmap/BRCA_fts_by_VALUE.tsv'
# datatype = 'N:METH'
# out = 'data/heatmap/BRCA_vals_meth.tsv'
######

# Format row:col:val - ft:sample:tarballval
df = pd.read_csv(args.file, sep='\t', index_col=0)
df = df.drop(['Labels'], axis=1) # rm col

# filter for datatype cols
keep =[]
for c in df.columns:
    if c.startswith(args.datatype):
        keep.append(c)
print(len(keep), '{} features kept'.format(args.datatype))
df = df[keep]

df.to_csv(args.outfile,sep='\t')
