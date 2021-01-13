#!/usr/bin/env python3

import pandas as pd
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-t", "--tumor", help ="cancer cohort", required=True, type=str)
    parser.add_argument("-f1", "--file1_transformedfts", help ="file of consolidated features from best model", required=True, type=str)
    parser.add_argument("-f2", "--file2_raw", help ="file of raw molecular values", required=True, type=str)
    parser.add_argument("-o1", "--out1", help ="output file 1 by team", required=True, type=str)
    parser.add_argument("-o2", "--out2", help ="output file 2 by value", required=True, type=str)
    return parser.parse_args()

args = get_arguments()
cancer = args.tumor
file = args.file1_transformedfts
file_tarball = args.file2_raw
out1 = args.out1
out2 = args.out2

# open file and rm fts no model used
df = pd.read_csv(file, sep='\t', index_col=0)
df = df[(df.T != 0).any()]
df['Total'] = df.sum(axis=1) # add total col
df.to_csv(out1, sep='\t')

# Get values of original matrix for features of interest
tarball = pd.read_csv(file_tarball, sep = '\t', index_col = 0)
fts = list(df.index)
mat = tarball[['Labels']+ fts]
mat.to_csv(out2, sep='\t')
