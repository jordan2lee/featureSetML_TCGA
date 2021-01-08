#!/usr/bin/python

import pandas as pd
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-t", "--tumor", help ="cancer cohort", required=True, type=str)
    parser.add_argument("-m", "--metric", help ="classification performance metric", required=True, type=str)
    parser.add_argument("-f1", "--file1_fts", help ="classification performance file with all groups", required=True, type=str)
    parser.add_argument("-f2", "--file2_perform", help ="feature set file with all groups", required=True, type=str)
    parser.add_argument("-o", "--out", help ="output file", required=True, type=str)
    return parser.parse_args()

args = get_arguments()
cancer = args.tumor
pmetric = args.metric
file_fts = args.file1_fts
file_preds = args.file2_perform
file_output = args.out

############# Hardcoded Object
groups = ['gnosis', 'CF|All', 'AKLIMATE', 'nn', ['rfe15', 'fbedeBIC']]
#############

performance_df = pd.read_csv(file_preds, sep = '\t', low_memory=False)

# For a group: select best model
ct = 1
best = []
for group in groups:
    # all models from that group for the 3 criteria
    if type(group)== list:
        subset = performance_df[performance_df['feature_list_method'].isin(['rfe15', 'fbedeBIC'])]
        subset = subset[subset['cohort'] == cancer]
        subset = subset[subset['performance_metric']== pmetric].reset_index(drop=True)
        ftID = subset.sort_values(by='Mean', ascending=False).reset_index(drop=True)['featureID'][0]
    else:
        subset = performance_df[performance_df['feature_list_method'] == group]
        subset = subset[subset['cohort'] == cancer]
        subset = subset[subset['performance_metric']== pmetric].reset_index(drop=True)


    # Grab the name of the model with highest MEAN performance metric
    ftID = subset.sort_values(by='Mean', ascending=False).reset_index(drop=True)['featureID'][0]

    ##
    # Fix naming of CloudForest (to match ft file)
    if "CF|" in ftID:
        ftID='Top '.join(ftID.split('Top_'))
    ##

    best.append(ftID)
    print(ftID, ' selected as best model for group')

assert len(best) == 5, 'best model not found for all 5 groups'


# Subset feature matrix for only best model per team
ft_df = pd.read_csv(file_fts, sep = '\t', index_col=0, low_memory=False)
ft_df = ft_df[best].drop(['feature_list_method','feature_list_cohort','feature_list_size']) # rm annotat rows
ft_df.to_csv(file_output, sep='\t')
