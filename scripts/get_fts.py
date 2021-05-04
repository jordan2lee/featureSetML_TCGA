#!/usr/bin/python

import pandas as pd
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-t", "--tumor", help ="cancer cohort", required=True, type=str)
    parser.add_argument("-m", "--metric", help ="classification performance metric", required=True, type=str)
    parser.add_argument("-fil", "--filters", help ="none or integer for feature set filter max size to be considered for best model", required=True, type=str)
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
filters = args.filters

############# Hardcoded Object
# groups = ['gnosis', 'CF|All', 'AKLIMATE', 'nn', ['rfe15', 'fbedeBIC']]
groups = ['jadbio', 'CF', 'AKLIMATE', 'subSCOPE', 'skgrid']
#############

performance_df = pd.read_csv(file_preds, sep = '\t', low_memory=False)

# For a group: select best model
ct = 1
best = []
for group in groups:
    print(group)
    subset = performance_df[performance_df['feature_list_method'] == group]
    subset = subset[subset['cohort'] == cancer]
    subset = subset[subset['performance_metric'] == pmetric].reset_index(drop=True)
    if filters != 'none':
        max_ft_size = int(filters)
        subset["total_features"] = pd.to_numeric(subset["total_features"])
        subset = subset[subset['total_features'] <= max_ft_size].reset_index(drop=True)
    subset = subset.sort_values(by='Mean', ascending=False).reset_index(drop=True)
    # Grab the name of the model with highest MEAN performance metric
    # if found at least one model
    if subset.shape[0] > 0:
        ftID = subset.sort_values(by='Mean', ascending=False).reset_index(drop=True)['featureID'][0]
    #else no models fitting the above filters, need to finish
    else:
        ftID = 'NO_MODEL_MATCH_' + group

    ##
    # Fix naming of CloudForest (to match ft file)
    if "CF_" in ftID:
        ftID='Top_'.join(ftID.split('Top '))
    ##

    best.append(ftID)
    print(ftID, ' selected as best model for group')


assert len(best) == 5, 'best model not found for all 5 groups'


# Subset feature matrix for only best model per team
ft_df = pd.read_csv(file_fts, sep = '\t', index_col=0, low_memory=False)

# likely want to refactor below to be more streamline and short #

# Handle instances where no model for a team found
best_orders = dict()
best_models_present = []
best_models_not_present = dict()
for i in range(0, len(best)):
    # Create Index for each model - to maintain order in final df
    best_orders[i]= best[i]
    m = best[i]
    # Add to list models that actually exist
    if not m.startswith('NO_MODEL_MATCH'):
        if m.startswith('subSCOPE'):
            best_models_present.append(m+'_'+cancer)
        else:
            best_models_present.append(m)
    else:
        # if subSCOPE then add nn_jg to ensure uniq naming
        if m == 'NO_MODEL_MATCH_nn':
            m = 'NO_MODEL_MATCH_nn_jg'
        # add it
        best_models_not_present[i]=m

# Create output df of teams with models found
ft_df = ft_df[best_models_present].drop(['feature_list_method','feature_list_cohort','feature_list_size']) # rm annotat rows

# if present, Add back in team without models
if len(best_models_not_present) != 0:
    for i, model in best_models_not_present.items():
        ft_df.insert(i, model, [0] * ft_df.shape[0])
ft_df.to_csv(file_output, sep='\t')
