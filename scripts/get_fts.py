#!/usr/bin/python

import pandas as pd

############# Hardcoded Values
cancer = 'BRCA'
pmetric = 'overall_weighted_f1'
groups = ['gnosis', 'CF|All', 'AKLIMATE', 'nn', ['rfe15', 'fbedeBIC']]

file_preds = 'src/feature_list_with_performance_with_subtype_names_20200828.tsv.gz'
file_fts = 'src/collected_features_matrix_20200722.tsv.gz'
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
ft_df.to_csv('data/exact_match/best_models_{}.tsv'.format(cancer), sep='\t')
