#!/usr/bin/python

import pandas as pd
import json
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

############# Hardcoded Values
cancer = 'BRCA'
pmetric = 'overall_weighted_f1'
groups = ['gnosis', 'CF|All', 'AKLIMATE', 'nn', ['rfe15', 'fbedeBIC']]

# File formatted by BK (derived from Chris Wong's master matrix, but with RFE combined)
performance_matrix = 'src/brca_ak_sk.tsv'
scoring_matrix = 'src/model_scoring.features'

output_path = 'data/correlation/CORR_'
#############

performance_df = pd.read_csv(performance_matrix,sep='\t')

# For a group: select best model
ct = 1
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

    if 'AKLIMATE' in ftID: #  fix aklimate formatting issue
        if cancer == ftID.split('_')[-1]:
            ftID = '_'.join(ftID.split('_')[:-1])
    if 'CF|' in ftID: # fix cloud forest mapping
        ftID = cancer + ':' + ('_'.join(ftID.split('_')[:-1]).replace('Top_', 'Top ') )
    if 'gnosis' in ftID: # fix gnosis
        ftID = ':'.join([ftID.split('_')[0], ftID.split('_')[2], ftID.split('_')[1]])
    if 'nn_jg' in ftID: # fix subscope
        ftID = '_'.join(ftID.split('_')[:-1])
    print(ftID, ' selected as best model for group')


    # Retrieve the feature list for this particular ftID
    with open(scoring_matrix, 'r') as fh:
        data = json.load(fh)


    if type(group) == list:
        if 'fbed' in ftID:
            ft_set = []
            for k, v in data.items():
                if 'fbed' in k:
                    ft_set = ft_set + v
            ft_set = list(set(ft_set))

        if 'rfe' in ftID:
            ft_set = []
            for k, v in data.items():
                if 'RFE' in k:
                    ft_set = ft_set + v
            ft_set = list(set(ft_set))
    else:
        ft_set = data.get(ftID)
    print(len(ft_set), 'length of feature set')

    # save output as df
    if ct == 1: # create df if first group
        ft_set1 = ft_set
        ftID1 = ftID
    elif ct == 2:
        ft_set2 = ft_set
        ftID2 = ftID
    elif ct == 3:
        ft_set3 = ft_set
        ftID3 = ftID
    elif ct == 4:
        ft_set4 = ft_set
        ftID4 = ftID
    else:
        ft_df = pd.DataFrame.from_dict({ftID1: ft_set1, ftID2: ft_set2, ftID3: ft_set3, ftID4: ft_set4, ftID: ft_set}, orient='index')
#     else:
#         ft_df[ftID]=ft_set
    ct+=1

ft_df = ft_df.T
ft_df.to_csv(output_path + cancer + '_featurelist.tsv', sep = '\t',index=False)
ft_df
