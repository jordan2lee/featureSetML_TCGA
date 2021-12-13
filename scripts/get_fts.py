#!/usr/bin/python

import pandas as pd
from collections import Counter
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-t", "--tumor", help ="cancer cohort", required=True, type=str)
    parser.add_argument("-f1", "--file_fts", help ="classification performance file with all groups", required=True, type=str)
    parser.add_argument("-f2", "--file_top", help ="file path", required=True, type=str)
    parser.add_argument("-o", "--out", help ="output file", required=True, type=str)
    return parser.parse_args()

args = get_arguments()
cancer = args.tumor
file_fts = args.file_fts
file_output = args.out
main = args.file_top

# Read in feature sets
ft_df = pd.read_csv(file_fts, sep = '\t', index_col=0, low_memory=False)

# Fix concat ft names (index)
if 'N:METH:cg20568322:LMOD1:TssD7:NA_N:METH:cg07671858:SYNPO2L:TssD3203:Island_N:GEXP::HNF4A:3172:' in ft_df.index:
    # Issue ft and the corrected split into multiple features
    issue = 'N:METH:cg20568322:LMOD1:TssD7:NA_N:METH:cg07671858:SYNPO2L:TssD3203:Island_N:GEXP::HNF4A:3172:'
    split_ft = ['N:METH:cg20568322:LMOD1:TssD7:NA','N:METH:cg07671858:SYNPO2L:TssD3203:Island','N:GEXP::HNF4A:3172:']

    # Grab model column headers that have this feature present - ignoring 'total_number_of_lists' col
    s2 = ft_df.loc[issue,]
    in_models = []
    for i in range(0, s2.shape[0]):
        included = s2[i]
        m = s2.index[i]
        if included == '1' and m != 'total_number_of_lists':
            in_models.append(m)

    # Update df with single fts (conver 0 to 1)
    for ft in split_ft:
        for model in in_models:
            ft_df.at[ft,model]='1'

    # Drop issue ft
    ft_df = ft_df.drop([issue])

if 'N:METH:cg16128701:RP11-363E6.3:TssD233:Island_N:GEXP::ILF3:3609:' in ft_df.index:
    # Issue ft and the corrected split into multiple features
    issue = 'N:METH:cg16128701:RP11-363E6.3:TssD233:Island_N:GEXP::ILF3:3609:'
    split_ft = ['N:METH:cg16128701:RP11-363E6.3:TssD233:Island', 'N:GEXP::ILF3:3609:']

    # Grab model column headers that have this feature present - ignoring 'total_number_of_lists' col
    s2 = ft_df.loc[issue,]
    in_models = []
    for i in range(0, s2.shape[0]):
        included = s2[i]
        m = s2.index[i]
        if included == '1' and m != 'total_number_of_lists':
            in_models.append(m)

    # Update df with single fts (conver 0 to 1)
    for ft in split_ft:
        for model in in_models:
            ft_df.at[ft,model]='1'

    # Drop issue ft
    ft_df = ft_df.drop([issue])


if 'N:METH:cg04703174:KHDRBS2:TssD610:Shore_N:GEXP::SLC2A4RG:56731:' in ft_df.index:
    # Issue ft and the corrected split into multiple features
    issue = 'N:METH:cg04703174:KHDRBS2:TssD610:Shore_N:GEXP::SLC2A4RG:56731:'
    split_ft = ['N:METH:cg04703174:KHDRBS2:TssD610:Shore','N:GEXP::SLC2A4RG:56731:']

    # Grab model column headers that have this feature present - ignoring 'total_number_of_lists' col
    s2 = ft_df.loc[issue,]
    in_models = []
    for i in range(0, s2.shape[0]):
        included = s2[i]
        m = s2.index[i]
        if included == '1' and m != 'total_number_of_lists':
            in_models.append(m)

    # Update df with single fts (conver 0 to 1)
    for ft in split_ft:
        for model in in_models:
            ft_df.at[ft,model]='1'

    # Drop issue ft
    ft_df = ft_df.drop([issue])




# Read in file
df = pd.read_csv(main, sep='\t')
df = df[['cohort','featureID', 'Mean', 'performance_metric', 'feature_list', 'feature_list_method']]

# Pull and reformat top models
s1 = df[df['cohort']==cancer].reset_index(drop=True)

# If a team tie pick the first and record that there were ties
best = list(s1['featureID'])
# If intra-team ties, pick the first on
if len(best) >5:
    new_best = []
    for t in set(s1['feature_list_method']):
        model = s1[s1['feature_list_method']==t]['featureID'].values[0]
        print(model)
        new_best.append(model)
    best = new_best
assert len(best) == 5, 'Error: more than one model per team is ranked as best'

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
