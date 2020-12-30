#!/usr/bin/env python

#####################################
# Purpose: run recursive feature elimination down to 15 features
  # svc with linear kernal
  # eliminate until 15 fts remain
  # each round of removing will remove 1% of features
  # MUST USE PYTHON THAT USES SORTED DICTIONARY
#####################################

import pandas as pd
import os
import re
from sklearn.feature_selection import RFE
from sklearn.svm import SVC
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-ip", "--input_path", help ="path where input file located (don't include filename)", required=True, type=str)
    parser.add_argument("-i_matrix", "--input_matrix", help ="name of input tumor matrix of ft and labels", required=True, type=str)
    parser.add_argument("-op", "--output_path", help ="path where output files created (don't include filename)", required=True, type=str)
    parser.add_argument("-t", "--tumor", help ="prefix of tumor type (THYM) for output file naming", required=True, type=str)
    return parser.parse_args()

args = get_arguments()
ip = args.input_path
i_matrix = args.input_matrix
op = args.output_path
tumor = args.tumor

################### Use for man. entry
# ip = "/Users/leejor/Ellrott_Lab/02_ML/02_more-tumtypes/data/raw/forupload"
# i_matrix = "THYM_20190814.tsv" #input feature/label file
# op = "/Users/leejor/Ellrott_Lab/02_ML/02_more-tumtypes/data/thym/ft-sel"
# tumor = "THYM"
###################

ft_method = "rfe15"
out_name = tumor + "_" + ft_method + "--combined_platform.tsv"

###############################
# SET UP
###############################
os.chdir(ip)
raw = pd.read_csv(i_matrix, delimiter= "\t", index_col=0)

bal_acc_scoresRFE15 = [] # list of balance accur scores

data = raw # rename to avoid potential alter to raw data

######################################
#  CREATE LABEL AND FEATURE MATRICES for all platforms
#####################################
names = data.columns
assert len(names) == data.shape[1]

# # Create list of col headers
# MUT_names = [item for item in names if re.match('[A-Za-z0-9]:MUTA:', item)]
# CNVR_names = [item for item in names if re.match('[A-Za-z0-9]:CNVR::', item)]
# GEXP_names = [item for item in names if re.match('[A-Za-z0-9]:GEXP::', item)]
# MI_names = [item for item in names if re.match('[A-Za-z0-9]:MIR::', item)]
# MET_names = [item for item in names if re.match('[A-Za-z0-9]:METH:', item)]
# assert (len(MUT_names)+len(CNVR_names)+len(GEXP_names)+len(MI_names)+len(MET_names)+1) == data.shape[1]

# Subset by PLATFORM = CREATE FEATURE MATRICES
ALL_ft_matrix = data
# MUT_ft_matrix = data[MUT_names]
# CNVR_ft_matrix = data[CNVR_names]
# GEXP_ft_matrix = data[GEXP_names]
# MI_ft_matrix = data[MI_names]
# MET_ft_matrix = data[MET_names]

# CREATE LABEL MATRIX
#same for all platform subsets because didn't alter rows
y = data["Labels"]
assert data.shape[0] == len(y)

# rm cols that aren't fts
# print('first few col names {}'.format(ALL_ft_matrix.columns[0:3])) #sanity check
ALL_ft_matrix = ALL_ft_matrix.iloc[:,1:]
# print('col names post filter {}'.format(ALL_ft_matrix.columns[0:3])) #sanity check

# SET UP ALGORITHM FOR ALL PLATFORMS
fs = RFE(SVC(kernel="linear",random_state=1), 15, 0.01)

######################################
# FEATURE SELECTION : ALL
#####################################
platform = 'ALL'
print('##### starting ALL: rfe15 #####')

# Create dictionary skeleton
k=['Feature', 'Rank', 'Ft_method']
v=[{},{}, {}]
res_dictionary= dict(zip(k,v))

if ALL_ft_matrix.shape[1] > 0:
    # Run algorithm
    fs.fit(ALL_ft_matrix, y)

    # ALL RESULTS:
    # Get entries to pop dictionary
    list_ranks = fs.ranking_
    list_feature = ALL_ft_matrix.columns
    assert len(list_ranks)==len(list_feature)
    # Populate dictionary
    for i in range(0, len(list_ranks)):
        res_dictionary['Feature'][i]= list_feature[i]
        res_dictionary['Rank'][i]= list_ranks[i]
        res_dictionary['Ft_method'][i]=str(fs.estimator).replace('\n','').replace(' ','')

    # Pause: Save all features+rankings to a .tsv
    df_all_FEATURES = pd.DataFrame.from_dict(res_dictionary)
    df_all_FEATURES = df_all_FEATURES.sort_values(by=['Rank']).reset_index(drop=True)
    df_all_FEATURES.to_csv(op+'summaryfullranks-'+tumor+'_rfe15--'+platform+"--combined_platform.txt", sep='\t', index=False)
    print('output saved at: ', op+'summaryfullranks-'+tumor+'_rfe15--'+platform+"--combined_platform.txt")

    # Create list of top 15 features
    ALL_select = list(df_all_FEATURES.iloc[0:15,0])

if ALL_ft_matrix.shape[1] == 0:
    print('Input matrix Did NOT have any ALL features')
    ALL_select = [] # create empty list

###############################
# Create selected features tsv
################################
# summary = pd.DataFrame.from_dict({"MUT": MUT_select, "CNVR": CNVR_select, "GEXP": GEXP_select,"MI": MI_select, "MET": MET_select}, orient='index').T
with open(op+out_name, 'w') as out:
    out.write('feature\n')
    for a in ALL_select:
        out.write(a + '\n')
# ALL_select.to_csv(op+out_name, sep="\t", index = False)
print('output saved at: ', op+out_name)
print("all complete")
