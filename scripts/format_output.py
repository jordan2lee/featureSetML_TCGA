#!/usr/bin/env python3


import pandas as pd

cancer_list = [
    "ACC",
    "BLCA",
    "BRCA",
    "CESC",
    "COADREAD",
    "ESCC",
    "GEA",
    "HNSC",
    "KIRCKICH",
    "KIRP",
    "LGGGBM",
    "LIHCCHOL",
    "LUAD",
    "LUSC",
    "MESO",
    "OV",
    "PAAD",
    "PCPG",
    "PRAD",
    "SARC",
    "SKCM",
    "TGCT",
    "THCA",
    "THYM",
    "UCEC",
    "UVM"
]


###
# function
###
def get_ft_str(file):
    '''create ft string - easier to write output'''
    df = pd.read_csv(file, sep='\t')
    fts = list(df.iloc[:, 0]) # grab fts
    fts = str(fts).replace("\'", "\"").replace(' ', '')
    return fts



###
# write output file
###
methods = ['fbedeBIC', 'rfe15']

for ftmethod in methods:

    outfile = 'output/02_ft-selection/temp/{}_scikit_features.tsv'.format(ftmethod)
    with open(outfile, 'w') as out:
        out.write('Feature_Set_ID\tTCGA_Projects\tFeatures\n')

        for cancer in cancer_list:
            ###
            # setup objects
            ###
            ftid = 'skgrid_{}_{}'.format(cancer, ftmethod)
            file  = 'output/02_ft-selection/{}/{}_{}--ALL120920.tsv'.format(cancer, cancer, ftmethod)
            out.write(ftid + '\t' + '["{}"]'.format(cancer) + '\t' + get_ft_str(file) + '\n')
