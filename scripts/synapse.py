#!/usr/bin/env python3
# Purpose convert file to group format for synapse

import argparse

########## hardcoded
team = 'skgrid'
list_ft_methods = ['fbedeBIC', 'rfe15']
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
##########

#####
# Functions
#####
def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-o", "--outfile", help ="output file name", required=True, type=str)
    return parser.parse_args()


def parse_combined_file(input_file, ft_method):
    with open(input_file, 'r') as fhc:
        header = True
        ftset = []
        for line in fhc:
            if header == False:
                # FtsetID
                fID = '_'.join([team, tumor, ''.join([ft_method, iteration] ) ])
                # cancer cohort
                cohort = '[\"'+ tumor + '\"]'
                # ft set
                ft = line.strip().split('\t')[1]
                ftset.append(ft)
            else:
                header = False
        # clean up list of fts
        ftset = str(ftset).replace("\'", "\"").replace(' ', '')
        # write output results
        out.write('\t'.join([fID, cohort, ftset]) + '\n')


def parse_per_file(file, ft_method):
    ftset = []
    with open(file, 'r') as fh:
        header = True
        for line in fh:
            if header == False:
                # FtsetID
                fID = '_'.join([team, tumor, ''.join([ft_method, iteration] ) ])
                # cancer cohort
                cohort = '[\"'+ tumor + '\"]'
                # ft set
                line = line.strip().split('\t')[1:]
                for s in line:
                    if s != "NA":
                        ftset.append(s)
            else:
                header = False
    # clean up list of fts
    ftset = str(ftset).replace("\'", "\"").replace(' ', '')
    # write output results
    out.write('\t'.join([fID, cohort, ftset]) + '\n')


####
# Main
####
args = get_arguments()
file_out = args.outfile

with open(file_out, 'w') as out:
    out.write('Feature_Set_ID\tTCGA_Projects\tFeatures\n')

    # For each feature selection model
    for ft_method in list_ft_methods:
        # For each cancer pull the feature list
        for tumor in cancer_list:
            # For the combined file
            iteration = 'combined'
            file = 'data/ft_selection_skgrid/{}/{}_{}--combined_platform.tsv'.format(tumor, tumor, ft_method)
            parse_combined_file(file, ft_method)
            # For the per platform file
            iteration = 'perplatform'
            file = 'data/ft_selection_skgrid/{}/{}_{}--per_platform.tsv'.format(tumor, tumor, ft_method)
            parse_per_file(file, ft_method)
