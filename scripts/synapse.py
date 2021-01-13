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
                if ft_method == 'fbedeBIC':
                    # FtsetID
                    fID = '_'.join([team, tumor, '_'.join([ft_method, iteration] ) ])
                    # cancer cohort
                    cohort = '[\"'+ tumor + '\"]'
                    # ft set
                    ft = line.strip().split('\t')[1]
                    ftset.append(ft)
                elif ft_method == 'rfe15':
                    # FtsetID
                    fID = '_'.join([team, tumor, '_'.join([ft_method, iteration] ) ])
                    # cancer cohort
                    cohort = '[\"'+ tumor + '\"]'
                    # ft set
                    ft = line.strip()
                    ftset.append(ft)
                else:
                    print('error, found a ft selection method not accounted for ', ft_method)
            else:
                header = False
        # clean up list of fts
        ftset = str(ftset).replace("\'", "\"").replace(' ', '')
        # write output results
        out.write('\t'.join([fID, cohort, ftset]) + '\n')


def parse_per_file(file, ft_method):
    ALL = [] # contains all 5 data types
    GEXP = [] # filters for gexp
    CNVR = [] # fitlers for cnvr
    MI = [] # filters for mirna
    METH = [] #filters for meth
    MUTA = []
    with open(file, 'r') as fh:
        header = True
        for line in fh:
            if header == False:
                # FtsetID
                fID = '_'.join([team, tumor, '_'.join([ft_method, iteration] ) ])
                # cancer cohort
                cohort = '[\"'+ tumor + '\"]'
                # ft set
                line = line.strip().split('\t')[1:]
                for s in line:
                    if s != "NA" and s != "":
                        # Add to each list
                        ALL.append(s)
                        if s.startswith('B:MUTA:'):
                            MUTA.append(s)
                        elif s.startswith('I:CNVR::'):
                            CNVR.append(s)
                        elif s.startswith('N:GEXP'):
                            GEXP.append(s)
                        elif s.startswith('N:METH:'):
                            METH.append(s)
                        elif s.startswith('N:MIR::'):
                            MI.append(s)
            else:
                header = False
    # clean up list of fts + write output results
    ALL = str(ALL).replace("\'", "\"").replace(' ', '')
    out.write('\t'.join([''.join([fID, 'ALL']), cohort, ALL]) + '\n')

    if len(MUTA)>0:
        MUTA = str(MUTA).replace("\'", "\"").replace(' ', '')
        out.write('\t'.join([''.join([fID, 'MUTA']), cohort, MUTA]) + '\n')
    if len(CNVR)>0:
        CNVR = str(CNVR).replace("\'", "\"").replace(' ', '')
        out.write('\t'.join([''.join([fID, 'CNVR']), cohort, CNVR]) + '\n')
    if len(GEXP)>0:
        GEXP = str(GEXP).replace("\'", "\"").replace(' ', '')
        out.write('\t'.join([''.join([fID, 'GEXP']), cohort, GEXP]) + '\n')
    if len(METH)>0:
        METH = str(METH).replace("\'", "\"").replace(' ', '')
        out.write('\t'.join([''.join([fID, 'METH']), cohort, METH]) + '\n')
    if len(MI)>0:
        MI = str(MI).replace("\'", "\"").replace(' ', '')
        out.write('\t'.join([''.join([fID, 'MI']), cohort, MI]) + '\n')


####
# Main
####
args = get_arguments()
file_out = args.outfile

with open(file_out, 'w') as out:
    out.write('Feature_Set_ID\tTCGA_Projects\tFeatures\n')

    # For each feature selection model
    for ft_method in list_ft_methods:
        print(ft_method)
        # For each cancer pull the feature list
        for tumor in cancer_list:
            print(tumor)
            # For the combined file
            iteration = 'combined'
            file = 'data/ft_selection_skgrid/{}/{}_{}--combined_platform.tsv'.format(tumor, tumor, ft_method)
            parse_combined_file(file, ft_method)
            # For the per platform file
            iteration = 'perplatform'
            file = 'data/ft_selection_skgrid/{}/{}_{}--per_platform.tsv'.format(tumor, tumor, ft_method)
            parse_per_file(file, ft_method)
