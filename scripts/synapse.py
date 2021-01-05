#!/usr/bin/env python3

# Purpose convert file to group format for synapse
########## hardcoded
team = 'skgrid'
ft_method = 'fbedeBIC'
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
# Function
#####
def parse_combined_file(input_file):
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


def parse_per_file(file):
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
file_out = 'data/ft_selection_skgrid/{}_{}.tsv'.format(ft_method, '010520')
with open(file_out, 'w') as out:
    out.write('Feature_Set_ID\tTCGA_Projects\tFeatures\n')

    # For each cancer pull the feature list
    for tumor in cancer_list:
        # for the combined file
        iteration = 'combined'
        file = 'data/ft_selection_skgrid/{}/{}_{}--combined_platform.tsv'.format(tumor, tumor, ft_method)
        parse_combined_file(file)


        # for the per platform file
        iteration = 'perplatform'
        file = 'data/ft_selection_skgrid/{}/{}_{}--per_platform.tsv'.format(tumor, tumor, ft_method)
        parse_per_file(file)
