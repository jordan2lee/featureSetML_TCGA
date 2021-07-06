#!usr/bin/python
import pandas as pd
import json
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--f_importance", help ="file with importance scores", required=True, type=str)
    parser.add_argument("-c", "--conversion_dictionary", help ="source dictionary of model matches", required=True, type=str)
    parser.add_argument("-o", "--outdir", help ="output dir with top prefix", required=True, type=str)
    return parser.parse_args()

args = get_arguments()
f_imp = args.f_importance
f_convert = args.conversion_dictionary
outdir = args.outdir

def get_team(name):
    '''get name ex CF, AKLIMATE, jadbio, subSCOPE, skgrip'''
    team = name.strip().split('_')[0].upper()
    # correct for trailing info for subscope
    if team.startswith('SUBSCOPE'):
        team = 'SUBSCOPE'
    cohort =name.strip().split('_')[-1]
    return(cohort+'_'+team)

# Read in Converter and importance matrix
with open(f_convert, 'r') as fh:
    dict1 = json.load(fh)
df = pd.read_csv(f_imp, sep ='\t')

# Pull out all models
main_list =[]
alt_list = []
for k,v in dict1.items():
    main_list.append(k)
    alt_list.append(v)

for i in range(0, len(main_list)):
    main_model_name = main_list[i]

    # Skip jadbio because have source files for it
    if 'jadbio' in main_model_name:
        continue
    else:
        # Else
        alt_name = alt_list[i]
        dict_importance = df[df['feature_importance_ID']==alt_name].reset_index(drop=True)['json_object'][0]
        dict_importance = json.loads(dict_importance)

        # Populate df and save
        features = []
        importance = []
        try:
            for k,v in dict_importance['feature_importance_scores'].items():
                features.append(k)
                importance.append(v)
            assert len(features)==len(importance)
        except:
            for k,v in dict_importance['feature_importance_score'].items():
                features.append(k)
                importance.append(v)
            assert len(features)==len(importance)
        out = pd.DataFrame(list(zip(features, importance)), columns = ['features', 'importance'])
        out.to_csv(outdir + '/TOP_{}.tsv'.format(get_team(main_model_name)), sep='\t', index=False)
