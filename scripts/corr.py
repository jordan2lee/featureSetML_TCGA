#!/usr/bin/env python
import pandas as pd
import json
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import argparse


def get_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", "--cancer", help ="cancer", required=True, type=str)
    parser.add_argument("-f", "--features", help ="path to feature file", required=True, type=str)
    parser.add_argument("-m", "--master", help ="raw tarball main file", required=True, type=str)
    parser.add_argument("-o", "--outdir", help ="outdir", required=True, type=str)
    parser.add_argument("-g1", "--group1", help ="group 1", required=True, type=str)
    parser.add_argument("-n1", "--name1", help ="group 1 for output file name", required=True, type=str)
    parser.add_argument("-g2", "--group2", help ="group 2", required=True, type=str)
    parser.add_argument("-n2", "--name2", help ="group 2 for output file name", required=True, type=str)
    return parser.parse_args()
args = get_arguments()


# Features for each model
ft_df = pd.read_csv(args.features, sep = '\t')

# Raw tarball - samples x [labels, feature values]
master_df = pd.read_csv(args.master, sep='\t', index_col=0)

# save venn diagram
venn2([set(ft_df[args.group1].dropna()), set(ft_df[args.group2].dropna())], set_labels=[args.group1,args.group2])
plt.savefig(args.outdir + '/venn_{}_{}.png'.format(args.name1, args.name2) )

# Get feature lists
grp1_fts = ft_df[args.group1].dropna()
grp2_fts = ft_df[args.group2].dropna()

# Dataframe prep
cols = ['corr', 'pval', 'gene1', 'gene2', 'group1', 'group2']
gene1 = []
gene2 = []
corr = []
pval = []
group1 = []
group2 = []

# calculate corr
for g1 in grp1_fts:
    for g2 in grp2_fts:
        # run for all - including if same feature
        # correlation matrix (coor coefficients)
        x = master_df[g1]
        y = master_df[g2]
        c, p = scipy.stats.spearmanr(x, y)

        gene1.append(g1)
        gene2.append(g2)
        corr.append(c)
        pval.append(p)
        group1.append(args.group1)
        group2.append(args.group2)

# save correlation result values
res = pd.DataFrame(list(zip( corr, pval,gene1, gene2, group1, group2)), columns = cols)
res = res.reindex(res['corr'].abs().sort_values(ascending=False).index).reset_index(drop=True) # sort by abs corr
res.to_csv(args.outdir + '/corrMatrix_' + args.name1 + '_' + args.name2 + '.tsv', sep='\t',index=False)

# write log file
with open(args.outdir + '/log_{}_{}.txt'.format(args.name1, args.name2), 'w') as out:
    n_correlated = len(res[res['corr']> 0.5]) + len(res[res['corr']< -0.5])
    out.write(' '.join(['#', args.group1, args.group2]) + '\n')
    out.write(str(n_correlated) + ' correlated genes pairs (' + str((n_correlated/res.shape[0] )*100) + '%)\n')
    out.write('excludes instances where both groups picked same gene\n')
