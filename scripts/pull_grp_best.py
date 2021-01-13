#!/usr/bin/env python3

##### Hardcoded
cancer = 'BRCA'
file = 'data/exact_match/best_models_BRCA.tsv'
file_tarball = 'src/tarball/{}_v9_20201029.tsv'.format(cancer)
out1 = 'data/heatmap/{}_fts_by_TEAM.tsv'.format(cancer)
out2 = 'data/heatmap/{}_fts_by_VALUE.tsv'.format(cancer)
#####

import pandas as pd

# open file and rm fts no model used
df = pd.read_csv(file, sep='\t', index_col=0)
df = df[(df.T != 0).any()]
df['Total'] = df.sum(axis=1) # add total col
df.to_csv(out1, sep='\t')

# Get values of original matrix for features of interest
tarball = pd.read_csv(file_tarball, sep = '\t', index_col = 0)
fts = list(df.index)
mat = tarball[['Labels']+ fts]
mat.to_csv(out2, sep='\t')
