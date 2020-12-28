#!/usr/bin/env python

# manually remove mirna
import pandas as pd
df = pd.read_csv('raw/KIRCKICH_v8_20200203.tsv', sep = '\t', index_col=0)
new = []

for c in df.columns:
    if c.startswith('N:MIR::'):
        skip = 'yes'
    else:
        new.append(c)

noMIR = df[new]
noMIR.to_csv('output/00_mi_rm/KIRCKICH_v8_20200203_nomir.tsv',sep='\t')


# manually remove mirna
import pandas as pd
df = pd.read_csv('raw/LIHCCHOL_v8_20200203.tsv', sep = '\t', index_col=0)
new = []

for c in df.columns:
    if c.startswith('N:MIR::'):
        skip = 'yes'
    else:
        new.append(c)

noMIR = df[new]
noMIR.to_csv('output/00_mi_rm/LIHCCHOL_v8_20200203_nomir.tsv',sep='\t')
