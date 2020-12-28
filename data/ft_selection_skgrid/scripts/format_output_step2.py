#!/usr/bin/env python3

file1 = 'output/02_ft-selection/temp/fbedeBIC_scikit_features.tsv'

file2 = 'output/02_ft-selection/temp/rfe15_scikit_features.tsv'


with open('output/02_ft-selection/skgrid_features_121420.tsv', 'w') as out, open(file1, 'r') as fh1, open(file2, 'r') as fh2:
    for line in fh1:
        out.write(line)
        print(line)
    ct = 0
    for line in fh2:
        if ct ==0:
            ct ==1
        else:
            out.write(line)
            print(line)
