#!/usr/bin/bash
# created on 12/14/20 to combine all output files and format it

mkdir output/02_ft-selection/temp/

python scripts/format_output.py

#combine into one file
cat output/02_ft-selection/temp/rfe15_scikit_features.tsv > output/02_ft-selection/skgrid_features_121420.tsv
cat output/02_ft-selection/temp/fbedeBIC_scikit_features.tsv | grep -v "^Feature_Set_ID" >> output/02_ft-selection/skgrid_features_121420.tsv


# cleanup
rm -r output/02_ft-selection/temp/
