#!/usr/bin/bash

# requires run of RUN.sh for step 00 as input into this RUN2.sh
#######
# Create heatmap of feature set - for heatmap shows one model per
#######
## 1. Get feature_list_cohort
# run heatmpa_topmodel.ipynb
# 2. then run heatmap_fts.Rmd


######
# Create heatmap of feature set - for heatmap shows best model per team
## 1. Get feature_list_cohort
python scripts/get_fts.py
echo 'completed ft list formatting'

## 2. Clean file and create 2 new files
# run pull_grp_best.ipynb

## 3. Create heatmap
# 2. then run heatmap_fts.Rmd
