#!/usr/bin/bash

echo '# Purpose: Create heatmap of feature set - for heatmap shows best model per team #'

# 0. Format feature lists of groups (best performing model)
python scripts/get_fts.py
echo 'completed ft list formatting'

# 1. Clean file and create 2 new files
python scripts/pull_grp_best.py
echo 'completed pulling group best and file cleaning'

# 2. Create heatmap
# then run scripts/heatmap_fts.Rmd
echo 'ACTION:'
echo '-- next step: create heatmap --'
echo '-- now manually run figures/heatmap_fts2.Rmd --'
