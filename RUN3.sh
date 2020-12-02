#!/usr/bin/bash

# 0. Format feature lists of groups (best performing model)
python scripts/get_fts.py
echo 'completed ft list formatting'

# 1. Clean
echo '-- run: notebooks/wip/step1_clean.ipynb --'

# 2. Correlations
echo '-- run: scripts/wip/corr_matrix.Rmd --'
