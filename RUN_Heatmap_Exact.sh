#!/usr/bin/bash

echo '# Purpose: Create heatmap of feature set - for heatmap shows best model per team #'

tumor_cohort='BRCA'

# 1. Clean file and create 2 new files
python scripts/pull_grp_best.py \
    --tumor ${tumor_cohort} \
    -f1 data/figure_panel_a/best_models_${tumor_cohort}.tsv \
    -f2 src/tarball/${tumor_cohort}_v8_20200203.tsv \
    --out1 data/figure_panel_b/${tumor_cohort}_fts_by_TEAM.tsv \
    --out2 data/figure_panel_b/${tumor_cohort}_fts_by_VALUE.tsv
echo 'completed pulling group best and file cleaning'

# 2. Create heatmap
# then run scripts/heatmap_fts.Rmd
echo 'ACTION:'
echo '-- next step: create heatmap --'
echo '-- now manually run figures/heatmap_fts2.Rmd --'
