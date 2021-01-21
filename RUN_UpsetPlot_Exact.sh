#!/usr/bin/bash

tumor_cohort='BRCA'
#
# # 0. Format feature lists of groups (best performing model)
# python scripts/get_fts.py \
#     --tumor ${tumor_cohort} \
#     -m overall_weighted_f1 \
#     -f1 src/collected_features_matrix_20200722.tsv.gz \
#     -f2 src/feature_list_with_performance_with_subtype_names_20200828.tsv.gz \
#     --out data/figure_panel_a/best_models_${tumor_cohort}.tsv
# echo 'completed ft list formatting'

# 1. Exact Feature Match - Overlaps
# Create upset plots
# Note that --headers must match order of --infile headers
Rscript scripts/upset_exactMatch.R \
    -c ${tumor_cohort} -m distinct --infile data/figure_panel_a/best_models_${tumor_cohort}.tsv \
    --outdir data/figure_panel_a --outname upsetplot_${tumor_cohort}.pdf \
    --headers JADBIO,CForest,AKLIMATE,SubSCOPE,SKGrid
echo 'completed upset plot - mode distinct'
