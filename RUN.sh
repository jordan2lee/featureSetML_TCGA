#!/usr/bin/bash

# 0. Format feature lists of groups (best performing model)
python scripts/get_fts.py
echo 'completed ft list formatting'

# 1. Exact Feature Match - Overlaps
# Create upset plots
Rscript scripts/upset_exactMatch.R \
    -c BRCA -m distinct --infile data/exact_match/best_models_BRCA.tsv \
    --outdir data/exact_match --outname upsetPlot_distinct_exact.pdf \
    --headers Gnosis,CloudForest,AKLIMATE,SubScope,SciKitGrid
echo 'completed upset plot - mode distinct'

Rscript scripts/upset_exactMatch.R \
    -c BRCA -m intersect --infile data/exact_match/best_models_BRCA.tsv \
    --outdir data/exact_match --outname upsetPlot_intersect_exact.pdf \
    --headers Gnosis,CloudForest,AKLIMATE,SubScope,SciKitGrid
echo 'completed upset plot - mode intersect'



# # 2. Correlations Exact Feature Match - Overlaps
# # Options 'gnosis:BRCA:1',
# #         'BRCA:CF|All_Top 100',
# #         'AKLIMATE_BRCA_reduced_model_1000_feature_set',
# #         'nn_jg_2020-03-20_top1kfreq:BRCA',
# #         'fbedeBIC_BRCA'
# python scripts/corr.py \
#     -c BRCA \
#     -n1 gnosis -g1 'gnosis:BRCA:1' \
#     -n2 cloudforest -g2 'BRCA:CF|All_Top 100' \
#     --master src/tarball/BRCA_v8_20200203.tsv \
#     --features data/correlation/CORR_BRCA_featurelist.tsv \
#     -o data/correlation_exact_match
# echo 'completed gnosis vs cloudforest'
#
# python scripts/corr.py \
#     -c BRCA \
#     -n1 gnosis -g1 'gnosis:BRCA:1' \
#     -n2 aklimate -g2 'AKLIMATE_BRCA_reduced_model_1000_feature_set' \
#     --master src/tarball/BRCA_v8_20200203.tsv \
#     --features data/correlation/CORR_BRCA_featurelist.tsv \
#     -o data/correlation_exact_match
# echo 'completed gnosis vs aklimate'
#
# python scripts/corr.py \
#     -c BRCA \
#     -n1 gnosis -g1 'gnosis:BRCA:1' \
#     -n2 subscope -g2 'nn_jg_2020-03-20_top1kfreq:BRCA' \
#     --master src/tarball/BRCA_v8_20200203.tsv \
#     --features data/correlation/CORR_BRCA_featurelist.tsv \
#     -o data/correlation_exact_match
# echo 'completed gnosis vs subscope'
#
# python scripts/corr.py \
#     -c BRCA \
#     -n1 cloudforest -g1 'BRCA:CF|All_Top 100' \
#     -n2 aklimate -g2 'AKLIMATE_BRCA_reduced_model_1000_feature_set' \
#     --master src/tarball/BRCA_v8_20200203.tsv \
#     --features data/correlation/CORR_BRCA_featurelist.tsv \
#     -o data/correlation_exact_match
# echo 'completed cloudforest vs aklimate'
#
# python scripts/corr.py \
#     -c BRCA \
#     -n1 cloudforest -g1 'BRCA:CF|All_Top 100' \
#     -n2 subscope -g2 'nn_jg_2020-03-20_top1kfreq:BRCA' \
#     --master src/tarball/BRCA_v8_20200203.tsv \
#     --features data/correlation/CORR_BRCA_featurelist.tsv \
#     -o data/correlation_exact_match
# echo 'completed cloudforest vs subscope'
#
# # ... run corr.py for remaining combinations
#
# # Between JADBIO and SciKit_grid # KeyError: 'I:CNVR::CIRH1A:84916:'
# # Between CloudForest and SciKit_grid
# # Between AKLIMATE and SubSCOPE
# # Between AKLIMATE and SciKit_grid
# # Between SubSCOPE and SciKit_grid
