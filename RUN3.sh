#!/usr/bin/bash

# # 0. Format feature lists of groups (best performing model)
# python scripts/get_fts.py
# echo 'completed ft list formatting'

# 1. Clean
python scripts/step1_clean.py \
    -f data/heatmap/BRCA_fts_by_VALUE.tsv \
    -o data/heatmap/BRCA_vals_meth.tsv \
    -d N:METH
python scripts/step1_clean.py \
    -f data/heatmap/BRCA_fts_by_VALUE.tsv \
    -o data/heatmap/BRCA_vals_gexp.tsv \
    -d N:GEXP
python scripts/step1_clean.py \
    -f data/heatmap/BRCA_fts_by_VALUE.tsv \
    -o data/heatmap/BRCA_vals_cnvr.tsv \
    -d I:CNVR::
python scripts/step1_clean.py \
    -f data/heatmap/BRCA_fts_by_VALUE.tsv \
    -o data/heatmap/BRCA_vals_muta.tsv \
    -d B:MUTA:
python scripts/step1_clean.py \
    -f data/heatmap/BRCA_fts_by_VALUE.tsv \
    -o data/heatmap/BRCA_vals_mir.tsv \
    -d N:MIR::
echo 'completed cleaning'

# 2. Correlations
echo '-- run: scripts/cluster_fts_<datatype>.Rmd --'
