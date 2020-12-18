#!/usr/bin/bash

# 0. Format feature lists of groups (best performing model)
python scripts/get_fts.py
echo 'completed ft list formatting'

# 1. Exact Feature Match - Overlaps
# Create upset plots
Rscript scripts/upset_exactMatch.R \
    -c BRCA -m distinct --infile data/exact_match/best_models_BRCA.tsv \
    --outdir data/exact_match --outname upsetPlot_distinct_exact.pdf \
    --headers Gnosis,CloudForest,AKLIMATE,SubSCOPE,SciKitGrid
echo 'completed upset plot - mode distinct'
