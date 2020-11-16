#!/usr/bin/bash

# Format feature lists of groups (best performing model)
python scripts/get_fts.py

# Create upset plots
Rscript scripts/upset_exactMatch.R \
    -c BRCA -m distinct --indir data/correlation/CORR_ \
    --outdir figures --outname upsetPlot_distinct_exact.pdf

Rscript scripts/upset_exactMatch.R \
    -c BRCA -m intersect --indir data/correlation/CORR_ \
    --outdir figures --outname upsetPlot_intersect_exact.pdf #optional

# Correlations
print('Now run notebooks/01_correlations.ipynb for correlations')
