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

Rscript scripts/upset_exactMatch.R \
    -c BRCA -m intersect --infile data/exact_match/best_models_BRCA.tsv \
    --outdir data/exact_match --outname upsetPlot_intersect_exact.pdf \
    --headers Gnosis,CloudForest,AKLIMATE,SubSCOPE,SciKitGrid
echo 'completed upset plot - mode intersect'

# 2. Correlations Exact Feature Match - Overlaps
python scripts/corr.py \
    -c BRCA --threshold 0.5 \
    -n1 gnosis -g1 'gnosis_1_BRCA' \
    -n2 cloudforest -g2 'CF|All_Top 100_BRCA' \
    --master src/tarball/BRCA_v8_20200203.tsv \
    --features data/exact_match/best_models_BRCA.tsv \
    -o data/exact_match/corr
echo 'completed gnosis vs cloudforest'

python scripts/corr.py \
    -c BRCA --threshold 0.5 \
    -n1 gnosis -g1 'gnosis_1_BRCA' \
    -n2 aklimate -g2 'AKLIMATE_BRCA_reduced_model_1000_feature_set_BRCA' \
    --master src/tarball/BRCA_v8_20200203.tsv \
    --features data/exact_match/best_models_BRCA.tsv \
    -o data/exact_match/corr
echo 'completed gnosis vs aklimate'

python scripts/corr.py \
    -c BRCA --threshold 0.5 \
    -n1 gnosis -g1 'gnosis_1_BRCA' \
    -n2 subscope -g2 'nn_jg_2020-03-20_top1kfreq:BRCA_BRCA' \
    --master src/tarball/BRCA_v8_20200203.tsv \
    --features data/exact_match/best_models_BRCA.tsv \
    -o data/exact_match/corr
echo 'completed gnosis vs subscope'

python scripts/corr.py \
    -c BRCA --threshold 0.5 \
    -n1 gnosis -g1 'gnosis_1_BRCA' \
    -n2 scikitgrid -g2 'fbedeBIC_BRCA' \
    --master src/tarball/BRCA_v8_20200203.tsv \
    --features data/exact_match/best_models_BRCA.tsv \
    -o data/exact_match/corr
echo 'completed gnosis vs sci-kit grid'

python scripts/corr.py \
    -c BRCA --threshold 0.5 \
    -n1 cloudforest -g1 'CF|All_Top 100_BRCA' \
    -n2 aklimate -g2 'AKLIMATE_BRCA_reduced_model_1000_feature_set_BRCA' \
    --master src/tarball/BRCA_v8_20200203.tsv \
    --features data/exact_match/best_models_BRCA.tsv \
    -o data/exact_match/corr
echo 'completed cloudforest vs aklimate'

python scripts/corr.py \
    -c BRCA --threshold 0.5 \
    -n1 cloudforest -g1 'CF|All_Top 100_BRCA' \
    -n2 subscope -g2 'nn_jg_2020-03-20_top1kfreq:BRCA_BRCA' \
    --master src/tarball/BRCA_v8_20200203.tsv \
    --features data/exact_match/best_models_BRCA.tsv \
    -o data/exact_match/corr
echo 'completed cloudforest vs subscope'

python scripts/corr.py \
    -c BRCA --threshold 0.5 \
    -n1 cloudforest -g1 'CF|All_Top 100_BRCA' \
    -n2 scikitgrid -g2 'fbedeBIC_BRCA' \
    --master src/tarball/BRCA_v8_20200203.tsv \
    --features data/exact_match/best_models_BRCA.tsv \
    -o data/exact_match/corr
echo 'completed cloudforest vs sci-kit grid'

python scripts/corr.py \
    -c BRCA --threshold 0.5 \
    -n1 aklimate -g1 'AKLIMATE_BRCA_reduced_model_1000_feature_set_BRCA' \
    -n2 subscope -g2 'nn_jg_2020-03-20_top1kfreq:BRCA_BRCA' \
    --master src/tarball/BRCA_v8_20200203.tsv \
    --features data/exact_match/best_models_BRCA.tsv \
    -o data/exact_match/corr
echo 'completed aklimate vs subscope'

python scripts/corr.py \
    -c BRCA --threshold 0.5 \
    -n1 aklimate -g1 'AKLIMATE_BRCA_reduced_model_1000_feature_set_BRCA' \
    -n2 scikitgrid -g2 'fbedeBIC_BRCA' \
    --master src/tarball/BRCA_v8_20200203.tsv \
    --features data/exact_match/best_models_BRCA.tsv \
    -o data/exact_match/corr
echo 'completed aklimate vs sci-kit grid'

python scripts/corr.py \
    -c BRCA --threshold 0.5 \
    -n1 subscope -g1 'nn_jg_2020-03-20_top1kfreq:BRCA_BRCA' \
    -n2 scikitgrid -g2 'fbedeBIC_BRCA' \
    --master src/tarball/BRCA_v8_20200203.tsv \
    --features data/exact_match/best_models_BRCA.tsv \
    -o data/exact_match/corr
echo 'completed subscope vs sci-kit grid'
