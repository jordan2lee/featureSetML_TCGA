#!/usr/bin/bash

# run after RUN_ft_select_skgrid.sh

# 3. To Synapse format for fbed
python scripts/synapse.py \
    --outfile data/ft_selection_skgrid/skgrid_featuresets_010520.tsv
