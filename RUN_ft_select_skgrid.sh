#!/bin/bash

tumor=${1}
################
# Feature selection
################
mkdir data/ft_selection_skgrid/${tumor}

# # Run feature selection, method=fbed
# ./scripts/RUN-ft-sel_fbed.sh ${tumor}

# 2. Run feature selection, method=rfe15
./scripts/RUN-ft-sel_rfe15.sh ${tumor}
