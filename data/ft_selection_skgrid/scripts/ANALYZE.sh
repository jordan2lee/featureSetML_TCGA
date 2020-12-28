#!/bin/bash

tumor=${1}
################
# Feature selection
################
mkdir ../output/02_ft-selection/${tumor}

# # Run feature selection, method=fbed
./RUN-ft-sel_fbed.sh ${tumor}

# Run feature selection, method=rfe15
./RUN-ft-sel_rfe15.sh ${tumor}
