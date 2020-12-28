#!/bin/bash

tumor=$1

echo "####################"
echo "Now running "${tumor}
echo "####################"

# # Feature Selection: FBED (thres 0.05, conditional_test testIndMultinom, kmax 2, crit eBIC)
# 	#input: original ft matirx, output: feature selection files in 02_ft-selection/
# echo "Feature selection: FBED"
# Rscript FUNCTION-ft-sel_fbed.R \
# 	--file ../raw/${tumor}_v8_20200203.tsv \
# 	--class ${tumor} \
# 	--output_path ../output/02_ft-selection/${tumor}/

# use this for latest run if not KIRCKICH or LIHCCHOL 12/11/20
# # Feature selection of all data platforms as input
# echo "Feature selection: FBED - all platforms"
# Rscript FUNCTION-ft-sel_fbed_ALLDTYPES.R \
# 	--file ../raw/${tumor}_v8_20200203.tsv \
# 	--class ${tumor} \
# 	--output_path ../output/02_ft-selection/${tumor}/


# use this for latest run if KIRCKICH or LIHCCHOL 12/11/20
echo "Feature selection: FBED - all platforms"
Rscript FUNCTION-ft-sel_fbed_ALLDTYPES.R \
	--file ../output/00_mi_rm/${tumor}_v8_20200203_nomir.tsv \
	--class ${tumor} \
	--output_path ../output/02_ft-selection/${tumor}/
