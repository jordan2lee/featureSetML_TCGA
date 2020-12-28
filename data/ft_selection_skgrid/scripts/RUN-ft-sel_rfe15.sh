#!/bin/bash

tumor=${1}

echo "####################"
echo "Now running "${tumor}
echo "####################"

# # Feature Selection: Recursive feature elimination 15 features per platform
# 	#input: original ft matirx, output: feature selection files in 02_ft-selection/
# echo "Feature selection: RFE15"
# python3.7 FUNCTION-ft-sel_rfe15.py \
#    --input_path ../raw/ \
#    --input_matrix ${tumor}_v8_20200203.tsv \
#    --output_path ../output/02_ft-selection/${tumor}/ \
#    --tumor ${tumor}

# # use this for latest run if not KIRCKICH or LIHCCHOL 12/11/20
# # Feature Selection: Recursive feature elimination 15 features for all platforms combined
# 	#input: original ft matirx, output: feature selection files in 02_ft-selection/
# echo "Feature selection: RFE15 all platforms"
# python3.7 FUNCTION-ft-sel_rfe15_ALLDTYPES.py \
#    --input_path ../raw/ \
#    --input_matrix ${tumor}_v8_20200203.tsv \
#    --output_path ../output/02_ft-selection/${tumor}/ \
#    --tumor ${tumor}


# use this for latest run if KIRCKICH or LIHCCHOL 12/11/20
# Feature Selection: Recursive feature elimination 15 features for all platforms combined
	#input: original ft matirx, output: feature selection files in 02_ft-selection/
echo "Feature selection: RFE15 all platforms"
python3.7 FUNCTION-ft-sel_rfe15_ALLDTYPES.py \
   --input_path ../raw/ \
   --input_matrix ../output/00_mi_rm/${tumor}_v8_20200203_nomir.tsv \
   --output_path ../output/02_ft-selection/${tumor}/ \
   --tumor ${tumor}
