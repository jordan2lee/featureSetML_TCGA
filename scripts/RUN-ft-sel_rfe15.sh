#!/bin/bash

tumor=${1}

echo "####################"
echo "Now running "${tumor}
echo "####################"
echo 'RFE15'

# 1. Feature Selection: Recursive feature elimination 15 features per platform
	#input: original ft matirx, output: feature selection files in 02_ft-selection/
echo "Feature selection: RFE15"
python3.7 scripts/FUNCTION-ft-sel_rfe15.py \
   --input_path src/ \
   --input_matrix ${tumor}_v8_20200203.tsv \
   --output_path ../data/ft_selection_skgrid/${tumor}/ \
   --tumor ${tumor}


# Feature Selection: Recursive feature elimination 15 features for all platforms combined
# 2A. Use for ONLY for KIRCKICH or LIHCCHOL
if [[ "${tumor}" = "KIRCKICH" ]] || [[ "${tumor}" = "LIHCCHOL" ]]; then
    echo "Feature selection: (for KIRCKICH/LIHCCHOL) RFE15 all platforms"
    python3.7 scripts/FUNCTION-ft-sel_rfe15_ALLDTYPES.py \
       --input_path data/ \
       --input_matrix manual_remove_mi/${tumor}_v8_20200203_nomir.tsv \
       --output_path ft_selection_skgrid/${tumor}/ \
       --tumor ${tumor}
# 2B. Use for all other cancer
else
    echo "Feature selection: (for any non- KIRCKICH/LIHCCHOL) RFE15 all platforms"
    python3.7 scripts/FUNCTION-ft-sel_rfe15_ALLDTYPES.py \
       --input_path src/ \
       --input_matrix ${tumor}_v8_20200203.tsv \
       --output_path ../data/ft_selection_skgrid/${tumor}/ \
       --tumor ${tumor}
fi
