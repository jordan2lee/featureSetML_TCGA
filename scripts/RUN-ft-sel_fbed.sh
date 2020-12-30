#!/bin/bash

tumor=$1

echo "####################"
echo "Now running "${tumor}
echo "####################"
echo "FBED (thres 0.05, conditional_test testIndMultinom, kmax 2, crit eBIC)"

# Feature Selection: FBED separate for each platform
echo "Feature selection: FBED separate per platform"
Rscript scripts/FUNCTION-ft-sel_fbed.R \
	--file src/${tumor}_v8_20200203.tsv \
	--class ${tumor} \
	--output_path data/ft_selection_skgrid/${tumor}/


# Feature selection: FBED for all platforms combined
# use this for latest run if KIRCKICH or LIHCCHOL 12/11/20
if [[ "${tumor}" = "KIRCKICH" ]] || [[ "${tumor}" = "LIHCCHOL" ]]; then
	echo "Feature selection: (for KIRCKICH/LIHCCHOL) FBED all platforms"
	Rscript scripts/FUNCTION-ft-sel_fbed_ALLDTYPES.R \
		--file data/manual_remove_mi/${tumor}_v8_20200203_nomir.tsv \
		--class ${tumor} \
		--output_path data/ft_selection_skgrid/${tumor}/
else
	# 2B. Use for all other cancers
	echo "Feature selection: (for any non- KIRCKICH/LIHCCHOL) FBED all platforms"
	Rscript scripts/FUNCTION-ft-sel_fbed_ALLDTYPES.R \
		--file src/${tumor}_v8_20200203.tsv \
		--class ${tumor} \
		--output_path data/ft_selection_skgrid/${tumor}/
fi
