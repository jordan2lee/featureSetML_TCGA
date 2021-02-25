#!/bin/bash

tumor=$1

echo "####################"
echo "Now running "${tumor}
echo "####################"
echo "FBED (thres 0.05, conditional_test testIndMultinom, kmax 2, crit eBIC)"

# # Feature Selection: FBED separate for each platform
# echo "Feature selection: FBED separate per platform"
# Rscript scripts/FUNCTION-ft-sel_fbed.R \
# 	--file src/tarball/${tumor}_v11_20210220.tsv \
# 	--class ${tumor} \
# 	--output_path data/ft_selection_skgrid/${tumor}/


# Feature selection: FBED for all platforms combined
# 2B. Use for all other cancers
echo "Feature selection: FBED all platforms"
Rscript scripts/FUNCTION-ft-sel_fbed_ALLDTYPES.R \
	--file src/tarball/${tumor}_v11_20210220.tsv \
	--cancer ${tumor} \
	--conditional_test testIndMultinom \
	--method eBIC \
	-k 2 \
	--thres 0.05 \
	--backphase TRUE \
	--output_path data/ft_selection_skgrid/${tumor}/
