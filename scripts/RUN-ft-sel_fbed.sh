#!/bin/bash

tumor=$1

echo "####################"
echo "Now running "${tumor}
echo "####################"
echo "FBED (thres 0.05, conditional_test testIndMultinom, kmax 2, crit eBIC)"

# Feature Selection: FBED separate for each platform
echo "Feature selection: FBED separate per platform"
Rscript scripts/FUNCTION-ft-sel_fbed.R \
	--file src/tarball/${tumor}_v9_20201029.tsv \
	--class ${tumor} \
	--output_path data/ft_selection_skgrid/${tumor}/


# Feature selection: FBED for all platforms combined
# 2B. Use for all other cancers
echo "Feature selection: FBED all platforms"
Rscript scripts/FUNCTION-ft-sel_fbed_ALLDTYPES.R \
	--file src/tarball/${tumor}_v9_20201029.tsv \
	--class ${tumor} \
	--output_path data/ft_selection_skgrid/${tumor}/
