#!/bin/bash

tumor=$1

timestamp() {
  date
}
timestamp

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


declare -a pval_array=(0.001 0.01 0.05)
declare -a k_array=(0 1 2)

for p in ${pval_array[@]}; do
	for kmax in ${k_array[@]}; do
		# Feature selection: FBED for all platforms combined
		# 2B. Use for all other cancers
		timestamp
		echo "Feature selection: FBED all platforms"
		echo 'kmax:'${kmax} ' p threshold:' ${p}
		Rscript scripts/FUNCTION-ft-sel_fbed_ALLDTYPES.R \
			--file src/tarball/${tumor}_v11_20210220.tsv \
			--cancer ${tumor} \
			--conditional_test testIndMultinom \
			--method eBIC \
			-k ${kmax} \
			--thres ${p} \
			--backphase TRUE \
			--output_path data/ft_selection_skgrid/${tumor}/
	done
done
