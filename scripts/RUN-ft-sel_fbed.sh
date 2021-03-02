#!/bin/bash

# tumor=$1

timestamp() {
  date
}
timestamp

# echo "####################"
# echo "Now running "${tumor}
# echo "####################"
# echo "FBED (thres 0.05, conditional_test testIndMultinom, kmax 2, crit eBIC)"

# # Feature Selection: FBED separate for each platform
# echo "Feature selection: FBED separate per platform"
# Rscript scripts/FUNCTION-ft-sel_fbed.R \
# 	--file src/tarball/${tumor}_v11_20210220.tsv \
# 	--class ${tumor} \
# 	--output_path data/ft_selection_skgrid/${tumor}/

declare -a cancer_array=(CESC)
declare -a p_array=(0.5 0.001)
declare -a method_array=(eBIC)
declare -a backward_array=(FALSE)

# declare -a cancer_array=(CESC COADREAD ESCC GEA HNSC KIRCKICH KIRP)
# declare -a method_array=(eBIC LR)
for back in ${backward_array[@]}; do
  for p in ${p_array[@]}; do
    for tumor in ${cancer_array[@]}; do
      mkdir data/ft_selection_skgrid/${tumor} # TODO dev
      for m in ${method_array[@]}; do
      	# Feature selection: FBED for all platforms combined
      	# 2B. Use for all other cancers
      	timestamp
      	echo "Feature selection: FBED all platforms"
      	echo 'cancer:'${tumor} ' method:' ${m}
        echo 'extra' ${back} ${p}
      	Rscript scripts/FUNCTION-ft-sel_fbed_ALLDTYPES.R \
      		--file src/tarball/${tumor}_v11_20210220.tsv \
      		--cancer ${tumor} \
      		--conditional_test testIndMultinom \
      		--method ${m} \
      		-k 0 \
      		--thres ${p} \
      		--backphase FALSE \
      		--output_path data/ft_selection_skgrid/${tumor}/
      done
    done
  done
done
