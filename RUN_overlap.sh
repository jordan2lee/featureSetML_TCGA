#!/usr/bin/bash
# Exact Feature Match - Overlaps

timestamp() {
  date
}
timestamp

declare -a StringArray=('ACC' 'BLCA' 'BRCA' 'CESC' 'COADREAD' 'ESCC' 'GEA' 'HNSC' 'KIRCKICH' 'KIRP' 'LGGGBM' 'LIHCCHOL' 'LUAD' 'LUSC' 'MESO' 'OV' 'PAAD' 'PCPG' 'PRAD' 'SARC' 'SKCM' 'TGCT' 'THCA' 'THYM' 'UCEC' 'UVM')

for tumor_cohort in ${StringArray[@]}; do
  echo $tumor_cohort

  # 1. Format feature lists of groups (best performing model)
  python scripts/get_fts.py \
      --tumor ${tumor_cohort} \
      -m overall_weighted_f1 \
      --filters 1000 \
      -f1 src/collected_features_matrix_20200722.tsv.gz \
      -f2 src/feature_list_with_performance_with_subtype_names_20200828.tsv.gz \
      --out data/figure_panel_a/best_models_${tumor_cohort}.tsv
  echo 'completed ft list formatting'

  # 2. Clean files and create 2 new files
  python scripts/pull_grp_best.py \
      --tumor ${tumor_cohort} \
      -f1 data/figure_panel_a/best_models_${tumor_cohort}.tsv \
      -f2 src/tarball/${tumor_cohort}_v9_20201029.tsv \
      --out1 data/figure_panel_b/${tumor_cohort}_fts_by_TEAM.tsv \
      --out2 data/figure_panel_b/${tumor_cohort}_fts_by_VALUE.tsv
  echo 'completed pulling group best and file cleaning'

  # 3. Create upset plots
  # Note that --headers must match order of --infile headers
  if [[ ${tumor_cohort} == 'LGGGBM' ]]; then
      msize='800'
  elif [[ ${tumor_cohort} == 'GEA' ]]; then
      msize='600'
  elif [[  ${tumor_cohort} == 'BRCA' ]]; then
      msize='1000'
  else
      msize='1000'
  fi

  Rscript scripts/upset.R \
      -c ${tumor_cohort} \
      --model_headers JADBIO,CForest,AKLIMATE,SubSCOPE,SKGrid \
      --max_ftsize ${msize} \
      --outdir data/figure_panel_a --outname upsetplot_${tumor_cohort}.pdf
  echo 'completed upset plot - mode distinct'

  # 4. Create heatmap
  Rscript scripts/heatmap_importance.R \
      --cancer ${tumor_cohort} \
      --min_n_team_overlap 2 \
      --supplemental data/figure_panel_b/supplemental/ \
      --outdir ../main/
  echo 'completed heatmap'

  # 5. Clean up workspace
  mv data/figure_panel_b/supplemental/*heatmap*.tiff data/figure_panel_b/heatmaps/

  echo ''
  echo ''
done
