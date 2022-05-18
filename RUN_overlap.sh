#!/usr/bin/bash
# Exact Feature Match - Overlaps

timestamp() {
  date
}
timestamp

# Which cohorts? 4 main or remaining 24? -- remaining 24 will have legends and no htmap, manually comment out these lines
declare -a StringArray=('BRCA' 'LGGGBM' 'COADREAD' 'SKCM')
# declare -a StringArray=('ACC' 'BLCA' 'CESC' 'ESCC' 'GEA' 'HNSC' 'KIRCKICH' 'KIRP' 'LIHCCHOL' 'LUAD' 'LUSC' 'MESO' 'OV' 'PAAD' 'PCPG' 'PRAD' 'SARC' 'TGCT' 'THCA' 'THYM' 'UCEC' 'UVM')

for tumor_cohort in ${StringArray[@]}; do
  echo $tumor_cohort
  # 1. Format feature lists of groups (best performing model)
  python scripts/get_fts.py \
      --tumor ${tumor_cohort} \
      --file_fts src/classifier_metrics_20220511/collected_features_matrix.tsv \
      --file_top src/classifier_metrics_20220511/top_performing_models_lte_100_features.tsv \
      --out data/figure_panel_a/best_models_${tumor_cohort}.tsv
  echo 'completed ft list formatting'

  # 2. Clean files and create 2 new files
  python scripts/pull_grp_best.py \
      --tumor ${tumor_cohort} \
      -f1 data/figure_panel_a/best_models_${tumor_cohort}.tsv \
      -f2 src/tarball/${tumor_cohort}_v12_20210228.tsv \
      --out1 data/figure_panel_b/${tumor_cohort}_fts_by_TEAM.tsv \
      --out2 data/figure_panel_b/${tumor_cohort}_fts_by_VALUE.tsv
  echo 'completed pulling group best and file cleaning'

# TWO CANCERS HAD ISSUE WITH THIS STEP AND SO WILL BE EXCLUDED FOR NOW ONWARD (UCEC LIHCCHOL)
  # 3. Clean up feature importance file
  python scripts/format_ft_imp.py \
    --f_importance src/classifier_metrics_20220511/feature_importance.tsv \
    --conversion_dictionary src/modelID_performance2importance_ALLCOHORTS.json \
    --outdir data/top_model_importances
  echo 'completed cleaning feature importance file'

  # 4. Create upset plots
  # Note that --headers must match order of --infile headers
  if [[ ${tumor_cohort} == 'SKCM' ]]; then
      msize='25'
  elif [[ ${tumor_cohort} == 'PAAD' ]]; then
      msize='35'
  elif [[ ${tumor_cohort} == 'SARC' ]]; then
      msize='75'
  elif [[ ${tumor_cohort} == 'TGCT'  || ${tumor_cohort} == 'THCA' ]]; then
      msize='65'
  else
      msize='110'
  fi

  # 4. Create upset and heatmap. once without and once with ft names displayed on ht
  Rscript scripts/figures.R \
      --show_features \
      --cancer ${tumor_cohort} \
      --min_n_team_overlap 2 \
      --max_ftsize ${msize} \
      --outdir_upset data/figure_panel_a \
      --outdir_ht ../figure_panel_b \
      --input_team_display ScikitGrid,JADBIO,CloudForest,SubSCOPE,AKLIMATE
  echo 'completed heatmap'

  # 5. Clean up workspace
  mv data/figure_panel_b/supplemental/*heatmap*.tiff data/figure_panel_b/heatmaps/
  echo ''
  echo ''
done
