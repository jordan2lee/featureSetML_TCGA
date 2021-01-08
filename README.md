# Feature Selection Analysis
## Setup Workspace

1. Install Python and R

2. Create virtual environment, install packages, and create repo folders

```
python3 -m venv venv
. venv/bin/activate

pip install -r requirements.txt
bash init.sh
```
## Run Feature Selection Algorithms

Last updated: 1/5/21

Purpose: perform feature selection

Methods: RFE and FBED. For each cancer cohort, each ran in two iterations - across all molecular features (gene expression, copy number variation, miRNA, methylation, mutation status) - or ran independently for each of the five molecular feature types then after feature selection reported all features in one file.

Rationale: likely the dynamic range of the five molecular data types will bias feature selection algorithms towards data types with larger numbers (example gene expression can be in the thousands vs mutation status at most can be 1.0).

```
# Run feature selection models
bash RUN_ft_select_skgrid.sh <cancer-cohort>

# Consolidate and format feature selection results
bash RUN_to_synapse.sh
```

## Exact Feature Name Overlap Upset Plot (Base Upset Plot)

Last updated: 1/8/20

Purpose: create an upset plot that shows how much overlap there is between the feature sets of the best models

Methods: pull the best model per team (based on mean `overall weighted F1` score). pull corresponding feature set for each model. look at raw overlap of features across teams

```
# Preprocessing and Upset Plot - in script/ it calls get_fts.py and upset_exactMatch.R
bash RUN_UpsetPlot_Exact.sh
```

Outputs two files in `data/figure_panel_a/`:

+ `best_models_<cancer>.tsv`
+ `upsetplot_<cancer>.pdf`

## Molecular Feature Overlap Heatmap

Last updated: 12/18/20

Purpose: create heatmap of all features of best models and cluster based on original molecular tarball signatures

Methods: select a cancer cohort. pull the best model per team (currently based on mean overall weighted F1 score). pull corresponding feature set for each model. scale if needed (depends on data platform). cluster features and show on heatmap

```
# Preprocessing - in scripts/ it calls get_fts.py and pull_grp_best.py
bash RUN_ht.sh

# Heatmaps - manually run
figures/heatmap_fts.Rmd
```

Note that figures generated in `.Rmd` file were then manually copied into slides for group presentation `notebooks/leejordan_fig4_presentation_121820.pdf`


## Cluster Based Exact Feature Name Overlap Upset Plot

Last updated: 12/18/20

Purpose: create an upset plot that shows how much overlap there is between the feature sets of the best models

Methods: pull the best model per team (currently based on mean overall weighted F1 score). pull corresponding feature set for each model. cluster and incorporates pairwise correlation calculations between features look at overlap of features across teams

```
# Preprocessing - in scripts/ it calls get_fts.py, pull_grp_best.py and step1_clean.py
bash RUN_uplot.sh

# Upset plots - manually run
figures/cluster_fts_<datatype>.Rmd
```

Note: Created a few secondary files to the `cluster_fts_<datatype>.Rmd` which includes these 2

+ `scripts/archive/cluster_fts_gexp-explore.Rmd` this script was built first to explore the options. so only one data type was analyzed (GEXP) to see if this route is worthwhile running for all 5 data types. Includes aspects included in the main `cluster_fts_<datatype>.Rmd` (optimize clustering structure, optimize number of clusters, different clustering methods) but also includes aspects not in `cluster_fts_<datatype>.Rmd` such as testing out cluster plots to showcase similar features in a component analysis space and potentially other options

+ `figures/archive/cluster_fts_gexp-explore2_PAM.Rmd` adds the list of PAM50 features to the "teams" list shown in upset plots. compares how many groups selected features that also belong within the PAM50 set

Note that figures generated in `.Rmd` file were then manually copied into slides for group presentation `notebooks/leejordan_fig4_upset_clusters_presentation_121120.pdf`
