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
## 1. Run Feature Selection Algorithms

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

## 2. Exact Feature Name Overlap Upset Plot (Base Upset Plot)

Last updated: 1/8/21

Purpose: create an upset plot that shows how much overlap there is between the feature sets of the best models

Methods: pull the best model per team (based on mean `overall weighted F1` score). pull corresponding feature set for each model. look at raw overlap of features across teams

```
# Preprocessing and Upset Plot - in script/ it calls get_fts.py and upset_exactMatch.R
bash RUN_UpsetPlot_Exact.sh
```

Outputs two files in `data/figure_panel_a/`:

+ `best_models_<cancer>.tsv`
+ `upsetplot_<cancer>.pdf`

## 3. Exact Feature Overlap Heatmap - Clustering features

Last updated: 1/20/21

Purpose: create heatmap of all features of best models and cluster based on original molecular tarball signatures

Methods: select a cancer cohort. pull the best model per team (based on mean `overall weighted F1` score). Pull corresponding feature set for each model. Scale if needed (depends on data platform). Cluster features and show on heatmap.

*Note: assumes you ran scripts from section* `2. Exact Feature Name Overlap Upset Plot (Base Upset Plot)`

```
# Preprocessing - in scripts/ it calls get_fts.py, pull_grp_best.py, extract_hallmarks
bash RUN_Heatmap_Exact.sh

# Create heatmaps
## run scripts/heatmap.Rmd
```
