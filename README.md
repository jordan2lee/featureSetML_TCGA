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
## 2. Run Classifiers and Download Group Results

Run feature sets through Sci-kit Grid classifiers. Combine prediction results with the results from other team members and download the merged feature set file and classifier performance file into `src/`

## 3. Create team conversions

Create team conversion dictionaries for the labeling and matching of model name and feature set. This is because some teams reported differently the naming convention of these two.

Need to incorporate this script indo team overlap RUN_overlap.sh (needs to be ran after figure_panel_a/best_models.tsv are ran).

```
# Create files per team
python scripts/link_model_ftset.py

# Combine files into one merged file
python scripts/combine_links.py
```

## 3. [Exploratory] Prediction Performance Loss with Feature Set Size

Last updated: 2/23/21

Purpose: Explore how much prediction performance is loss with smaller feature set size.

This will help inform step 4 on the threshold cutoff for feature set size when determining the best model from each team.

```
# Create figure showing loss of prediction performance
notebooks/wip/Feature_set_gain_loss.ipynb
```

This outputs 3 file per team per cancer in `notebooks/wip/wip_figures/`

## 4. Team Overlap Analysis

Last updated: 4/2/21

Two parts:

1. Exact Feature Name Overlap Upset Plot (Base Upset Plot)

2. Exact Feature Overlap Heatmap - Clustering features

Purpose 4.1: create an upset plot that shows how much overlap there is between the feature sets of the best models

+ Methods: pull the best model per team (based on mean `overall weighted F1` score). pull corresponding feature set for each model. look at raw overlap of features across teams

Purpose 4.2: create heatmap of all features of best models and cluster based on original molecular tarball signatures

+ Methods: select a cancer cohort. pull the best model per team (based on mean `overall weighted F1` score). Pull corresponding feature set for each model. Scale if needed (depends on data platform). Cluster features and show on heatmap.

Cancers: All 26 cancer cohorts

**Implements a feature set threshold of <=100 features. Feature sets with larger sizes will not be considered for "top model"**

First, save a local file containing hallmark gene sets

```
# Run once.
# Does not need to be recalled for each cancer type
bash RUN_Extract_Hallmark_file.sh
```

Create upset plots and heatmaps

```
# Upset Plot + Heatmaps
bash RUN_overlap.sh > log.txt
```

Outputs two files in `data/figure_panel_a/` and `data/figure_panel_b/`:
