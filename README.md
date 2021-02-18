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

## 3. Team Overlap Analysis

Last updated: 2/18/21

Two parts:

1. Exact Feature Name Overlap Upset Plot (Base Upset Plot)

2. Exact Feature Overlap Heatmap - Clustering features

Purpose 3A: create an upset plot that shows how much overlap there is between the feature sets of the best models

+ Methods: pull the best model per team (based on mean `overall weighted F1` score). pull corresponding feature set for each model. look at raw overlap of features across teams

Purpose 3B: create heatmap of all features of best models and cluster based on original molecular tarball signatures

+ Methods: select a cancer cohort. pull the best model per team (based on mean `overall weighted F1` score). Pull corresponding feature set for each model. Scale if needed (depends on data platform). Cluster features and show on heatmap.

Cancers: BRCA, LGGGBM, GEA

**Implements a feature set threshold of <1K features. Feature sets with larger sizes will not be considered for "top model"**

First, save a local file containing hallmark gene sets

```
# Run once.
# Does not need to be recalled for each cancer type
bash RUN_Extract_Hallmark_file.sh
```

Create upset plots and heatmaps

```
# Upset Plot
bash RUN_overlap.sh <cancer>
```

Outputs two files in `data/figure_panel_a/` and `data/figure_panel_b/`:
