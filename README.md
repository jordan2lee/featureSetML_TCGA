# TODO update this script and others to finalized versions

# Set up workspace

```
. venv/bin/activate
bash init.sh
```
# Molecular Feature Overlap Heatmap

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

# TO CLEAN UP ORGANIZATION - Analysis

```
bash RUN.sh
```

Additional:

+ Distribution of correlation coefficients `scripts/corr_distributions.Rmd`


# TO CLEAN UP ORGANIZATION - WIP

Then run `notebooks/gene_overlaps` (TODO update hardcoded paths to files)
