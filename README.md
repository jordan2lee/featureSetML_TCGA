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

# Exact Feature Name Overlap Upset Plot (Base Upset Plot)

Last updated: 12/18/20

Purpose: create an upset plot that shows how much overlap there is between the feature sets of the best models

Methods: pull the best model per team (currently based on mean overall weighted F1 score). pull corresponding feature set for each model. look at raw overlap of features across teams

```
# Preprocessing and Upset Plot - in script/ it calls get_fts.py and upset_exactMatch.R
bash RUN_base_uplot.sh
```

Creates upset plot pdf files in `figures/upsetPlot_distinct_exact.pdf`

Note the default method used for overlap plotting is `distinct` but did explore using method `intersect`. Ultimately decided against using intersect because didn't match the narrative we were going for. If want to run `intersect` then just update input parameter fields when calling `scripts/upset_exactMatch.R`. The actual output figure of `intersect` will be moved to archive


# Cluster Based Exact Feature Name Overlap Upset Plot

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
