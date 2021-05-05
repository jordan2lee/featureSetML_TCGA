# feature_importance/

Contains team feature importance info. Each provided from the teams

```
2020v8_aggr_preds_filtcolumn.txt # Subscope
JADBio_TCGA_results_fig6_update_2020_01_allcohorts_MCastro.zip # JADBio
TMP_GEXP_v8_toplistsForPathwayAnalysis_20200203.RData ##Parsed GEXP lists with Top10, Top50, Top100, Top200, and Top500 genes by MCastro
```

# aklimate_feature_importance_scores_20200807.tar.gz

(syn22294880)

Contains AKLIMATE feature importance scores for v9 tarball. will be needed to update once final results come in

# tarball/

This dir contains the unzipped contents of raw molecular files (Theo's files). As of 12/29/20 this is for the V8 tarball. Note the paper will use V9 tarball that is identical to the V8 tarball except in that it has all miRNA features removed from `KIRCKICH` and  `LIHCCHOL` cancers.

# feature_list_with_performance_with_subtype_names_20200828.tsv.gz

(syn22337110)

Master feature list matrix of all group models. Created by Chris W. Downloaded from synapse on 11/22/20

Pre-miRNA fix

```
This is an update to the previous file (syn22322484).

This version adds datatype counts for the feature lists.

The columns are:
1 featureID
2 cohort
3 model
4 performance_metric
5 Mean
6 Median
7 Std
8 StdMed
9 SEM
10 Count
11 Max
12 Min
13 Sum
14 feature_list_method
15 feature_list_cohort
16 feature_list_size
17 subtypeID
18 subtype_size
19 TMP_name
20 Cluster_name
21 size_in_marker_paper
22 subtypes_defined_by_mRNA
23 subtypes_defined_by_CN
24 subtypes_defined_by_Mutations
25 subtypes_defined_by_Meth
26 subtypes_defined_by_miRNA
27 subtypes_defined_by_Histology
28 subtypes_defined_by_iCluster
29 subtypes_defined_by_COCA
30 subtypes_defined_by_Other
31 total_features
32 GEXP_features
33 CNVR_features
34 METH_features
35 MIR_features
36 MUTA_features
```

# collected_features_matrix_20200722.tsv.gz

(syn22271992)

Master feature list matrix of all group models. Created by Chris W. Downloaded from synapse on 11/22/20

Pre-miRNA fix

```
The current feature lists are collected in this file. Each row is a feature. Each column is a feature list. "0" or "1" indicates membership status in a feature list with "1" signifying membership. There is a special column, "total_number_of_lists", to record how many lists a feature participates in. Also, there are 3 special rows: "feature_list_method", "feature_list_cohort", "feature_list_size". The rows are sorted in descending order based on "total_number_of_lists" column.

The file records the data for:
24664 features
1272 feature lists
```
# brca_pam50_hits.tsv

BK created on 12/3/20. takes list of PAM50 genes. convert these into the TCGA naming system for GEXP features ONLY (so that have list of the PAM50 gexp features as named in Theo's tarball)

# classifier_metrics_20210430

(syn25653984)

Downloaded from Synapse on 5/3/21. This is the final combined results matrices. Created by Chris Wong. Contains model names with feature lists, feature importance scores, etc.

This should be the final set of matrices for the manuscript

# TMP_v12_20210228.tar.gz

(syn25007874)

Most recent tarball version. Downloaded on 5/3/21

# classifier_metrics_20210504

(syn25672419)

most recent results matrix from chris

Downloaded on 5/5/21
