# 5/10/21

Warning messages below when run on all 26 cancers. there are 3 cancers with errors. checked and all except one error is due to how teams report their features


# Message from SubSCOPE

the two features are ending up in the file because I concatenate CNV features across all cancer types, so even if the CNV feature is not in the UCEC file, it is in >1 of the other cancer types (quick check showed SKCM has both of these, for example). If the feature is not reported for a cancer type it gets filled with '0' instead. What seems to be happening is that at a high level the pan-cancer subtype classifier is finding that (the lack of) these two features is important for predicting a sample is part of the UCEC grouping vs other cancers, even though they would obviously not be important for subtype distinction within UCEC. This doesn't surprise me since the CNV model in subSCOPE showed an unfortunate trend for overfitting due to the large feature space in the input. You can keep/leave these two features out for UCEC, whichever works - or I can regenerate the UCEC and ENSEMBLE feature files excluding these two (earliest I can do this is Monday).  

# Messages below:

### Errors Pre-figure Building

```
LIHCCHOL
jadbio
jadbio_LIHCCHOL_MULTIDATATYPE_cumulative_feature_set24_LIHCCHOL  selected as best model for group
CF
CF_LIHCCHOL_All_Top_100_LIHCCHOL  selected as best model for group
AKLIMATE
AKLIMATE_GEXP_ONLY_LIHCCHOL_reduced_model_100_feature_set_LIHCCHOL  selected as best model for group
subSCOPE
subSCOPE-GEXP_2021-04-21_bootstrapfeatures_LIHCCHOL  selected as best model for group
skgrid
skgrid_LIHCCHOL_fbedeBIC_combined_LIHCCHOL  selected as best model for group
completed ft list formatting
Traceback (most recent call last):
  File "scripts/pull_grp_best.py", line 31, in <module>
    mat = tarball[['Labels']+ fts]
  File "/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/venv/lib/python3.7/site-packages/pandas/core/frame.py", line 2934, in __getitem__
    raise_missing=True)
  File "/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/venv/lib/python3.7/site-packages/pandas/core/indexing.py", line 1354, in _convert_to_indexer
    return self._get_listlike_indexer(obj, axis, **kwargs)[1]
  File "/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/venv/lib/python3.7/site-packages/pandas/core/indexing.py", line 1161, in _get_listlike_indexer
    raise_missing=raise_missing)
  File "/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/venv/lib/python3.7/site-packages/pandas/core/indexing.py", line 1252, in _validate_read_indexer
    raise KeyError("{} not in index".format(not_found))
KeyError: "['N:METH:cg20568322:LMOD1:TssD7:NA_N:METH:cg07671858:SYNPO2L:TssD3203:Island_N:GEXP::HNF4A:3172:', 'N:METH:cg16128701:RP11-363E6.3:TssD233:Island_N:GEXP::ILF3:3609:', 'N:METH:cg04703174:KHDRBS2:TssD610:Shore_N:GEXP::SLC2A4RG:56731:'] not in index"
completed pulling group best and file cleaning

UCEC
jadbio
jadbio_UCEC_MULTIDATATYPE_cumulative_feature_set20_UCEC  selected as best model for group
CF
CF_UCEC_All_Top_50_UCEC  selected as best model for group
AKLIMATE
AKLIMATE_UCEC_reduced_model_50_feature_set_UCEC  selected as best model for group
subSCOPE
subSCOPE-ENSEMBLE_2021-04-21_bootstrapfeatures_UCEC  selected as best model for group
skgrid
skgrid_UCEC_fbedeBIC_perplatformMUTA_UCEC  selected as best model for group
completed ft list formatting
Traceback (most recent call last):
  File "scripts/pull_grp_best.py", line 31, in <module>
    mat = tarball[['Labels']+ fts]
  File "/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/venv/lib/python3.7/site-packages/pandas/core/frame.py", line 2934, in __getitem__
    raise_missing=True)
  File "/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/venv/lib/python3.7/site-packages/pandas/core/indexing.py", line 1354, in _convert_to_indexer
    return self._get_listlike_indexer(obj, axis, **kwargs)[1]
  File "/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/venv/lib/python3.7/site-packages/pandas/core/indexing.py", line 1161, in _get_listlike_indexer
    raise_missing=raise_missing)
  File "/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/venv/lib/python3.7/site-packages/pandas/core/indexing.py", line 1252, in _validate_read_indexer
    raise KeyError("{} not in index".format(not_found))
KeyError: "['N:CNVR::TBXA2R:6915:', 'N:CNVR::ZBTB7A:51341:', 'N:CNVR::MRPL54:116541:', 'N:CNVR::MIR637:693222:', 'N:CNVR::MAP2K2:5605:', 'I:CNVR::ORAOV1:220064:', 'N:CNVR::GNG7:2788:', 'I:CNVR::GLDC:2731:', 'N:CNVR::CREB3L3:84699:'] not in index"
completed pulling group best and file cleaning
```

### Errors During Heatmap Building

```
KIRCKICH
Warning message:
Ignoring unknown parameters: inherit.blank
Error in if (all(is.finite(continuous_range_coord)) && diff(continuous_range_coord) <  :
  missing value where TRUE/FALSE needed
Calls: <Anonymous> ... expand_limits_continuous -> expand_limits_continuous_trans
Execution halted
Warning message:
Ignoring unknown parameters: inherit.blank
Error in if (all(is.finite(continuous_range_coord)) && diff(continuous_range_coord) <  :
  missing value where TRUE/FALSE needed
Calls: <Anonymous> ... expand_limits_continuous -> expand_limits_continuous_trans
Execution halted
mv: rename data/figure_panel_b/supplemental/*heatmap*.tiff to data/figure_panel_b/heatmaps/*heatmap*.tiff: No such file or directory
LIHCCHOL
Error in fread(paste("data/figure_panel_b/", args$cancer, "_fts_by_VALUE.tsv",  :
  File 'data/figure_panel_b/LIHCCHOL_fts_by_VALUE.tsv' does not exist or is non-readable. getwd()=='/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA'
Calls: %>% -> as.data.frame -> fread
Execution halted
Error in fread(paste("data/figure_panel_b/", args$cancer, "_fts_by_VALUE.tsv",  :
  File 'data/figure_panel_b/LIHCCHOL_fts_by_VALUE.tsv' does not exist or is non-readable. getwd()=='/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA'
Calls: %>% -> as.data.frame -> fread
Execution halted
mv: rename data/figure_panel_b/supplemental/*heatmap*.tiff to data/figure_panel_b/heatmaps/*heatmap*.tiff: No such file or directory
UCEC
Error in fread(paste("data/figure_panel_b/", args$cancer, "_fts_by_VALUE.tsv",  :
  File 'data/figure_panel_b/UCEC_fts_by_VALUE.tsv' does not exist or is non-readable. getwd()=='/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA'
Calls: %>% -> as.data.frame -> fread
Execution halted
Error in fread(paste("data/figure_panel_b/", args$cancer, "_fts_by_VALUE.tsv",  :
  File 'data/figure_panel_b/UCEC_fts_by_VALUE.tsv' does not exist or is non-readable. getwd()=='/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA'
Calls: %>% -> as.data.frame -> fread
Execution halted
mv: rename data/figure_panel_b/supplemental/*heatmap*.tiff to data/figure_panel_b/heatmaps/*heatmap*.tiff: No such file or directory
```
