# TODO update this script and others to finalized versions

# Set up workspace

```
. venv/bin/activate
bash init.sh
```
# Exact Feature Name Overlap Upset Plot

Last updated: 12/18/20

Purpose: create an upset plot that shows how much overlap there is between the feature sets of the best models

Methods: pull the best model per team (currently based on mean overall weighted F1 score). pull corresponding feature set for each model. look at overlap of features across teams


```

```


# Analysis

```
bash RUN.sh
```

Additional:

+ Distribution of correlation coefficients `scripts/corr_distributions.Rmd`


# WIP

Then run `notebooks/gene_overlaps` (TODO update hardcoded paths to files)
