#!/usr/bin/Rscript

# Packages
.packages = c("BiocManager","remotes", "tidyr", "dplyr", "gridtext", "magick", "stringr", "data.table", "testit", "knitr", "ggplot2", "circlize", "docstring", "argparse", "ComplexUpset")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst], repos="http://cran.us.r-project.org")


# these will need a and Yes prompts
BiocManager::install("ComplexHeatmap")
remotes::install_github("jokergoo/ComplexHeatmap")
