#!usr/bin/Rscript

# Purpose: Extract hallmark gene sets

# install.packages("msigdbr")
suppressPackageStartupMessages(library("msigdbr"))

# Pull Hallmark to gene
hallmark_mapping <- msigdbr(species = "Homo sapiens", category=c("H")) # grab Hallmarks for homosap
hallmark_mapping <- hallmark_mapping[,c('human_gene_symbol','gs_name')] # select
hallmark_mapping <- hallmark_mapping[order(hallmark_mapping$human_gene_symbol),] # sort
hallmark_mapping$gs_name<- gsub("HALLMARK_","",hallmark_mapping$gs_name) # remove prefix

write.table(hallmark_mapping, file='data/figure_panel_b/hallmarks.tsv', quote=FALSE, sep='\t', row.names = FALSE)
