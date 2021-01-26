#!/usr/bin/Rscript

# 11/22/20 JordanL
# Purpose: create upset plot of the overlap of features between groups (best model) for Figure 4 panel A

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-c", "--cancer", type='character', help='cancer type')
parser$add_argument('-i', '--infile', type='character', help='input file')
parser$add_argument('-o', '--outdir', type='character', help='output directory')
parser$add_argument('-on', '--outname', type='character', help='output figure name')
parser$add_argument('-m', '--mode', type='character', help='upset plot grouping mode (ex. distinct)')
parser$add_argument('-headers', '--headers', type='character', help='Order matters, new headers')
args <- parser$parse_args()

# 1. Get group feature sets
df <- read.csv(args$infile, header=TRUE, sep='\t', row.names=1) %>% as.data.frame()

# Preprocess
df <- df[apply(df[,-1], 1, function(x) !all(x==0)),] # rm fts never found in any models
model_names <- as.vector(strsplit(args$headers, ",")[[1]])
colnames(df) <- model_names # rename
m <- df %>% as.matrix()


# Make Upset Plot
setwd(args$outdir)
pdf(args$outname, onefile=TRUE)
m <- make_comb_mat(m, mode = args$mode)
#m<-m[comb_degree(m) >= 2] #only show two or more overlaps
cs = comb_size(m)

ht = UpSet(
  m,
  pt_size = unit(3, "mm"),
  lwd = 3,
  top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),
  comb_col = c('black','orange','orange', "orange","orange")[comb_degree(m)], # overlap =2 groups orange
  set_order = c("AKLIMATE", "SubSCOPE", "CForest", "JADBIO","SKGrid")
)

ht = draw(ht)
co = column_order(ht)

nc = ncol(m)
decorate_annotation("Intersection\nsize", {
    grid.text(cs[co],
        x = 1:nc,
        y = unit(cs[co], "native") + unit(1, "mm"),
        gp = gpar(fontsize = 8),
        just = "bottom",
        default.units = "native")
})

dev.off()
