#!/usr/bin/Rscript

# 11/16/20 JordanL
# Purpose: create upset plots of the overlap of features between groups (best model)

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-c", "--cancer", type='character', help='cancer type')
parser$add_argument('-i', '--indir', type='character', help='input directory')
parser$add_argument('-o', '--outdir', type='character', help='output directory')
parser$add_argument('-on', '--outname', type='character', help='output figure name')
parser$add_argument('-m', '--mode', type='character', help='upset plot grouping mode (ex. distinct)')
args <- parser$parse_args()


# 1. Get group feature sets
file = paste(args$indir, args$cancer, '_featurelist.tsv', sep='')
df <- read.csv(file, header=TRUE, sep='\t') %>% as.data.frame()

cforest <- df$BRCA.CF.All_Top.100[!df$BRCA.CF.All_Top.100 %in% ''] %>% as.vector()
nn <- df$nn_jg_2020.03.20_top1kfreq.BRCA[!df$nn_jg_2020.03.20_top1kfreq.BRCA %in% ''] %>% as.vector()
ak <- df$AKLIMATE_BRCA_reduced_model_1000_feature_set[!df$AKLIMATE_BRCA_reduced_model_1000_feature_set %in% ''] %>% as.vector()
gn <- df$gnosis.BRCA.1[!df$gnosis.BRCA.1 %in% ''] %>% as.vector()
skit <- df$fbedeBIC_BRCA[!df$fbedeBIC_BRCA %in% ''] %>% as.vector()

# 2. Save upset plot
setwd(args$outdir)
pdf(args$outname, onefile=TRUE)

lt = list(CloudForest = cforest,
          SubSCOPE = nn,
          AKLIMATE = ak,
          JADBIO=gn,
          Scikit_grid=skit
)

m= make_comb_mat(
  list_to_matrix(lt),
  mode=args$mode
)

#m<-m[comb_degree(m) >= 2] #only show two or more overlaps
cs = comb_size(m)

ht = UpSet(
  m,
  pt_size = unit(4, "mm"),
  lwd = 3,
  top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),
  comb_col = c('black','orange','orange', "orange")[comb_degree(m)] # overlap =2 groups yellow
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
