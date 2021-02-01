#!/usr/bin/Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexUpset))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-c", "--cancer", type='character', help='cancer type')
parser$add_argument('-o', '--outdir', type='character', help='output directory')
parser$add_argument('-on', '--outname', type='character', help='output figure name')
parser$add_argument('-headers', '--model_headers', type='character', help='Order matters, new headers')
parser$add_argument('-m', '--max_ftsize', type='integer', help='output figure name')
args <- parser$parse_args()

# Read in file
df_fts <- fread(paste('data/figure_panel_b/', args$cancer, '_fts_by_TEAM.tsv', sep=''))%>% as.data.frame()

# Move index col and rm non model cols
row.names(df_fts) <- df_fts$featureID
df_fts <- df_fts[,!names(df_fts) %in% c('featureID', 'Total')]
model_headers <- unlist(as.vector(strsplit(args$model_headers, ",")))
colnames(df_fts) <- model_headers
# Transform matrix for input into upset
df_fts[model_headers] = df_fts[model_headers] == 1

# Plot
setwd(args$outdir)
pdf(args$outname, onefile=TRUE)
upset(
  # Main plot
  data = df_fts,
  intersect = model_headers,
  mode = 'distinct',
  name='Top Model',
  width_ratio=0.2,
  height_ratio = 0.3,
  wrap=TRUE,
  # Set Size plot
  set_sizes=(
    upset_set_size(position = 'right') +
    theme(axis.ticks.x=element_line()) +
    geom_text(aes(label=..count..), hjust = -0.25,  size=rel(2),stat='count') + # Bar counts
    theme(axis.text.x=element_text(size=rel(0.75))) +  # set size x axis
    expand_limits(y=args$max_ftsize) # set max x value
  ),
  # Intersection matrix
  encode_sets = FALSE,
  matrix=(
    intersection_matrix((geom=geom_point(size = 2)))
  )
  ) +
  ggtitle(paste('Feature Overlap Between Top ', args$cancer, ' Models', sep = ''))
dev.off()
