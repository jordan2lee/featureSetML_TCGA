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

# Add col with datatype abrev
col_vals <- c()
for (ft in rownames(df_fts)){
  shorten <- unlist(strsplit(ft, ':'))[2]
  col_vals <- c(col_vals, shorten)
}
df_fts['Platform']<- col_vals


# Plot
setwd(args$outdir)
# pdf(args$outname, onefile=TRUE, width=10)
tiff(
  args$outname,
  width = 1000,
  height = 1200,
  res = 200,
  compression = "none"
)
upset(
  # Main plot
  data = df_fts,
  intersect = model_headers,
  mode = 'distinct',
  name='',
  width_ratio=0.3,
  height_ratio = 0.75,
  wrap=TRUE,
  guides = 'over',
  sort_intersections_by = 'degree',
  sort_intersections = 'ascending',
  # Set Size plot
  set_sizes=(
    upset_set_size(
      # Color set size plot
      geom=geom_bar(
          aes(fill=Platform, x=group),
          width=0.8,
      )
      ,
      position = 'right'
    ) +
    theme(axis.ticks.x=element_line()) +
    geom_text(aes(label=..count..), hjust = -0.25,  size=rel(3),stat='count') + # Bar counts
    theme(axis.text.x=element_text(size=rel(1.125))) +  # set size x axis
    expand_limits(y=args$max_ftsize) # set max x value
  ),
  # Intersection matrix
  encode_sets = FALSE,
  matrix=(
    intersection_matrix(geom=geom_point(size = 1.75)) # Circle size
  ),
  annotations = list(
        'Data Platform'=(
            ggplot(mapping=aes(fill=Platform))
            + geom_bar(stat='count', position='fill')
            + scale_y_continuous(labels=scales::percent_format())
  + scale_fill_manual(values=c(
      'CNVR'='#00688b', 'GEXP'='#FFA500',
      'METH'='#43CD80', 'MIR'='#FF7F00',
      'MUTA' = '#00BFFF'
            ))
            + ylab('Data Platform')
        )
  ),
) +
ggtitle(paste('Feature Overlap Between Top ', args$cancer, ' Models', sep = ''))
dev.off()
