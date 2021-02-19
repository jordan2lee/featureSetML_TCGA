#!usr/bin/Rscript

#####
# Packages and Functions
#####
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(gridtext))
suppressPackageStartupMessages(library(magick))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(testit))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(docstring))
suppressPackageStartupMessages(library(argparse))
source('func/cohort_config.R')
source('func/prep_for_heatmap.R')
source('func/draw_base_heatmap.R')
parser <- ArgumentParser()
parser$add_argument("-c", "--cancer", type='character', help='cancer cohort, all capitalized')
parser$add_argument('-op', '--outdir', type='character', help='output dir')
parser$add_argument('-os', '--supplemental', type='character', help='supplemental dir')
parser$add_argument('-min', '--min_n_team_overlap', type='character', help='min number of teams to show overlaps on all heatmaps')
args <- parser$parse_args()

######
# Hardcoded
######
file_imp_aklimate <- paste(
  'src/aklimate_feature_importance_scores_20200807/AKLIMATE_TEMPLATE_FULLY_TRAINED_',
  args$cancer,
  '_20200423_aklimate_ranked_feature_importance.tsv',
  sep=''
)
yes_scale <- c('N:METH', 'N:GEXP','N:MIR') # which fts to scale


####
# Read in files
####
# A. base files
df <- fread(
  paste(
    'data/figure_panel_b/',
    args$cancer,
    '_fts_by_VALUE.tsv',
    sep=''
  )
) %>% as.data.frame()

df_fts <- fread(
  paste(
    'data/figure_panel_b/',
    args$cancer,
    '_fts_by_TEAM.tsv',
    sep=''
  )
)%>% as.data.frame()

mappings <- fread(
  'data/figure_panel_b/hallmarks.tsv'
) %>% as.data.frame()

# B. Importance Score files
imp_aklimate <- fread(
  file_imp_aklimate
) %>% as.data.frame()

# C. Load Mauro's hallmark ranks for all cancer cohorts
load(file='src/mauro_files/Hallmark_nes_space_20210212.RData')

#####
# Set up
#####
# Define platform for hallmark heatmap
platform_of_interest <- get_platform_of_interest(args$cancer)

# Build list of data types present
platforms <- get_platforms_present(args$cancer)

#####
# Filter for ft rows that occur in >= 2 teams
#####
#first for df_fts
df_fts <- df_fts[df_fts['Total']>=args$min_n_team_overlap,]
rownames(df_fts) <- 1:nrow(df_fts) #reset row index names
cols_to_keep <- df_fts$featureID
cols_to_keep <- c(args$cancer, 'Labels', cols_to_keep)
#second for df
df <- subset(df, select=cols_to_keep)

##########
# Section 1: Heatmaps for each data type
##########
setwd(args$supplemental)
# Get models
models <- model2team(df_fts)
# Set up saving fig packet
for (prefix in platforms){
  # A. Base Heatmap
  results_list <- get_base_heatmap(
    prefix,
    args$cancer,
    models['JADBIO'],
    models['CForest'],
    models['AKLIMATE'],
    models['SubSCOPE'],
    models['SKGrid']
  )

  # B. If platform for Hallmarks Heatmap then save global variables - input for Section 2
  if (prefix == platform_of_interest){
    print('Using this as main platform of interest:')
    print(prefix)
    mat2 <- results_list[['results_matrix']]
    ftnames_order <- results_list[['results_ft_order']]
    subtype_annotation <- results_list[['subtype_annotation']]
    n_fts <- length(ftnames_order)
  }
}

######
# Section 2: Explore the hallmarks associated with GEXP features
######
# 1. Prep: look up if in hallmark mappings table
n_hallmarks <- c()
pooled_hallmarks <- c()

if (platform_of_interest == 'N:GEXP'){
  for (feature in ftnames_order){
    # 1. Preprocess - to gene symbol
    GENE <- unlist(strsplit(feature, '::'))[2]
    GENE <- unlist(strsplit(GENE, ':'))[1]
    # 2. Hallmark Mapping
    halls <- mappings[mappings$human_gene_symbol==GENE,]$gs_name %>% as.vector()
    # Multiple Hallmarks
    n_hallmarks <- c(n_hallmarks, length(halls))
    pooled_hallmarks <- c(pooled_hallmarks, halls)
  }
} else if (platform_of_interest == 'N:METH'){
  for (feature in ftnames_order){
    # 1. Preprocess - to gene symbol
    GENE <- unlist(strsplit(feature, ':'))[4]
    # 2. Hallmark Mapping
    halls <- mappings[mappings$human_gene_symbol==GENE,]$gs_name %>% as.vector()
    # Multiple Hallmarks
    n_hallmarks <- c(n_hallmarks, length(halls))
    pooled_hallmarks <- c(pooled_hallmarks, halls)
  }
}

# 2. Plot fig: How many hallmarks is a feature associated with?
tiff(
  paste(
    args$cancer,
    '_hallmark_fts',
    unlist(strsplit(platform_of_interest, ':'))[2],
    '.tiff',
    sep=''
  ),
  width = 1400,
  height = 1200,
  res = 200,
  compression = "none"
)
df2 <- table(n_hallmarks) %>% as.data.frame()
colnames(df2)<- c('nHallmarks', 'Freq')
ggplot(data=df2, aes(x=nHallmarks, y=Freq, fill='blue')) +
  geom_bar(stat='identity') +
  geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
  ggtitle('Histogram - How many hallmarks is a feature associated with?') +
  theme(legend.position = "none")
dev.off()

# 3. Plot fig: What hallmarks are fts most often associated with?
tiff(
  paste(
    args$cancer,
    '_hallmark_hist',
    unlist(strsplit(platform_of_interest, ':'))[2],
    '.tiff',
    sep=''
  ),
  width = 1400,
  height = 1200,
  res = 200,
  compression = "none"
)
df2 <- sort(table(pooled_hallmarks), decreasing=T) %>% as.data.frame()
colnames(df2)<- c('Hallmarks', 'Freq')
ggplot(data=df2, aes(x=Hallmarks, y=Freq, fill='blue')) +
  geom_bar(stat='identity') +
  geom_text(aes(label=Freq), hjust=-0.5, size=1.5) +
  ggtitle('What hallmarks are fts most often associated with') +
  theme(text = element_text(size=6)) +
  theme(legend.position = "none") +
  coord_flip()
dev.off()



######
# Section 3: Heatmap with Hallmarks and Feature Importance (Top 5)
######
setwd(args$outdir)
# Find top hallmarks from Pathway NES score
importance <- Hallmark.nes.space[,args$cancer]
importance <- sort(importance, decreasing = TRUE)
top_NES <- names(importance[1:5])

header_jadbio <- models['JADBIO']
header_cforest <- models['CForest']
header_aklimate <- models['AKLIMATE']
header_subscope <- models['SubSCOPE']
header_skgrid <- models['SKGrid']

subtype_ha <- subtype_annotation

# Hallmarks by Pathway NES score
vals_1_NES <- build_hallmark_vect(top_NES[1],ftnames_order, platform_of_interest)
vals_2_NES <- build_hallmark_vect(top_NES[2],ftnames_order, platform_of_interest)
vals_3_NES <- build_hallmark_vect(top_NES[3],ftnames_order, platform_of_interest)
vals_4_NES <- build_hallmark_vect(top_NES[4],ftnames_order, platform_of_interest)
vals_5_NES <- build_hallmark_vect(top_NES[5],ftnames_order, platform_of_interest)

# Build annotation bars of teams feature sets.
# 1A. df of all teams. match ft order in heatmap
team_df<- df_fts %>% filter(featureID %in% ftnames_order) %>% arrange(match(featureID, ftnames_order))

# 2. Create MinMax Values where appropriate
jadbio <- team_df %>% pull(header_jadbio) %>% as.character()
cforest <- team_df %>% pull(header_cforest) %>% as.character()
aklimate_minmax <- normalize_data(imp_aklimate, team_df)
subscope <- team_df %>% pull(header_subscope) %>% as.character()
skgrid <- team_df %>% pull(header_skgrid) %>% as.character()

# Build annotation
col_annot <- HeatmapAnnotation(
  # Names of Annot Bars
  annotation_label  = gt_render(
    c(
      'AKLIMATE', "SubSCOPE", "CloudForest", "JADBio", "SciKitGrid",
      "nTeams",
      paste(top_NES[1], ' (n=',gene_set_size(top_NES[1]), ')', sep = ''),
      paste(top_NES[2], ' (n=',gene_set_size(top_NES[2]), ')', sep = ''),
      paste(top_NES[3], ' (n=',gene_set_size(top_NES[3]), ')', sep = ''),
      paste(top_NES[4], ' (n=',gene_set_size(top_NES[4]), ')', sep = ''),
      paste(top_NES[5], ' (n=',gene_set_size(top_NES[5]), ')', sep = '')
    )
  ),

  # A. ft binary membership
  AKLIMATE = aklimate_minmax,
  SubSCOPE = subscope,
  CloudForest = cforest,
  JADBio = jadbio,
  SciKitGrid = skgrid,

  # B. N teams selected
  nTeams= anno_barplot (
    team_df$Total,
    bar_width=1,
    axis_param = list(side = "right", facing='outside')
  ),

  # C. Version 2: Hallmarks by NES
  hallmark1 = vals_1_NES,
  hallmark2 = vals_2_NES,
  hallmark3 = vals_3_NES,
  hallmark4 = vals_4_NES,
  hallmark5 = vals_5_NES,

  col = list(
    AKLIMATE =  colorRamp2(c(0, 1), c("#333333", "cadetblue1")),
    SubSCOPE =  c('0' = "#333333", '1' = "palegreen2"),
    CloudForest =  c('0' = "#333333", '1' = "mediumpurple1"),
    JADBio = c('0' = "#333333", '1' = "#D55B5B"),
    SciKitGrid =  c('0' = "#333333", '1' = "#EFA9A9"),
    hallmark1 = c('0' = "#333333", '1' = "darkgoldenrod3"),
    hallmark2 = c('0' = "#333333", '1' = "darkgoldenrod3"),
    hallmark3 = c('0' = "#333333", '1' = "darkgoldenrod3"),
    hallmark4 = c('0' = "#333333", '1' = "darkgoldenrod3"),
    hallmark5 = c('0' = "#333333", '1' = "darkgoldenrod3")
  ),
  show_legend = FALSE,
  gap = unit(c(0,0,0,0,1,1,0,0,0,0), 'mm')
)

# Plot
ht_rows <- nrow(mat2)
ht_cols <- ncol(mat2)
plat <- unlist(strsplit(platform_of_interest, ':'))[2]
#if z scores > add to heatmap legend
if ( plat == 'METH' || plat == 'GEXP' || plat == 'MIR' ){
  fig <- Heatmap(
    mat2, #each col will have mean 0, sd 1
    # width = unit(10, 'cm'),
    # height = unit(10, 'cm'),
    name = paste('Z-scores', plat, sep='\n'),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = paste('Selected Features (n=', ht_cols, ')', sep=''),
    column_title_gp = gpar(fontsize = 12),
    row_title = paste('Samples (n=', ht_rows, ')', sep=''),
    row_title_gp = gpar(fontsize = 12),
    right_annotation = subtype_ha,
    bottom_annotation = col_annot,
    row_title_side = "right",
    use_raster = TRUE,
    na_col = 'white',
    col = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red'))
  )
} else {
  fig <- Heatmap(
    mat2, #each col will have mean 0, sd 1
    # width = unit(10, 'cm'),
    # height = unit(10, 'cm'),
    name = str_to_title(plat),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = paste('Selected Features (n=', ht_cols, ')', sep=''),
    column_title_gp = gpar(fontsize = 12),
    row_title = paste('Samples (n=', ht_rows, ')', sep=''),
    row_title_gp = gpar(fontsize = 12),
    right_annotation = subtype_ha,
    bottom_annotation = col_annot,
    row_title_side = "right",
    use_raster = TRUE,
    na_col = 'white',
    col = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red'))
  )
}

# Set up saving fig packet
tiff(
  paste(
    args$cancer,
    '_ht_ftimportance_top5',
    unlist(strsplit(platform_of_interest, ':'))[2],
    '.tiff',
    sep=''
  ),
  width = 1600,
  height = 1400,
  res = 200,
  compression = "none"
)
draw(fig,heatmap_legend_side = c('right'))
dev.off()
