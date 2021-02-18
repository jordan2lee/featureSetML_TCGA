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

# create parser object
parser <- ArgumentParser()
parser$add_argument("-c", "--cancer", type='character', help='cancer cohort, all capitalized')
parser$add_argument('-op', '--outdir', type='character', help='output dir')
parser$add_argument('-os', '--supplemental', type='character', help='supplemental dir')
args <- parser$parse_args()

######
# Hardcoded
######
# cancer <- 'GEA' # 'GEA', 'LGGGBM', 'BRCA'
# outdir ='data/figure_panel_b/'
file_imp_aklimate <- paste(
  'src/aklimate_feature_importance_scores_20200807/AKLIMATE_TEMPLATE_FULLY_TRAINED_',
  args$cancer,
  '_20200423_aklimate_ranked_feature_importance.tsv',
  sep=''
)
yes_scale <- c('N:METH', 'N:GEXP','N:MIR') # which fts to scale




model2team <- function(df){
  #' Read team model names and map to team name.
#' Returns named vector of specific team models.
#' Names: "JADBIO", "CForest", "AKLIMATE", "SubSCOPE", "SKGrid"
headers_vector <- c()
names_vector <- c()
models <- colnames(df)[!colnames(df) %in% c('featureID', 'Total')]

for (header in models){
  if (grepl('gnosis', header, fixed=TRUE)){
    headers_vector <- c(headers_vector, header)
    names_vector <- c(names_vector, 'JADBIO')
  }
  else if (grepl('CF', header, fixed=TRUE)){
    headers_vector <- c(headers_vector, header)
    names_vector <- c(names_vector, 'CForest')
  }
  else if (grepl('nn_jg', header, fixed=TRUE)){
    headers_vector <- c(headers_vector, header)
    names_vector <- c(names_vector, 'SubSCOPE')
  }
  else if (grepl('AKLIMATE', header, fixed=TRUE)){
    headers_vector <- c(headers_vector, header)
    names_vector <- c(names_vector, 'AKLIMATE')
  }
  else {
    headers_vector <- c(headers_vector, header)
    names_vector <- c(names_vector, 'SKGrid')
  }
}
names(headers_vector) <- names_vector
return(headers_vector)
}


get_colors <- function(df){
  #' Pull color codes based on total number of subtypes found in "Labels" column
  #' Input df that contains column "Labels" where subtypes are
  #' Returns vector of set heatmap color codes for subtypes
  nsubs <- df %>%
    select(Labels) %>%
    unique() %>%
    nrow() %>%
    as.integer()
  if (nsubs == 4){
    color_codes <- c(
      "1" = 'salmon4',
      "2"='red3',
      "3"='orangered',
      "4"='orange1'
    )
  } else if (nsubs == 7){
    color_codes <- c(
      "1" = '#F58748',
      "2" ='#A95757',
      "3"='#4A2918',
      "4"='red3',
      "5"='orangered',
      "6"= 'salmon',
      "7" = 'orange1'
    )
  }
  return(color_codes)
}



get_base_heatmap <- function(prefix, cancer, header_jadbio, header_cforest, header_aklimate, header_subscope, header_skgrid){
  #' Create base heatmap - no hallmark info. For exploratory purposes
  ######
  # Preprocess
  #####
  # A. Order by subtype
  df_transform <- df %>% arrange(Labels)
  # B. Column annotation
  s_matrix <- pull(df_transform, Labels) %>% as.vector()
  s_matrix <- sapply(strsplit(s_matrix, '_'), "[", 2) %>% as.matrix()
  subtype_ha <- rowAnnotation(
    Subtype = s_matrix,
    na_col = 'grey',
    col = list(
      Subtype = get_colors(df)
    ),
    show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm") # width
  )
  # C. Select data type
  df_transform <- df_transform %>%
    select(-Labels) %>%
    select(-all_of(cancer)) %>%
    select(starts_with(prefix))
  if (ncol(df_transform) != 0){

    mat <- df_transform %>%
      as.matrix() %>%
      t()
    print(prefix)

    #####
    # 1. Generate temp Heatmap
    # and pull row/col order
    #####
    # Scale if appropriate
    #z-scores == each ft row will have mean 0, sd 1. omit NAs
    if (prefix %in% yes_scale){
      mat <- scale(t(mat), center=TRUE, scale=TRUE)
    } else {
      mat <- t(mat) #flip for heatmap looks
    }
    # Heatmap
    ht_rows <- nrow(mat)
    ht_cols <- ncol(mat)
    fig <- Heatmap(
      mat,
      name = 'first heatmap',
      cluster_rows = FALSE,
      clustering_distance_columns = "euclidean",
      clustering_method_columns = "ward.D",
      cluster_columns = TRUE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_title = paste('Selected Features (n=', ht_rows, ')', sep=''),
      row_title = paste('Samples (n=', ht_cols, ')', sep=''),
      right_annotation = subtype_ha
    )
    ####
    # Add team annotation bar
    ####
    # Ordering
    # 1. Get order of features post heatmap clustering
    heatmap_order <- column_order(fig) # index vector
    ftnames_order <- c() # featurename vector
    for (i in heatmap_order){
      add_ft <- colnames(mat)[i]
      ftnames_order <- c(ftnames_order, add_ft)
    }
    # 2. Get new matrix that is ordered by heatmap clustering
    mat2 <- mat[,match(ftnames_order, colnames(mat))]
    #####
    # Build annotation bars of teams feature sets
    #####
    # 1. df of all teams. match ft order in heatmap
    team_df<- df_fts %>%
      filter(featureID %in% ftnames_order) %>%
      arrange(match(featureID, ftnames_order))
    # 2. Pull just the team of interest
    jadbio <- team_df %>%
      pull(header_jadbio) %>%
      as.character()
    cforest <- team_df %>%
      pull(header_cforest) %>%
      as.character()
    aklimate <- team_df %>%
      pull(header_aklimate) %>%
      as.character()
    subscope <- team_df %>%
      pull(header_subscope) %>%
      as.character()
    skgrid <- team_df %>%
      pull(header_skgrid) %>%
      as.character()
    team_list <- HeatmapAnnotation(
      JADBio = jadbio,
      CloudForest = cforest,
      AKLIMATE = aklimate,
      SubSCOPE = subscope,
      SciKitGrid = skgrid,
      col = list(
        JADBio = c('0' = "#333333", '1' = "#D55B5B"),
        CloudForest = c('0' = "#333333", '1' = "mediumpurple1"),
        AKLIMATE = c('0' = "#333333", '1' = "cadetblue1"),
        SubSCOPE = c('0' = "#333333", '1' = "palegreen2"),
        SciKitGrid = c('0' = "#333333", '1' = "#EFA9A9")
      ),
      show_legend = FALSE,
      nTeams= anno_barplot (
        team_df$Total,
        bar_width=1,
        axis_param = list(side = "right", facing='outside')
      )
      # simple_anno_size = unit(2, 'mm') # height
    )
    #####
    # 3. Output Heatmap
    #####
    # Set up saving fig packet
    tiff(
      paste(
        cancer,
        '_ht_base_',
        unlist(strsplit(prefix, ':'))[2],
        '.tiff',
        sep=''
      ),
      width = 1400,
      height = 1200,
      res = 200,
      compression = "none"
    )
    # Draw
    # handle if only 1 feature in heatmap
    if (length(ftnames_order)==1){
      # Sanity check
      ht_rows <- nrow(mat2)
      ht_cols <- ncol(mat2)
      assert('Assertion Error: in if loop for 1 ft but mat2 object dim() does not match',
             is.null(ht_rows) && is.null(ht_cols)
      )

      # If only one feature set ht_rows/cols manually
      ht_rows <- length(mat2)
      ht_cols <- 1

      fig <- Heatmap(
        mat2,
        name = str_to_title(prefix),
        # width = unit(12, 'cm'),
        # height = unit(12, 'cm'),
        cluster_rows = FALSE,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D",
        # column_order = ftnames_order, # NO ORDERING NEEDED
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_title = paste('Features (n=', ht_cols, ')', sep=''),
        row_title = paste('Samples (n=', ht_rows, ')', sep=''),
        right_annotation = subtype_ha,
        bottom_annotation = team_list,
        na_col = 'white'
      )
      draw(fig)
      dev.off()
    } else {
      ht_rows <- nrow(mat2)
      ht_cols <- ncol(mat2)
      fig <- Heatmap(
        mat2,
        name = str_to_title(prefix),
        # width = unit(12, 'cm'),
        # height = unit(12, 'cm'),
        cluster_rows = FALSE,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D",
        column_order = ftnames_order,
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_title = paste('Features (n=', ht_cols, ')', sep=''),
        row_title = paste('Samples (n=', ht_rows, ')', sep=''),
        right_annotation = subtype_ha,
        bottom_annotation = team_list,
        na_col = 'white'
      )
      draw(fig)
      dev.off()
    }
    print(
      paste(
        'Distance metric = ',
        fig@column_dend_param$distance,
        '. Method = ',
        fig@column_dend_param$method,
        sep=' '
      )
    )
    # Assign to 'output' variables
    return(
      list(
        'results_matrix' = mat2,
        'results_ft_order' = ftnames_order,
        'subtype_annotation' = subtype_ha
      )
    )
  } else {
    return(
      list(
        'results_matrix' = NULL,
        'results_ft_order' = NULL,
        'subtype_annotation' = NULL
      )
    )
  }
}


build_hallmark_vect <- function(hallmark, ftnames_order, platform_of_interest){
  #' Create hallmark vector
  #' for builidng hallmark heatmap
  fts_checked <- c() # sanity check
  hallmark_present <- c()
  for (feature in ftnames_order){
    # 1. Preprocess - to gene symbol
    if (platform_of_interest == 'N:GEXP'){
      GENE <- unlist(strsplit(feature, '::'))[2]
      GENE <- unlist(strsplit(GENE, ':'))[1]
    } else if (platform_of_interest == 'N:METH'){
      GENE <- unlist(strsplit(feature, ':'))[4]
    }
    # 2. Hallmark Mapping
    halls <- mappings[mappings$human_gene_symbol==GENE,]$gs_name %>%
      as.vector()
    # If hallmark present
    if (hallmark %in% halls){
      i <- match(hallmark, halls)
      hallmark_present <- c(hallmark_present, 1)
      fts_checked <- c(fts_checked, feature)
    }
    else {
      hallmark_present <- c(hallmark_present, 0)
      fts_checked <- c(fts_checked, feature)
    }
  }
  return(hallmark_present)
}


mm_normal <- function(x) {
  #' Calculate min-max normalization of input values
  return ((x - min(x)) / (max(x) - min(x)))
}


normalize_data <- function(df, team_df){
  #' Purpose: min-max normalization and order features to match heatmap
  #' Input df that contains columns c("features","importance")
  #' Input df of all teams
  #' Return vector of normalized values

  # Min max normalize importance scores
  subset_df <- df[,c('features', 'importance')]
  minmax_vals <- mm_normal(subset_df$importance)
  subset_df <- cbind(subset_df, minmax_vals)

  # Map and order ft importances to pooled ft order
  minmax_norm <- c()
  test_ft_order <- c()
  all_features <- subset_df$features
  for (pooled_ft in team_df$featureID){
    if (pooled_ft %in% all_features){
      mm <- subset_df[subset_df$'features'==pooled_ft,]$minmax_vals
      minmax_norm <- c(minmax_norm, mm)
      test_ft_order <- c(test_ft_order, pooled_ft)
    } else {
      minmax_norm <- c(minmax_norm, 0)
      test_ft_order <- c(test_ft_order, pooled_ft)
    }
  }
  return(minmax_norm)
}

gene_set_size <- function(input_hallmark){
  # subset for hallmark rows
  tab<- mappings[mappings['gs_name']==input_hallmark,]
  # count unique gene symbols
  genes <- unique(unlist(tab['human_gene_symbol']))
  return(length(genes))
}


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

####
# Define platform for hallmark heatmap
####
if (args$cancer == 'BRCA'){
  platform_of_interest <- 'N:GEXP'
} else if (args$cancer == 'LGGGBM'){
  platform_of_interest <- 'N:METH'
} else if (args$cancer == 'GEA'){
  platform_of_interest <- 'N:METH'
}

####
# Build list of data types present
####
if (args$cancer == 'BRCA' | args$cancer == 'GEA' ){
  platforms <- c('N:METH', 'N:MIR', 'I:CNVR', 'B:MUTA', 'N:GEXP')
} else if (args$cancer == 'LGGGBM'){
  platforms <- c('N:METH', 'I:CNVR', 'B:MUTA', 'N:GEXP')
}


########
# Filter for ft rows that occur in >= 2 teams
########
#first for df_fts
df_fts <- df_fts[df_fts['Total']>=2,]
rownames(df_fts) <- 1:nrow(df_fts) #reset row index names
cols_to_keep <- df_fts$featureID
cols_to_keep <- c(args$cancer, 'Labels', cols_to_keep)
#second for df
df <- subset(df, select=cols_to_keep)

##########
### Section 1: Heatmaps for each data type
##########
setwd(args$supplemental)

####
# Create Base Heatmap for ALL DATA TYPES
####
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

########
### Section 2: Explore the hallmarks associated with GEXP features
###########
# One annotation bar (top hallmarks only, all others are not colored)
######
# Explore Hallmarks
######
# 1. Prep
#look up if in hallmark mappings table
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
    'hallmark_fts',
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
    'hallmark_hist',
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
### Section 3: Heatmap with Hallmarks and Feature Importance (Top 5)
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
    'ht_ftimportance_top5',
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
