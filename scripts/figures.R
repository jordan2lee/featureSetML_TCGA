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
source('func/draw_upset.R')
source('func/save_fig.R')
source('func/draw_main_heatmap.R')
parser <- ArgumentParser()
parser$add_argument('-min', '--min_n_team_overlap', type='character', help='min number of teams to show overlaps on all heatmaps')
parser$add_argument("-c", "--cancer", type='character', help='cancer cohort, all capitalized')
parser$add_argument('-m', '--max_ftsize', type='integer', help='upset plot max value for team ft set size plot')
parser$add_argument('-upset', '--outdir_upset', type='character', help='output dir for upset plots')
parser$add_argument('-os', '--outdir_ht', type='character', help='supplemental dir for heatmap plots')
parser$add_argument('-t', '--input_team_display', type='character', help='desired team order in plots, but in reverse order')
parser$add_argument('-s', '--show_features', action='store_true', default=FALSE, help='boolean if show features on heatmap, include flag if want to see figures')
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
yes_scale <- c('N:GEXP','N:MIR') # which fts to scale

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

# D. PAM50 for breast cancer
f <- '/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/src/brca_pam50_hits.tsv'
pam <- fread(f) %>% as.data.frame(row.names=1)
pam <- pam %>% select(-V1) %>% colnames()

###### PART 1: UPSET PLOT ######
upset_fig <- get_upset(args$cancer, args$input_team_display, args$max_ftsize, get_ymax_upset(args$cancer))
setwd(args$outdir_upset)
image_name <- paste('upsetplot_', args$cancer, '.tiff', sep='')
image_capture(image_name)
upset_fig
dev.off()
print('completed upset plot - mode distinct')









###### PART 2: HEATMAP #########
#####
# Set up
#####
# # Define platform for hallmark heatmap
# platform_of_interest <- get_platform_of_interest(args$cancer)

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
setwd(paste(args$outdir_ht, '/supplemental/', sep=''))
# Get models
models <- model2team(df_fts)
# print('##################')
# print(models)
# print('##################')
# Set up saving fig packet
for (prefix in platforms){

  # Create only base heatmap without hallmarks for miRNA - due to symbols not in hallmark db
  if (prefix == 'N:MIR'){
    print(paste('WORKING ON:', prefix, sep=' '))
    # A. Base Heatmap
    results_list <- get_base_heatmap(
      prefix,
      args$cancer,
      models['JADBIO'],
      models['CloudForest'],
      models['AKLIMATE'],
      models['SubSCOPE'],
      models['ScikitGrid'],
      yes_scale
    )
    # Save figure and exit
    figure <- results_list[['figure']]
    print('in mir loop')
    # If not null fig then save
    if (is.null(figure) == FALSE){
      image_name <- paste(args$cancer,'_heatmap_basic_',unlist(strsplit(prefix, ':'))[2],'.tiff',sep='')
      image_capture(image_name)
      draw(figure, merge_legend = TRUE,legend_grouping ='original', heatmap_legend_side = c('bottom'))
      # print(figure)
      dev.off()
    }
    print('in mir loop DONE')

  } else {
    print(paste("WORKING ON", prefix, sep=' '))
    # A. Base Heatmap
    results_list <- get_base_heatmap(
      prefix,
      args$cancer,
      models['JADBIO'],
      models['CloudForest'],
      models['AKLIMATE'],
      models['SubSCOPE'],
      models['ScikitGrid'],
      yes_scale
    )
    mat2 <- results_list[['results_matrix']]
    ftnames_order <- results_list[['results_ft_order']]
    subtype_annotation <- results_list[['subtype_annotation']]
    n_fts <- length(ftnames_order)
    print(paste('n fts are', n_fts, sep =' '))

    # Exit if have no features to plot
    if (n_fts != 0){


      ######
      # Section 2: Explore the hallmarks associated with GEXP features
      ######
      # 1. Prep: look up if in hallmark mappings table
      n_hallmarks <- c()
      pooled_hallmarks <- c()
      symbols <- c()
      if (prefix == 'N:GEXP'){
        for (feature in ftnames_order){
          # 1. Preprocess - to gene symbol
          GENE <- unlist(strsplit(feature, '::'))[2]
          GENE <- unlist(strsplit(GENE, ':'))[1]
          # 2. Hallmark Mapping
          halls <- mappings[mappings$human_gene_symbol==GENE,]$gs_name %>% as.vector()
          # Multiple Hallmarks
          n_hallmarks <- c(n_hallmarks, length(halls))
          pooled_hallmarks <- c(pooled_hallmarks, halls)
          # TODO dev for symbol reporting on ht
          symbols <- c(symbols, GENE)
        }
      } else if (prefix == 'N:METH' || prefix == 'B:MUTA' || prefix == 'I:CNVR'){
        for (feature in ftnames_order){
          # 1. Preprocess - to gene symbol
          GENE <- unlist(strsplit(feature, ':'))[4]
          # 2. Hallmark Mapping
          halls <- mappings[mappings$human_gene_symbol==GENE,]$gs_name %>% as.vector()
          # Multiple Hallmarks
          n_hallmarks <- c(n_hallmarks, length(halls))
          pooled_hallmarks <- c(pooled_hallmarks, halls)
        }
      } else if (prefix == 'N:MIR'){
        print('MIR fts will not be mapped to hallmarks')
      }
      # 2. Plot fig: How many hallmarks is a feature associated with?
      image_name <- paste(args$cancer, '_bars_ft2hall_', unlist(strsplit(prefix, ':'))[2],'.tiff',sep='')
      image_capture(image_name)
      df2 <- table(n_hallmarks) %>% as.data.frame()
      colnames(df2)<- c('nHallmarks', 'Freq')
      p <- ggplot(data=df2, aes(x=nHallmarks, y=Freq, fill='blue')) +
        geom_bar(stat='identity') +
        geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
        ggtitle('Histogram - How many hallmarks is a feature associated with?') +
        theme(legend.position = "none")
      print(p)
      dev.off()

      # 3. Plot fig: What hallmarks are fts most often associated with?
      image_name <- paste(args$cancer, '_bars_hall_', unlist(strsplit(prefix, ':'))[2], '.tiff', sep='')
      image_capture(image_name)
      df2 <- sort(table(pooled_hallmarks), decreasing=T) %>% as.data.frame()
      # If there are no hallmarks
      if (nrow(df2) == 0 ){
        print('no hallmarks found to features, skipping creation of hallmark_hist fig')
        p <- ggplot(data.frame(0)) +
          ggtitle('No Hallmarks mapped to features (or no features)')
        print(p)
        dev.off()
      } else {
        # If there is only one hallmark force df to have 2 cols instead of 1
        if (nrow(df2) == 1){
          h <- rownames(df2)
          df2<- cbind(h, df2)
          rownames(df2) <- 1
        }
        # Create fig
        colnames(df2)<- c('Hallmarks', 'Freq')
        p <- ggplot(data=df2, aes(x=Hallmarks, y=Freq, fill='blue')) +
          geom_bar(stat='identity') +
          geom_text(aes(label=Freq), hjust=-0.5, size=1.5) +
          ggtitle('What hallmarks are fts most often associated with') +
          theme(text = element_text(size=6)) +
          theme(legend.position = "none") +
          coord_flip()
        print(p)
        dev.off()
      }

      ######
      # Set up for PAM if BRCA
      #####
      # otherwise in_pam == NULL
      # Mark if ft is in PAM50
      if (args$cancer == 'BRCA'){
        in_pam <- c()
        for (i in seq(1, length(ftnames_order))){
          f <- ftnames_order[i]
          if (f %in% pam == TRUE){
            in_pam <- c(in_pam, "+")
          } else {
            in_pam <- c(in_pam, '')
          }
        }
      } else {
        in_pam = NULL
      }


      ######
      # Section 3: Heatmap with Hallmarks and Feature Importance (Top 5)
      ######
      # Find top hallmarks from Pathway NES score
      importance <- Hallmark.nes.space[,args$cancer]
      importance <- sort(importance, decreasing = TRUE)
      top_NES <- names(importance[1:5])

      header_aklimate <- models['AKLIMATE']
      header_subscope <- models['SubSCOPE']
      header_cforest <- models['CloudForest']
      header_jadbio <- models['JADBIO']
      header_skgrid <- models['ScikitGrid']

      subtype_ha <- subtype_annotation

      # # Hallmarks by Pathway NES score
      # vals_1_NES <- build_hallmark_vect(top_NES[1],ftnames_order, prefix)
      # vals_2_NES <- build_hallmark_vect(top_NES[2],ftnames_order, prefix)
      # vals_3_NES <- build_hallmark_vect(top_NES[3],ftnames_order, prefix)
      # vals_4_NES <- build_hallmark_vect(top_NES[4],ftnames_order, prefix)
      # vals_5_NES <- build_hallmark_vect(top_NES[5],ftnames_order, prefix)

      # Build annotation bars of teams feature sets.
      # 1A. df of all teams. match ft order in heatmap
      team_df<- df_fts %>% filter(featureID %in% ftnames_order) %>% arrange(match(featureID, ftnames_order))

      # 2. Create MinMax Values where appropriate
      aklimate_minmax <- normalize_data(imp_aklimate, team_df)
      subscope <- team_df %>% pull(header_subscope) %>% as.character()
      cforest <- team_df %>% pull(header_cforest) %>% as.character()
      jadbio <- team_df %>% pull(header_jadbio) %>% as.character()
      skgrid <- team_df %>% pull(header_skgrid) %>% as.character()
      # Build annotation
      col_annot <- HeatmapAnnotation(
        # Names of Annot Bars
        annotation_label  = gt_render(
          c(
            'Model Overlap', 'AKLIMATE', "SubSCOPE", "Cloud Forest", "JADBio", "SciKitGrid",
            paste(top_NES[1], ' (n=',gene_set_size(top_NES[1]), ')', sep = ''),
            paste(top_NES[2], ' (n=',gene_set_size(top_NES[2]), ')', sep = ''),
            paste(top_NES[3], ' (n=',gene_set_size(top_NES[3]), ')', sep = ''),
            paste(top_NES[4], ' (n=',gene_set_size(top_NES[4]), ')', sep = ''),
            paste(top_NES[5], ' (n=',gene_set_size(top_NES[5]), ')', sep = '')
          )
        ),
        # A. N teams selected
        nTeams= anno_barplot(
          team_df$Total,
          bar_width=1,
          gp = gpar(fill = 'darkgray', col = 'azure4'),
          border = FALSE,
          rot = 45,
          axis_param = list(side = "right", facing='outside', gp=gpar(fontsize=5)) #yaxis size
        ),

        # B. ft binary membership
        "AKLIMATE\nmin-max" = aklimate_minmax,
        "SubSCOPE" = subscope,
        "Cloud Forest" = cforest,
        "JADBio" = jadbio,
        "SciKitGrid" = skgrid,

        annotation_name_rot = 0,

        # # C. Version 2: Hallmarks by NES
        # hallmark1 = vals_1_NES,
        # hallmark2 = vals_2_NES,
        # hallmark3 = vals_3_NES,
        # hallmark4 = vals_4_NES,
        # hallmark5 = vals_5_NES,

        col = list(
          'AKLIMATE\nmin-max' =  colorRamp2(c(0, 0.05, 1), c("#333333", "cadetblue4", "#BFFEFF")),
          "SubSCOPE" =  c('0' = "#333333", '1' = "#AEFEB0"),
          "Cloud Forest" =  c('0' = "#333333", '1' = "#BFBFFF"),
          "JADBio" = c('0' = "#333333", '1' = "#FBBD91"),
          "SciKitGrid" =  c('0' = "#333333", '1' = "#FCC0BF"),
          hallmark1 = c('0' = "#333333", '1' = "azure4"),
          hallmark2 = c('0' = "#333333", '1' = "azure4"),
          hallmark3 = c('0' = "#333333", '1' = "azure4"),
          hallmark4 = c('0' = "#333333", '1' = "azure4"),
          hallmark5 = c('0' = "#333333", '1' = "azure4")
        ),
        show_legend = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
        gp = gpar(fontsize = 1), # grid all col annot
        annotation_name_gp= gpar(fontsize = 8),
        gap = unit(c(1,0,0,0,0), 'mm')
      )

      # prepare for plotting
      ht_rows <- nrow(mat2)
      ht_cols <- ncol(mat2)
      plat <- unlist(strsplit(prefix, ':'))[2]
      col_title = paste(title_info(plat), ' Features Selected by ≥ 2 Teams (n=', ht_cols, ')', sep='')

      # Create fig object
      if ( prefix %in% yes_scale ){
        main_ht_name = paste(platform_display_text(plat), 'z-score', sep=' ')
        if (prefix == 'N:MIR' ){
          fig <- get_main_heatmap(plat, main_ht_name, args$cancer)
        } else if (prefix == 'N:GEXP') {
          fig <- get_main_heatmap(plat, main_ht_name, args$cancer)
        }
      } else if (plat == 'MUTA' || plat == 'METH' || plat == 'CNVR'){
        fig <- get_main_heatmap(plat, platform_display_text(plat), args$cancer)
      }

      # Set up saving fig packet
      if (args$show_features == TRUE){
        image_name <- paste(args$cancer, '_heatmap_', unlist(strsplit(prefix, ':'))[2], '_NAMES', '.tiff', sep='')
        image_capture(image_name)
        draw(fig,merge_legend = TRUE,legend_grouping ='original', heatmap_legend_side = c('bottom'))
        dev.off()
      } else {
        image_name <- paste(args$cancer, '_heatmap_', unlist(strsplit(prefix, ':'))[2], '.tiff', sep='')
        image_capture(image_name)
        draw(fig,merge_legend = TRUE,legend_grouping ='original', heatmap_legend_side = c('bottom'))
        dev.off()
    }
      # Save tsv of heatmap data
      write.table(mat2, file =paste(args$cancer, '_matrix_heatmap_', unlist(strsplit(prefix, ':'))[2], '.tsv', sep = ''), sep='\t', row.names=TRUE, col.names = TRUE)
    } else {

      print(paste('No features shared by at least', args$min_n_team_overlap, 'groups for', prefix, sep = ' '))
    }
  }
}
