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
source('func/config_analysis.R')
source('func/prep_for_heatmap.R')
source('func/draw_base_heatmap.R')
source('func/draw_upset.R')
source('func/save_fig.R')
source('func/draw_main_heatmap.R')
source('func/build_ht_func.R')
source('func/ht_col_annot.R')
source('func/ht_legend.R')
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
######################
# # TODO update from BRCA hardcoded here
file_imp_scikitgrid <- paste(
  'data/ft_importances/TOP_skgrid_BRCA_fbedeBIC_perplatformALL_BRCA.tsv',
  sep='\t'
)
file_imp_subscope <- paste(
  'data/ft_importances/TOP_subSCOPE-GEXP_2021-04-21_bootstrapfeatures_BRCA_BRCA.tsv',
  sep='\t'
)
file_imp_cf <- paste(
  'data/ft_importances/TOP_CF_BRCA_All_Top_50_BRCA.tsv',
  sep='\t'
)
file_imp_jadbio <- paste(
  'data/ft_importances/TOP_jadbio_BRCA_GEXP_cumulative_feature_set25_BRCA.tsv',
  sep='\t'
)
######################
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
imp_scikitgrid <- fread(
  file_imp_scikitgrid
) %>% as.data.frame()
imp_subscope <- fread(
  file_imp_subscope
) %>% as.data.frame()
imp_cf <- fread(
  file_imp_cf
) %>% as.data.frame()
imp_jadbio <- fread(
  file_imp_jadbio
) %>% as.data.frame()

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

  # Create only base heatmap for miRNA
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
      # Section 3: Heatmap with Feature Importance (Top 5)
      ######
      header_aklimate <- models['AKLIMATE']
      header_subscope <- models['SubSCOPE']
      header_cforest <- models['CloudForest']
      header_jadbio <- models['JADBIO']
      header_skgrid <- models['ScikitGrid']

      subtype_ha <- subtype_annotation

      # Build annotation bars of teams feature sets.
      # 1A. df of all teams. match ft order in heatmap
      team_df<- df_fts %>% filter(featureID %in% ftnames_order) %>% arrange(match(featureID, ftnames_order))

      # 2. Create MinMax Values where appropriate
      aklimate_minmax <- normalize_data(imp_aklimate, team_df)
      subscope <- normalize_data(imp_subscope, team_df)
      cforest <- normalize_data(imp_cf, team_df)
      jadbio <- normalize_data(imp_jadbio, team_df)
      write.table(imp_jadbio, file ='/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/data/TEST_imp_jadbio.tsv', sep='\t', row.names=FALSE, col.names = TRUE)
      write.table(jadbio, file ='/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/data/TEST_jadbio.tsv', sep='\t', row.names=TRUE, col.names = TRUE)
      skgrid <- normalize_data(imp_scikitgrid, team_df)

      # Build annotation
      col_annot <- HeatmapAnnotation(
        # Names of Annot Bars
        annotation_label  = gt_render(
          c(
            'Model Overlap', 'AKLIMATE', "SubSCOPE", "Cloud Forest", "JADBio", "SciKitGrid"
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

        col = list(
          'AKLIMATE\nmin-max' =  colorRamp2(c(0, 0.05, 1), c("#333333", "cadetblue4", "#BFFEFF")),
          "SubSCOPE" =  c('0' = "#333333", '1' = "#AEFEB0"),
          "Cloud Forest" =  c('0' = "#333333", '1' = "#BFBFFF"),
          "JADBio" = c('0' = "#333333", '1' = "#FBBD91"),
          "SciKitGrid" =  c('0' = "#333333", '1' = "#FCC0BF")
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
        draw(fig,merge_legend = TRUE, heatmap_legend_side = c('bottom')) # opt add: legend_grouping ='original'
        dev.off()
      } else {
        image_name <- paste(args$cancer, '_heatmap_', unlist(strsplit(prefix, ':'))[2], '.tiff', sep='')
        image_capture(image_name)
        draw(fig,merge_legend = TRUE,heatmap_legend_side = c('bottom')) # opt add: legend_grouping ='original'
        dev.off()
    }
      # Save tsv of heatmap data
      write.table(mat2, file =paste(args$cancer, '_matrix_heatmap_', unlist(strsplit(prefix, ':'))[2], '.tsv', sep = ''), sep='\t', row.names=TRUE, col.names = TRUE)
    } else {

      print(paste('No features shared by at least', args$min_n_team_overlap, 'groups for', prefix, sep = ' '))
    }
  }
}
