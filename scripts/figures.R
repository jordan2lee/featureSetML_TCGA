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
sourceRecursive <- function(path) {
  #' Input dir and will source all files and subdirs within it
  dirs <- list.dirs(path, recursive = FALSE)
  files <- list.files(path, pattern = "^.*[Rr]$", include.dirs = FALSE, full.names = TRUE)
  for (f in files)
    source(f)
  for (d in dirs)
    sourceRecursive(d)
}
sourceRecursive("./func")
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

file_imp_scikitgrid <- paste(
  'data/top_model_importances/TOP', args$cancer, 'SKGRID.tsv',
  sep='_'
)
file_imp_subscope <- paste(
  'data/top_model_importances/TOP', args$cancer, 'SUBSCOPE.tsv',
  sep='_'
)
file_imp_cf <- paste(
  'data/top_model_importances/TOP', args$cancer, 'CF.tsv',
  sep='_'
)

jadbio_ft_file <- function(cancer){
  #' Input cancer and output the file containing jadbio ft importances
  jadbio_files <- list(
    'LGGGBM' = 'src/jadbio_ft_importances_f1/LGGGBM_MULTIDATATYPE_markers_outcome_association_multisignature.csv',
    'BRCA' = 'src/jadbio_ft_importances_f1/BRCA_GEXP_markers_outcome_association.csv',
    'COADREAD' = 'src/jadbio_ft_importances_f1/COADREAD_MULTIDATATYPE_markers_outcome_association_multisignature.csv',
    'SKCM' = 'src/jadbio_ft_importances_f1/SKCM_MUTA_markers_outcome_association.csv'
  )
  return(jadbio_files[[cancer]])
}
file_imp_jadbio <- jadbio_ft_file(args$cancer)

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
colnames(imp_jadbio)<- c('features'	,'importance')


# D. PAM50 for breast cancer
f <- '/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/src/brca_pam50_hits.tsv'
pam <- fread(f) %>% as.data.frame(row.names=1)
pam <- pam %>% select(-V1) %>% colnames()

# Methlyation lit support
f <- '/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/src/toshi_meth_lit_search/20220425_TMP_DNA_methylation_features_analysis_COAD.tsv'
coad_lit <- fread(f) %>% as.data.frame(row.names=1)
coad_lit <- coad_lit[,1:8] # if read in null cols
names(coad_lit)[8] <- 'InLit'

f <- '/Users/leejor/Ellrott_Lab/02_ML/08_manuscript/featureSetML_TCGA/src/toshi_meth_lit_search/20220425_TMP_DNA_methylation_features_analysis_LGGGBM.tsv'
lgggbm_lit <- fread(f) %>% as.data.frame(row.names=1)
lgggbm_lit <- lgggbm_lit[,1:8] # if read in null cols
names(lgggbm_lit)[8] <- 'InLit'

###### PART 1: UPSET PLOT ######
upset_fig <- draw_upset(args$cancer, args$input_team_display, args$max_ftsize, get_ymax_upset(args$cancer))
setwd(args$outdir_upset)
image_name <- paste('upsetplot_', args$cancer, '.tiff', sep='')
image_capture_upset(image_name)
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
      yes_scale,
      df,
      Labels
    )
    # Save figure and exit
    figure <- results_list[['figure']]
    print('in mir loop')
    # If not null fig then save
    if (is.null(figure) == FALSE){
      image_name <- paste(args$cancer,'_heatmap_basic_',unlist(strsplit(prefix, ':'))[2],'.tiff',sep='')
      image_capture_ht(image_name)
      draw(figure, merge_legend = TRUE,legend_grouping ='original', heatmap_legend_side = c('bottom'))
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
      yes_scale,
      df,
      Labels
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
      image_capture_ht(image_name)
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
      image_capture_ht(image_name)
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

      ####
      # Set up for meth if COADREAD OR LGGGBM
      ####
      # Mark if ft is in lit support
      # Grab fts based on filter: present in at least 1 lit
      print(paste('working on ', prefix))
      if (args$cancer == 'COADREAD' && prefix == 'N:METH'){
        coad_lit_ft <- coad_lit[coad_lit$InLit > 0,]$DNA_methylation_feature
        lit_support <- get_lit_vector(args$cancer, prefix, coad_lit_ft, ftnames_order)
        print(paste('literature support is: '))
        print(lit_support)
      } else if (args$cancer == 'LGGGBM' && prefix == 'N:METH'){
        lgggbm_lit_ft <- lgggbm_lit[lgggbm_lit$InLit > 0,]$DNA_methylation_feature
        lit_support <- get_lit_vector(args$cancer, prefix, lgggbm_lit_ft, ftnames_order)
        print(paste('literature support is: '))
        print(lit_support)
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
      aklimate_minmax <- normalize_data(imp_aklimate, team_df, header_aklimate)
      subscope <- normalize_data(imp_subscope, team_df, header_subscope)
      cforest <- normalize_data(imp_cf, team_df, header_cforest)
      jadbio <- normalize_data(imp_jadbio, team_df, header_jadbio)
      skgrid <- normalize_data(imp_scikitgrid, team_df, header_skgrid)

      # Build annotation
      col_annot <- HeatmapAnnotation(
        # Names of Annot Bars
        annotation_label  = gt_render(
          c(
            'Model Overlap', 'AKLIMATE', "SubSCOPE", "CloudForest", "JADBio", "SK Grid"
          )
        ),
        # A. N teams selected
        nTeams= anno_barplot(
          team_df$Total,
          bar_width=1,
          gp = gpar(
            fill = 'black',
            col = 'black'
          ),
          border = FALSE,
          rot = 45,
          axis_param = list(
            side = "right",
            facing='outside',
            gp=gpar(
            fontsize=get_gpar('model_overlap_size'),
            fontfamily = get_gpar('font_fam')
            )
          ) #yaxis size
        ),

        # B. ft binary membership
        "AKLIMATE" = aklimate_minmax,
        "SubSCOPE" = subscope,
        "CloudForest" = cforest,
        "JADBio" = jadbio,
        "SK Grid" = skgrid,

        annotation_name_rot = 0,

        col = list(
        'AKLIMATE' =  colorRamp2(c(0, 0.05, 1), c("#BFFEFF", "#1CBAB9", "#085250")),
        "SubSCOPE" = colorRamp2(c(0, 0.05, 1), c("#AEFEB0", "#0DBF59", "#054621")),
        "CloudForest" =  colorRamp2(c(0, 0.05, 1), c("#BFBFFF", "#B2A0EC", "#5138A1")),
        "JADBio" = colorRamp2(c(0, 0.05, 1), c("#FBBD91", "#F8BE99", "#E57228")),
        "SK Grid" =  colorRamp2(c(0, 0.05, 1), c("#FCC0BF", "#ED94B4", "#99355A"))
        ),
        na_col = "white", # color of NA in bottom annot
        show_legend = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
        gp = gpar(fontsize = 1, col = "black"), # grid all col annot
        annotation_name_gp= gpar(fontsize = get_gpar('annot_size'), fontfamily = get_gpar('font_fam')),
        annotation_legend_param = list(
          direction= 'horizontal',
          title_position = "lefttop",  # legend title location
          legend_width = unit(4, "cm"),
          title_gp = gpar(fontsize = get_gpar('legend_size_title'), fontfamily = get_gpar('font_fam')),
          labels_gp = gpar(fontsize = get_gpar('legend_size'), fontfamily = get_gpar('font_fam'))),
        gap = unit(c(2,0,0,0,0), 'mm')
      )

      # prepare for plotting
      ht_rows <- nrow(mat2)
      ht_cols <- ncol(mat2)
      plat <- unlist(strsplit(prefix, ':'))[2]
      col_title = paste(title_info(plat), ' Core Feature Set Selected by ≥ 2 Methods (n=', ht_cols, ')', sep='')
      # col_title = paste(update_cohort_name(args$cancer), ' ', title_info(plat), ' Features Selected by ≥ 2 Methods (n=', ht_cols, ')', sep='')

      # Create fig object
      if ( prefix %in% yes_scale ){
        main_ht_name = paste(platform_display_text(plat), 'z-score', sep=' ')
        if (prefix == 'N:MIR' ){
          fig <- get_main_heatmap(plat, main_ht_name, args$cancer)
        } else if (prefix == 'N:GEXP') {
          fig <- get_main_heatmap(plat, main_ht_name, args$cancer)
        } else {
          print('ERROR. yes_scale obj contains unrecognized content')
        }
      } else if (plat == 'MUTA' || plat == 'METH' || plat == 'CNVR'){
        main_ht_name = platform_display_text(plat)
        fig <- get_main_heatmap(plat, platform_display_text(plat), args$cancer)
      }

      ########### TESTING
      # Purpose: see if can save different image size

      ########### TESTING

      # Set up saving fig packet
      if (args$show_features == TRUE){
        image_name <- paste(args$cancer, '_heatmap_', unlist(strsplit(prefix, ':'))[2], '_NAMES', '.tiff', sep='')
        image_capture_ht(image_name)
        draw(fig,merge_legend = TRUE, heatmap_legend_side = c('bottom')) # opt add: legend_grouping ='original'
        dev.off()
      } else {
        image_name <- paste(args$cancer, '_heatmap_', unlist(strsplit(prefix, ':'))[2], '.tiff', sep='')
        image_capture_ht(image_name)
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
