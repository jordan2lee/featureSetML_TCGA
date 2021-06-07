get_top_annot <- function(k){
  #' Gene indcies to show on heatmap
  #' If input cancer not in list then will return NULL
  #' and heatmap will plot without gene names shown
  l1 <- list(
    'BRCA_GEXP' = HeatmapAnnotation(
      highlight_fts = anno_mark(
        at = c(4,6,12,17,28,31,37,38,43,49,54),
        labels = symbols[c(4,6,12,17,28,31,37,38,43,49,54)],
        labels_gp = gpar(
          fontsize = get_gpar('symbol_size'),
          fontfamily = get_gpar('font_fam'),
          which = 'column',
          link_gp = gpar(
            fontsize = get_gpar('symbol_size'),
            fontfamily = get_gpar('font_fam')
          )
        )
      ),
      annotation_width=unit(1, 'mm')
    )
  )
  return(l1[[k]])
}

dev_bottom_annot <- function(k, in_pam){
  #' DEV
  l1 <- list(
    'BRCA_GEXP' = HeatmapAnnotation(
      annotation_label  = gt_render(
        c(
          'Model Overlap', 'AKLIMATE', "SubSCOPE", "Cloud Forest", "JADBio", "SciKitGrid"
        )
      ),
      # A. N teams selected
      nTeams= anno_barplot(
        team_df$Total,
        bar_width=1,
        gp = gpar(
          fill = 'darkgray',
          col = 'azure4'
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
      "Cloud Forest" = cforest,
      "JADBio" = jadbio,
      "SciKitGrid" = skgrid,

      annotation_name_rot = 0,

      col = list(
        'AKLIMATE' =  colorRamp2(c(0, 0.05, 1), c("#333333", "cadetblue4", "#BFFEFF")),
        "SubSCOPE" = colorRamp2(c(0, 0.05, 1), c("#333333", "#7ea07e", "#AEFEB0")),
        "Cloud Forest" =  colorRamp2(c(0, 0.002, 1), c("#333333", "#858599", "#BFBFFF")),
        # "JADBio" = colorRamp2(c(0, 0.002, 1), c("#333333", "#9d856c", "#FBBD91")),
        "JADBio" = colorRamp2(c(0, 0.002, 1), c("orange", "#9d856c", "#FBBD91")),
        "SciKitGrid" =  colorRamp2(c(0, 0.05, 1), c("#333333", "#957575", "#FCC0BF"))
      ),
      na_col = "white", # color of NA in bottom annot
      show_legend = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
      gp = gpar(fontsize = 1), # show gridlines, but font size doesn't impact border size
      annotation_name_gp= gpar(fontsize = get_gpar('annot_size'), fontfamily = get_gpar('font_fam')),
      annotation_legend_param = list(
        direction= 'horizontal',
        title_position = "lefttop",  # legend title location
        legend_width = unit(4, "cm"),
        title_gp = gpar(fontsize = get_gpar('legend_size_title'), fontfamily = get_gpar('font_fam')),
        labels_gp = gpar(fontsize = get_gpar('legend_size'), fontfamily = get_gpar('font_fam'))),
      gap = unit(c(2,0,0,0,0), 'mm')

    )
  )
  return(l1[[k]])
}


get_main_heatmap <- function(plat, ht_name, cancer){
  #' Create main heatmap (non-miRNA heatmaps)
  # Set up
  list_legend <- list(
    'CNVR' = structure(c('blue', 'white', 'red'), names = c(-1, 0, 1)),
    'MIR' = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
    'METH' = NULL,
    'GEXP' = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
    'MUTA' = structure(c('blue','red'), names = c(0, 1))
  )
  list_l_param <- list(
    # TODO add font rules for all data types but update downstream to handle nonNULL
    'CNVR' = NULL,
    'MIR' = NULL,
    'METH' = list(
      direction= 'horizontal',
      title_position = "lefttop",  # legend title location
      legend_width = unit(4, "cm"),
      title = ht_name,
      title_gp = gpar(fontsize = get_gpar('legend_size_title'), fontfamily = get_gpar('font_fam')),
      labels_gp = gpar(fontsize = get_gpar('legend_size'), fontfamily = get_gpar('font_fam'))
    ),
    'GEXP' = list(
      direction= 'horizontal',
      title_position = "lefttop", # legend title location
      legend_width = unit(4, "cm"),
      title = ht_name,
      title_gp = gpar(fontsize = get_gpar('legend_size_title'), fontfamily = get_gpar('font_fam')),
      labels_gp = gpar(fontsize = get_gpar('legend_size'), fontfamily = get_gpar('font_fam'))
    ),
    'MUTA' = NULL
  )
  legend_colors <- list_legend[[plat]]
  legend_l_param <- list_l_param[[plat]]


  # Draw
  if (cancer == 'BRCA' && plat == 'GEXP'){
    colnames(mat2) <- ft2gene_gexp(colnames(mat2)) # full ft to gene symbol only
    print('#### FEATURES:')
    print(colnames(mat2))

    symb_colors <- color_pam(in_pam)
    fig <- Heatmap(
      mat2,
      name = ht_name,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = args$show_features,
      column_names_gp = gpar(col = symb_colors),
      column_title = col_title,
      column_title_gp = gpar(
        col = get_colors_platform(plat),
        fontface = 'bold',
        fontfamily = get_gpar('font_fam'),
        fontsize = get_gpar('axis_size')
      ),
      row_title = paste('Samples (n=', ht_rows, ')', sep=''),
      row_title_gp = gpar(
        fontface = 'bold',
        fontfamily = get_gpar('font_fam'),
        fontsize = get_gpar('axis_size')
      ),
      right_annotation = subtype_ha,
      # bottom_annotation = col_annot,
      bottom_annotation = dev_bottom_annot(paste(cancer, plat, sep='_'), in_pam),
      row_title_side = "right",
      use_raster = TRUE,
      na_col = 'white',
      col = legend_colors,
      heatmap_legend_param = legend_l_param
    )
  } else if (plat == 'GEXP'){
    # temp fix to handle error for when only 1 ft (SKCM)
    if (length(colnames(mat2)<1)){
      colnames(mat2) <- ft2gene_gexp(colnames(mat2)) # full ft to gene symbol only
    }
    print('#### FEATURES:')
    print(colnames(mat2))
    # Standard ht
    fig <- Heatmap(
      mat2,
      name = ht_name,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = args$show_features,
      column_title = col_title,
      column_title_gp = gpar(fontfamily = get_gpar('font_fam'), fontsize = get_gpar('axis_size')),
      row_title = paste('Samples (n=', ht_rows, ')', sep=''),
      row_title_gp = gpar(fontfamily = get_gpar('font_fam'), fontsize = get_gpar('axis_size')),
      right_annotation = subtype_ha,
      bottom_annotation = col_annot,
      row_title_side = "right",
      use_raster = TRUE,
      na_col = 'white'
    )
  } else if (is.null(legend_colors) == TRUE){
    print('#### FEATURES:')
    print(colnames(mat2))
    # Standard ht
    fig <- Heatmap(
      mat2,
      name = ht_name,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = args$show_features,
      column_title = col_title,
      column_title_gp = gpar(fontfamily = get_gpar('font_fam'), fontsize = get_gpar('axis_size')),
      row_title = paste('Samples (n=', ht_rows, ')', sep=''),
      row_title_gp = gpar(fontfamily = get_gpar('font_fam'), fontsize = get_gpar('axis_size')),
      right_annotation = subtype_ha,
      bottom_annotation = col_annot,
      row_title_side = "right",
      use_raster = TRUE,
      na_col = 'white'
    )
  } else if (is.null(legend_l_param) == FALSE && is.null(legend_colors) == FALSE){
    print('#### FEATURES:')
    print(colnames(mat2))
    fig <- Heatmap(
      mat2,
      name = ht_name,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = args$show_features,
      column_title = col_title,
      column_title_gp = gpar(fontfamily = get_gpar('font_fam'), fontsize = get_gpar('axis_size')),
      row_title = paste('Samples (n=', ht_rows, ')', sep=''),
      row_title_gp = gpar(fontfamily = get_gpar('font_fam'), fontsize = get_gpar('axis_size')),
      right_annotation = subtype_ha,
      bottom_annotation = col_annot,
      row_title_side = "right",
      use_raster = TRUE,
      na_col = 'white',

      col = legend_colors,
      # top_annotation = legend_top,
      heatmap_legend_param = legend_l_param
    )
  } else {
    print('#### FEATURES:')
    print(colnames(mat2))
    # Ht with col field
    fig <- Heatmap(
      mat2,
      name = ht_name,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = args$show_features,
      column_title = col_title,
      column_title_gp = gpar(fontfamily = get_gpar('font_fam'), fontsize = get_gpar('axis_size')),
      row_title = paste('Samples (n=', ht_rows, ')', sep=''),
      row_title_gp = gpar(fontfamily = get_gpar('font_fam'), fontsize = get_gpar('axis_size')),
      right_annotation = subtype_ha,
      bottom_annotation = col_annot,
      row_title_side = "right",
      use_raster = TRUE,
      na_col = 'white',
      col = legend_colors
    )
  }
  return(fig)
}
