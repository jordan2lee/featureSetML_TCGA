get_top_annot <- function(k, in_pam){
  #' Gene indcies to show on heatmap
  #' If input cancer not in list then will return NULL
  #' and heatmap will plot without gene names shown
  l1 <- list(
    'BRCA_GEXP' = HeatmapAnnotation(
      highlight_fts = anno_mark(at = c(4,6,12,17,28,31,37,38,43,49,54), labels = symbols[c(4,6,12,17,28,31,37,38,43,49,54)]),
      in_pam = anno_text(in_pam, gp = gpar(fontfamily = 'sans')),
      annotation_width=unit(1, 'mm'))
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
    'CNVR' = NULL,
    'MIR' = NULL,
    'METH' = NULL,
    'GEXP' = list(title = ht_name),
    'MUTA' = NULL
  )
  legend_colors <- list_legend[[plat]]
  legend_l_param <- list_l_param[[plat]]

  # Draw
  if (cancer == 'BRCA' && plat == 'GEXP'){
    fig <- Heatmap(
      mat2,
      name = ht_name,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = args$show_features,
      column_title = col_title,
      column_title_gp = gpar(fontfamily = 'sans', fontsize = 11, fontface = 'bold'),
      row_title = paste('Samples (n=', ht_rows, ')', sep=''),
      row_title_gp = gpar(fontfamily = 'sans', fontsize = 11, fontface = 'bold'),
      right_annotation = subtype_ha,
      bottom_annotation = col_annot,
      row_title_side = "right",
      use_raster = TRUE,
      na_col = 'white',
      col = legend_colors,
      heatmap_legend_param = legend_l_param,
      top_annotation = get_top_annot(paste(cancer, plat, sep='_'), in_pam)
    )
  } else if (is.null(legend_colors) == TRUE){
    # Standard ht
    fig <- Heatmap(
      mat2,
      name = ht_name,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = args$show_features,
      column_title = col_title,
      column_title_gp = gpar(fontfamily = 'sans', fontsize = 11, fontface = 'bold'),
      row_title = paste('Samples (n=', ht_rows, ')', sep=''),
      row_title_gp = gpar(fontfamily = 'sans', fontsize = 11, fontface = 'bold'),
      right_annotation = subtype_ha,
      bottom_annotation = col_annot,
      row_title_side = "right",
      use_raster = TRUE,
      na_col = 'white'
    )
  } else if (is.null(legend_colors) == FALSE && is.null(list_l_param) == FALSE){
    fig <- Heatmap(
      mat2,
      name = ht_name,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = args$show_features,
      column_title = col_title,
      column_title_gp = gpar(fontfamily = 'sans', fontsize = 11, fontface = 'bold'),
      row_title = paste('Samples (n=', ht_rows, ')', sep=''),
      row_title_gp = gpar(fontfamily = 'sans', fontsize = 11, fontface = 'bold'),
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
    # Ht with col field
    fig <- Heatmap(
      mat2,
      name = ht_name,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = args$show_features,
      column_title = col_title,
      column_title_gp = gpar(fontfamily = 'sans', fontsize = 11, fontface = 'bold'),
      row_title = paste('Samples (n=', ht_rows, ')', sep=''),
      row_title_gp = gpar(fontfamily = 'sans', fontsize = 11, fontface = 'bold'),
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
