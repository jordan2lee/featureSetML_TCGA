get_main_heatmap <- function(plat, ht_name, cancer){
  # Set up
  list_legend <- list(
    'CNVR' = structure(c('blue', 'white', 'red'), names = c(-1, 0, 1)),
    'MIR' = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
    'METH' = NULL,
    'GEXP' = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
    'MUTA' = structure(c('blue','red'), names = c(0, 1))
  )
  # list_top <- list(
  #   'CNVR' = NULL,
  #   'MIR' = NULL,
  #   'METH' = NULL,
  #   'GEXP' = get_top_annot(paste(cancer, plat, sep='_'), pam_vector),
  #   'MUTA' = NULL
  # )
  list_l_param <- list(
    'CNVR' = NULL,
    'MIR' = NULL,
    'METH' = NULL,
    'GEXP' = list(title = ht_name),
    'MUTA' = NULL
  )
  legend_colors <- list_legend[[plat]]
  # legend_top <- list_top[[plat]]
  legend_l_param <- list_l_param[[plat]]

  # Draw, top_annotation,heatmap_legend_param
  if (is.null(legend_colors) == TRUE){
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
