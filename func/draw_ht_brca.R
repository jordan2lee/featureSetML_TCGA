brca_gexp <- function(
  cancer,
  platform,
  data_matrix,
  heatmap_name,
  annot_subtype_object,
  heatmap_body_colors,
  legend_parameters,
  color_of_pam50
  ){
  fig <- Heatmap(
    data_matrix,
    name = heatmap_name,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = args$show_features,
    column_names_gp = gpar(col = color_of_pam50),
    column_title = col_title,
    column_title_gp = gpar(
      col = get_colors_platform(platform),
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
    right_annotation = annot_subtype_object,
    # bottom_annotation = col_annot,
    bottom_annotation = dev_bottom_annot(paste(cancer, platform, sep='_'), in_pam),
    row_title_side = "right",
    use_raster = TRUE,
    na_col = 'white',
    col = heatmap_body_colors,
    heatmap_legend_param = legend_parameters
  )
  return(fig)
}
