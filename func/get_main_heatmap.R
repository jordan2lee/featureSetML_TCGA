get_main_heatmap <- function(plat, ht_name, cancer){
  #' Create main heatmap (non-miRNA heatmaps)
  # Setup
  legend_colors <- ht_matrix_colors(plat)
  legend_l_param <- ht_combined_legend(plat, main_ht_name)

  # Draw
  if (cancer == 'BRCA' && plat == 'GEXP'){
    # Quick abbrev of feature gene symbols
    colnames(mat2) <- ft2symb(colnames(mat2), plat)
    # Logging
    print('#### FEATURES:')
    print(colnames(mat2))
    # Create final heatmap
    fig <- draw_ht_brca(
      cancer,
      plat,
      mat2,
      ht_name,
      subtype_ha,
      ht_matrix_colors(plat),
      ht_combined_legend(plat, main_ht_name),
      color_pam(in_pam)
    )
  } else if (cancer == 'LGGGBM' || cancer == 'COADREAD' && plat == 'METH'){
    # Quick abbrev of feature gene symbols
    colnames(mat2) <- ft2symb(colnames(mat2), plat)
    # Logging
    print('#### FEATURES:')
    # print(colnames(mat2))
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
      heatmap_legend_param = ht_combined_legend(plat, main_ht_name) # legend params all
    )
  } else if (cancer == 'SKCM' && plat == 'MUTA'){
    # Quick abbrev of feature gene symbols
    colnames(mat2) <- ft2symb(colnames(mat2), plat)
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
      bottom_annotation = col_annot,
      row_title_side = "right",
      use_raster = TRUE,
      na_col = 'white',
      col = legend_colors,
      heatmap_legend_param = ht_combined_legend(plat, main_ht_name) # legend params all
    )
  } else if (plat == 'GEXP'){
    # temp fix to handle error for when only 1 ft (SKCM)
    if (length(colnames(mat2)<1)){
      colnames(mat2) <- ft2symb(colnames(mat2), plat) # full ft to gene symbol only
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
    # used this for skcm muta
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
