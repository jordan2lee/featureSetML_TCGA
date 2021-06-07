row_annot_cleanup_1 <- function(df, Labels){
  # A. Order by subtype
  df_transform <- df %>% arrange(Labels)
  # B. Column annotation
  s_matrix <- pull(df_transform, Labels) %>% as.vector()
  s_matrix <- sapply(strsplit(s_matrix, '_'), "[", 2) %>% as.matrix()
  return(
    list(
      'subtype_ordered_matrix' = df_transform,
      'subtype_ordered_annotat_matrix' = s_matrix
    )
  )
}


row_annot_obj <- function(subtype_matrix, df_by_values){
  #' Create the Complex Heatmap row annotation object
  #' for base heatmap and basis of final heatmap
  subtype_annot_obj <- rowAnnotation(
    Subtype = subtype_matrix,
    na_col = 'grey',
    col = list(
      Subtype = get_colors(df_by_values)
    ),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = get_gpar('legend_size_title'), fontfamily = get_gpar('font_fam')),
      labels_gp = gpar(fontsize = get_gpar('legend_size'), fontfamily = get_gpar('font_fam')),
      nrow = 1, # legend rotated into one row
      title_position = "lefttop"  # legend title location
    ),
    show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm") # width
  )
  return(subtype_annot_obj)
}

row_cleanup_2 <- function(df, Labels, cancer){
  df_cleaned <- df %>%
    select(-Labels) %>%
    select(-all_of(cancer)) %>%
    select(starts_with(prefix))
  return(df_cleaned)
}
