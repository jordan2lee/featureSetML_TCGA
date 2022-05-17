row_annot_cleanup_1 <- function(df, Labels){
  #' Clean up used in get_base_heatmap.R
  #' Orders by subtype then column annotation updated
  # Subtype
  df_transform <- df %>% arrange(Labels)
  # Column annotation
  s_matrix <- pull(df_transform, Labels) %>% as.vector()
  s_matrix <- sapply(strsplit(s_matrix, '_'), "[", 2) %>% as.matrix()
  return(
    list(
      'subtype_ordered_matrix' = df_transform,
      'subtype_ordered_annotat_matrix' = s_matrix
    )
  )
}
