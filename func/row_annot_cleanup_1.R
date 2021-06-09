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
