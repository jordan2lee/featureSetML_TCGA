row_cleanup_2 <- function(df, Labels, cancer){
  #' Clean up to pull only certain cols for get_base_heatmap
  df_cleaned <- df %>%
    select(-Labels) %>%
    select(-all_of(cancer)) %>%
    select(starts_with(prefix))
  return(df_cleaned)
}
