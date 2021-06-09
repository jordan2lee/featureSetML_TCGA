row_cleanup_2 <- function(df, Labels, cancer){
  df_cleaned <- df %>%
    select(-Labels) %>%
    select(-all_of(cancer)) %>%
    select(starts_with(prefix))
  return(df_cleaned)
}
