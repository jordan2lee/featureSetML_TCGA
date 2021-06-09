normalize_data <- function(df, team_df){
  #' Purpose: min-max normalization and order features to match heatmap
  #' Input df that contains columns c("features","importance")
  #' Input df of all teams
  #' Return vector of normalized values
  #' runs min max norm across all data types together

  # Min max normalize importance scores
  subset_df <- df[,c('features', 'importance')]
  minmax_vals <- mm_normal(subset_df$importance)
  subset_df <- cbind(subset_df, minmax_vals)

  # Map and order ft importances to pooled ft order
  minmax_norm <- c()
  test_ft_order <- c()
  all_features <- subset_df$features
  for (pooled_ft in team_df$featureID){
    if (pooled_ft %in% all_features){
      mm <- subset_df[subset_df$'features'==pooled_ft,]$minmax_vals
      minmax_norm <- c(minmax_norm, mm)
      test_ft_order <- c(test_ft_order, pooled_ft)
    # else if not present than mark with NA
    } else {
      minmax_norm <- c(minmax_norm, NA)
      test_ft_order <- c(test_ft_order, pooled_ft)
    }
  }
  return(minmax_norm)
}
