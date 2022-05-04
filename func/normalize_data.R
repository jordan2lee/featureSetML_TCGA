normalize_data <- function(df, team_df, model_name){
  #' Purpose: min-max normalization and order features to match heatmap
  #' Input df that contains columns c("features","importance")
  #' Input df of all teams
  #' Input model name (must match header of team_df)
  #' Return vector of normalized values
  #' runs min max norm across all data types together

  mm_normal <- function(x) {
    #' Calculate min-max normalization of input values
    return ((x - min(x)) / (max(x) - min(x)))
  }

  # Min max normalize importance scores
  subset_df <- df[,c('features', 'importance')]
  # Handle for JADBio flipped ft importance. small = most import and large = least
  if (grepl('jadbio', model_name, fixed = TRUE)){
    print('log transforming:')
    print(model_name)
    # neg log transform
    vals <- -log10(subset_df$importance)
  } else {
    vals <- subset_df$importance
  }
  minmax_vals <- mm_normal(vals)
  subset_df <- cbind(subset_df, minmax_vals)
  # Map and order ft importances to pooled ft order
  minmax_norm <- c()
  test_ft_order <- c()
  all_features <- subset_df$features
  # Grab featureIDs present in input model
  fts <- team_df[,c('featureID',model_name)]
  for ( i in seq(1, nrow(fts)) ){
    member <- fts[i,model_name]
    # if in set then add ft name and look up importance score
    if (member == 1){
      test_ft_order <- c(test_ft_order, fts$featureID[i])
      mm <- subset_df[subset_df$'features'== fts$featureID[i],]$minmax_vals
      minmax_norm <- c(minmax_norm, mm)
    # if not in set then add ft name and add NA as importance score
    } else {
      test_ft_order <- c(test_ft_order, fts$featureID[i])
      minmax_norm <- c(minmax_norm, NA)
    }
  }
  return(minmax_norm)
}
