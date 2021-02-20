model2team <- function(df){
  #' Read team model names and map to team name.
#' Returns named vector of specific team models.
#' Names: "JADBIO", "CForest", "AKLIMATE", "SubSCOPE", "SKGrid"
headers_vector <- c()
names_vector <- c()
models <- colnames(df)[!colnames(df) %in% c('featureID', 'Total')]
for (header in models){
  if (grepl('gnosis', header, fixed=TRUE)){
    headers_vector <- c(headers_vector, header)
    names_vector <- c(names_vector, 'JADBIO')
  }
  else if (grepl('CF', header, fixed=TRUE)){
    headers_vector <- c(headers_vector, header)
    names_vector <- c(names_vector, 'CForest')
  }
  else if (grepl('nn_jg', header, fixed=TRUE)){
    headers_vector <- c(headers_vector, header)
    names_vector <- c(names_vector, 'SubSCOPE')
  }
  else if (grepl('AKLIMATE', header, fixed=TRUE)){
    headers_vector <- c(headers_vector, header)
    names_vector <- c(names_vector, 'AKLIMATE')
  }
  else {
    headers_vector <- c(headers_vector, header)
    names_vector <- c(names_vector, 'SKGrid')
  }
}
names(headers_vector) <- names_vector
return(headers_vector)
}


mm_normal <- function(x) {
  #' Calculate min-max normalization of input values
  return ((x - min(x)) / (max(x) - min(x)))
}


normalize_data <- function(df, team_df){
  #' Purpose: min-max normalization and order features to match heatmap
  #' Input df that contains columns c("features","importance")
  #' Input df of all teams
  #' Return vector of normalized values

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
    } else {
      minmax_norm <- c(minmax_norm, 0)
      test_ft_order <- c(test_ft_order, pooled_ft)
    }
  }
  return(minmax_norm)
}


build_hallmark_vect <- function(hallmark, ftnames_order, platform_of_interest){
  #' Create hallmark vector
  #' for builidng hallmark heatmap
  fts_checked <- c() # sanity check
  hallmark_present <- c()
  for (feature in ftnames_order){
    # 1. Preprocess - to gene symbol
    if (platform_of_interest == 'N:GEXP'){
      GENE <- unlist(strsplit(feature, '::'))[2]
      GENE <- unlist(strsplit(GENE, ':'))[1]
  } else if (platform_of_interest == 'N:METH' || platform_of_interest == 'B:MUTA' || platform_of_interest == 'I:CNVR'){
      GENE <- unlist(strsplit(feature, ':'))[4]
    }
    # 2. Hallmark Mapping
    halls <- mappings[mappings$human_gene_symbol==GENE,]$gs_name %>%
      as.vector()
    # If hallmark present
    if (hallmark %in% halls){
      i <- match(hallmark, halls)
      hallmark_present <- c(hallmark_present, 1)
      fts_checked <- c(fts_checked, feature)
    }
    else {
      hallmark_present <- c(hallmark_present, 0)
      fts_checked <- c(fts_checked, feature)
    }
  }
  return(hallmark_present)
}


gene_set_size <- function(input_hallmark){
  # subset for hallmark rows
  tab<- mappings[mappings['gs_name']==input_hallmark,]
  # count unique gene symbols
  genes <- unique(unlist(tab['human_gene_symbol']))
  return(length(genes))
}
