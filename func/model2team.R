model2team <- function(df){
  #' Read team model names and map to team name.
  #' Returns named vector of specific team models.
  #' Names: "JADBIO", "Cloud Forest", "AKLIMATE", "SubSCOPE", "SKGrid"
  headers_vector <- c()
  names_vector <- c()
  models <- colnames(df)[!colnames(df) %in% c('featureID', 'Total')]
  for (header in models){
    if (grepl('jadbio', header, fixed=TRUE)){
      headers_vector <- c(headers_vector, header)
      names_vector <- c(names_vector, 'JADBIO')
    }
    else if (grepl('CF', header, fixed=TRUE)){
      headers_vector <- c(headers_vector, header)
      names_vector <- c(names_vector, 'CloudForest')
    }
    else if (grepl('subSCOPE', header, fixed=TRUE)){
      headers_vector <- c(headers_vector, header)
      names_vector <- c(names_vector, 'SubSCOPE')
    }
    else if (grepl('AKLIMATE', header, fixed=TRUE)){
      headers_vector <- c(headers_vector, header)
      names_vector <- c(names_vector, 'AKLIMATE')
    }
    else {
      headers_vector <- c(headers_vector, header)
      names_vector <- c(names_vector, 'ScikitGrid')
    }
  }
  names(headers_vector) <- names_vector
  return(headers_vector)
}
