get_lit_vector <- function(cancer, platform_longform, literature_fts, feature_order_ht){
  #' Convert to simple vector for heatmaps to read if there is literature support
  #' Returns vector where + if feature present, '' if feature absent

  # Extract gene symbol name
  symbols <- ft2symb(feature_order_ht, 'METH')
  # Search for literature support
  lit_support <- c()
  sanitycheck <- c()
  for (i in seq(1, length(symbols))){
    f <- symbols[i]
    if (f %in% literature_fts == TRUE){
      lit_support <- c(lit_support, "+")
    } else {
      lit_support <- c(lit_support, '')
    }
  }
  return(lit_support)
}
