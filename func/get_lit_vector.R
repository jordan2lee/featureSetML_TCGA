get_lit_vector <- function(cancer, platform_longform, literature_fts, feature_order_ht){
  #' Returns vector where + if feature present, '' if feature absent

  # Convert ft full name to ft gene symbol name
  symbols <- ft2symb(feature_order_ht, 'METH')
  # Create literature vector
  lit_support <- c()
  sanitycheck <- c()
  for (i in seq(1, length(symbols))){
    f <- symbols[i]
    if (f %in% literature_fts == TRUE){
      lit_support <- c(lit_support, "+")
      #sanitycheck <- c(sanitycheck, f) #used if want to see symbol order, matches lit_support
    } else {
      lit_support <- c(lit_support, '')
      #sanitycheck <- c(sanitycheck, f) #used if want to see symbol order, matches lit_support
    }
  }
  return(lit_support)
}
