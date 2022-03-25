update_cohort_name <- function(cancer){
  #' Display name of cancer for figures.
  l1 <- list(
    'LGGGBM' = 'LGG/GBM',
    'COADREAD' = 'COAD/READ',
    'KIRCKICH' = 'KIRC/KICH',
    'LIHCCHOL' = 'LIHC/CHOL'
  )
  if (is.null(l1[[cancer]]) == TRUE){
    return(cancer)
  }
  return(l1[[cancer]])
}
