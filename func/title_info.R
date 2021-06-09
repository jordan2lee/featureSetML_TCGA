title_info <- function(plat){
  #' Display name of platform for figures
  l1 <- list(
    'GEXP' = 'Gene Expression',
    'MIR' = 'miRNA',
    'METH' = 'Methylation',
    'CNVR' = 'Copy Number Variation',
    'MUTA' = 'Mutation Status'
  )
  return(l1[[plat]])
}
