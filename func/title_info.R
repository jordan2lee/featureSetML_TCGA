title_info <- function(plat){
  #' Display name of platform for figures
  l1 <- list(
    'GEXP' = 'Gene Expression',
    'MIR' = 'microRNA',
    'METH' = 'DNA Methylation',
    'CNVR' = 'Copy Number',
    'MUTA' = 'Mutation'
  )
  return(l1[[plat]])
}
