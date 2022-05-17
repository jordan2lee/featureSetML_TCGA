platform_display_text <- function(platform){
  #' Input platform and convert to the text for figure display
  l1 <- list(
    'GEXP' = 'mRNA',
    'CNVR' = 'CN',
    'METH' = 'DNA Methylation',
    'MIR' = 'miR',
    'MUTA' = 'Mutation'
  )
  return(l1[[platform]])
}
