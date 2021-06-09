platform_display_text <- function(platform){
  #input platform and convert to the text for figure display
  l1 <- list(
    'GEXP' = 'mRNA',
    'CNVR' = 'CN',
    'METH' = 'DNA Methylation',
    'MIR' = 'miR'
  )
  return(l1[[platform]])
}
