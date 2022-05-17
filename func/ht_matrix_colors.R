ht_matrix_colors <- function(platform){
  #' Heatmap colors
  #' Explicit colors for distrete values (ex. Mutations are binary)
  #' and continuous values are on a spectrum (colorRamp2)
  #' Also copy number shows instances of 0 copy number variations
  l1 <- list(
    'CNVR' = structure(c('blue', 'white', 'red'), names = c(-1, 0, 1)),
    'MIR' = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
    'METH' = NULL,
    'GEXP' = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
    'MUTA' = structure(c('blue','red'), names = c(0, 1))
  )
  return(l1[[platform]])
}
