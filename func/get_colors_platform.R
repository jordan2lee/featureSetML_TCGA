get_colors_platform <- function(platform){
  #' Universal colors for each data platform. Used in figures
  colors <- list(
    'MUTA' = '#00BFFF',
    'CNVR' = '#00688b',
    'METH' = '#43CD80',
    'GEXP' = '#FFA500',
    'MIR' = '#FF7F00'
  )
  return(colors[[platform]])
}
