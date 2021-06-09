get_colors_platform <- function(platform){
  #' Input the data platform. Example string "GEXP"
  #' Returns color code
  colors <- list(
    'MUTA' = '#00BFFF',
    'CNVR' = '#00688b',
    'METH' = '#43CD80',
    'GEXP' = '#FFA500',
    'MIR' = '#FF7F00'
  )
  return(colors[[platform]])
}
