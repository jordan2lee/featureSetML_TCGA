get_gpar <- function(key){
  #' Font properties used in figures
  l1 <- list(
    'c' = 'black', # font color
    'font_fam' = 'sans',
    'font_fam_ggplot' = 'Arial',
    'main_title_size' = 18,
    'axis_size' = 16,
    'minor_axis_size' = 5,
    'legend_size_title'= 11,
    'legend_size'= 9,
    'annot_size' = 12, # ht bottom annot
    'model_overlap_size' = 9, # axis ticks
    'symbol_size' = 9,
    'pam_size' = 11,
    'small_size' = 1.5,
    'med_size' = 7
  )
  return(l1[[key]])
}
