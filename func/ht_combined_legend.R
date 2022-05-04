ht_combined_legend <- function(platform, ht_name){
  #' Get Complex Heatmap legend params for main platforms for paper
  #' These platforms are GEXP, MUTA, and METH
  l1 <- list(
    'CNVR' = NULL,
    'MIR' = NULL,
    'standard' = list(
      direction= 'horizontal', # legend rotation - continuous vals
      title_position = "lefttop", 
      legend_width = unit(4, "cm"),
      title = ht_name,
      title_gp = gpar(fontsize = get_gpar('legend_size_title'), fontfamily = get_gpar('font_fam')),
      labels_gp = gpar(fontsize = get_gpar('legend_size'), fontfamily = get_gpar('font_fam'))
    ),
    'MUTA' = list(
      nrow = 1, # legend rotation - discrete vals
      title_position = "lefttop",  # legend rotation - continuous vals
      legend_width = unit(4, "cm"),
      title = ht_name,
      title_gp = gpar(fontsize = get_gpar('legend_size_title'), fontfamily = get_gpar('font_fam')),
      labels_gp = gpar(fontsize = get_gpar('legend_size'), fontfamily = get_gpar('font_fam'))
    )
  )
  if (platform == 'GEXP' || platform == 'METH'){
    return(l1[['standard']])
  } else {
    return(l1[[platform]])
  }
}
