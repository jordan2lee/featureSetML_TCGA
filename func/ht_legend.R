ht_combined_legend <- function(platform, ht_name){
  l1 <- list(
    # TODO add font rules for all data types but update downstream to handle nonNULL
    'CNVR' = NULL,
    'MIR' = NULL,
    'standard' = list(
      direction= 'horizontal',
      title_position = "lefttop",  # legend title location
      legend_width = unit(4, "cm"),
      title = ht_name,
      title_gp = gpar(fontsize = get_gpar('legend_size_title'), fontfamily = get_gpar('font_fam')),
      labels_gp = gpar(fontsize = get_gpar('legend_size'), fontfamily = get_gpar('font_fam'))
    ),
    'MUTA' = NULL
  )
  if (platform == 'GEXP' || platform == 'METH'){
    return(l1[['standard']])
  } else {
    return(l1[[platform]])
  }
}
