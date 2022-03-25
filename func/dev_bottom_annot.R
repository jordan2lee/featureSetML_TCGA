dev_bottom_annot <- function(k, in_pam){
  #' DEV
  l1 <- list(
    'BRCA_GEXP' = HeatmapAnnotation(
      annotation_label  = gt_render(
        c(
          'Model Overlap', 'AKLIMATE', "SubSCOPE", "CloudForest", "JADBio", "SK Grid"
        )
      ),
      # A. N teams selected
      nTeams= anno_barplot(
        team_df$Total,
        bar_width=1,
        gp = gpar(
          fill = 'black',
          col = 'black'
        ),
        border = FALSE,
        rot = 45,
        axis_param = list(
          side = "right",
          facing='outside',
          gp=gpar(
            fontsize=get_gpar('model_overlap_size'),
            fontfamily = get_gpar('font_fam')
          )
        ) #yaxis size
      ),

      # B. ft binary membership
      "AKLIMATE" = aklimate_minmax,
      "SubSCOPE" = subscope,
      "CloudForest" = cforest,
      "JADBio" = jadbio,
      "SK Grid" = skgrid,

      annotation_name_rot = 0,

      col = list(
        'AKLIMATE' =  colorRamp2(c(0, 0.05, 1), c("#BFFEFF", "#1CBAB9", "#085250")),
        "SubSCOPE" = colorRamp2(c(0, 0.05, 1), c("#AEFEB0", "#0DBF59", "#054621")),
        "CloudForest" =  colorRamp2(c(0, 0.05, 1), c("#BFBFFF", "#B2A0EC", "#5138A1")),
        "JADBio" = colorRamp2(c(0, 0.05, 1), c("#FBBD91", "#F8BE99", "#E57228")),
        "SK Grid" =  colorRamp2(c(0, 0.05, 1), c("#FCC0BF", "#ED94B4", "#99355A"))
      ),
      na_col = "white", # color of NA in bottom annot
      show_legend = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
      gp = gpar(fontsize = 1, col = "black"), # show gridlines, but font size doesn't impact border size
      annotation_name_gp= gpar(fontsize = get_gpar('annot_size'), fontfamily = get_gpar('font_fam')),
      annotation_legend_param = list(
        direction= 'horizontal',
        title_position = "lefttop",  # legend title location
        legend_width = unit(4, "cm"),
        title_gp = gpar(fontsize = get_gpar('legend_size_title'), fontfamily = get_gpar('font_fam')),
        labels_gp = gpar(fontsize = get_gpar('legend_size'), fontfamily = get_gpar('font_fam'))),
      gap = unit(c(2,0,0,0,0), 'mm')

    )
  )
  return(l1[[k]])
}
