dev_bottom_annot <- function(k, in_pam){
  #' DEV
  l1 <- list(
    'BRCA_GEXP' = HeatmapAnnotation(
      annotation_label  = gt_render(
        c(
          'Model Overlap', 'AKLIMATE', "SubSCOPE", "Cloud Forest", "JADBio", "SciKitGrid"
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
      "Cloud Forest" = cforest,
      "JADBio" = jadbio,
      "SciKitGrid" = skgrid,

      annotation_name_rot = 0,

      col = list(
        'AKLIMATE' =  colorRamp2(c(0, 0.05, 1), c("#BFFEFF", "cadetblue4", "#085250")),
        "SubSCOPE" = colorRamp2(c(0, 0.05, 1), c("#AEFEB0", "#7ea07e", "#2A5F31")),
        "Cloud Forest" =  colorRamp2(c(0, 0.002, 1), c("#BFBFFF", "#858599", "#45384B")),
        "JADBio" = colorRamp2(c(0, 0.05, 1), c("#FBBD91", "#e1b589", "#BB6F10")),
        "SciKitGrid" =  colorRamp2(c(0, 0.05, 1), c("#FCC0BF", "#957575", "#70364B"))
      ),
      na_col = "snow2", # color of NA in bottom annot
      show_legend = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
      gp = gpar(fontsize = 1), # show gridlines, but font size doesn't impact border size
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
