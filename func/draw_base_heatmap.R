get_base_heatmap <- function(prefix, cancer, header_jadbio, header_cforest, header_aklimate, header_subscope, header_skgrid, to_scale_types, df, Labels){
  #' Create base heatmap - no hallmark info. For exploratory purposes
  ######
  # Preprocess
  #####
  # A. Order by subtype and Column annotation
  annot_out <- row_annot_cleanup_1(df, Labels)

  # Get heatmap annot object
  subtype_ha <- row_annot_obj(annot_out[['subtype_ordered_annotat_matrix']], df)

  # C. Select data type
  df_transform <- row_cleanup_2(annot_out[['subtype_ordered_matrix']], Labels, cancer)
  if (ncol(df_transform) != 0){
    mat <- df_transform %>%
      as.matrix() %>%
      t()

    #####
    # 1. Generate temp Heatmap
    # and pull row/col order
    #####
    # Scale if appropriate
    #z-scores == each ft row will have mean 0, sd 1. omit NAs
    if (prefix %in% to_scale_types){
      print(paste('Calculating z-score for', prefix, sep=' '))
      mat <- scale(t(mat), center=TRUE, scale=TRUE)
    } else {
      print(paste('SKIP calculating z-score for', prefix, sep=' '))
      mat <- t(mat) #flip for heatmap looks
    }
    # Heatmap
    ht_rows <- nrow(mat)
    ht_cols <- ncol(mat)
    plat <- unlist(strsplit(prefix, ':'))[2]
    col_title = paste(title_info(plat), ' Features Selected by ≥2 Teams (n=', ht_cols, ')', sep='')

    suppressMessages(fig <- Heatmap(
      mat,
      name = 'first heatmap',
      cluster_rows = FALSE,
      clustering_distance_columns = "euclidean",
      clustering_method_columns = "ward.D",
      cluster_columns = TRUE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_title = paste('Selected Features (n=', ht_rows, ')', sep=''),
      row_title = paste('Samples (n=', ht_cols, ')', sep=''),
      right_annotation = subtype_ha
    ))
    ####
    # Add team annotation bar
    ####
    # Ordering
    # 1. Get order of features post heatmap clustering
    heatmap_order <- column_order(fig) # index vector
    ftnames_order <- c() # featurename vector
    for (i in heatmap_order){
      add_ft <- colnames(mat)[i]
      ftnames_order <- c(ftnames_order, add_ft)
    }
    # 2. Get new matrix that is ordered by heatmap clustering
    mat2 <- mat[,match(ftnames_order, colnames(mat))]
    #####
    # Build annotation bars of teams feature sets
    #####
    # 1. df of all teams. match ft order in heatmap
    team_df<- df_fts %>%
      filter(featureID %in% ftnames_order) %>%
      arrange(match(featureID, ftnames_order))
    # 2. Pull just the team of interest
    aklimate <- team_df %>%
      pull(header_aklimate) %>%
      as.character()
    subscope <- team_df %>%
      pull(header_subscope) %>%
      as.character()
    cforest <- team_df %>%
      pull(header_cforest) %>%
      as.character()
    jadbio <- team_df %>%
      pull(header_jadbio) %>%
      as.character()
    skgrid <- team_df %>%
      pull(header_skgrid) %>%
      as.character()
    team_list <- HeatmapAnnotation(
      # Names of Annot Bars
      annotation_label  = gt_render(
        c(
          'Model Overlap', 'AKLIMATE', "SubSCOPE", "Cloud Forest", "JADBio", "SciKitGrid"
        )
      ),

      annotation_name_rot = 0,
      nTeams= anno_barplot (
        team_df$Total,
        gp = gpar(fill = 'darkgray', col = 'azure4'),
        border = FALSE,
        rot = 45,
        bar_width=1,
        axis_param = list(side = "right", facing='outside', gp=gpar(fontsize=5)) #yaxis size
      ),
      "AKLIMATE\nmin-max" = aklimate,
      "SubSCOPE" = subscope,
      "JADBio" = jadbio,
      "Cloud Forest" = cforest,
      "SciKitGrid" = skgrid,
      col = list(
        "JADBio" = c('0' = "#333333", '1' = "#FBBD91"),
        "Cloud Forest" = c('0' = "#333333", '1' = "#BFBFFF"),
        "AKLIMATE\nmin-max" = c('0' = "#333333", '1' = "#BFFEFF"),
        "SubSCOPE" = c('0' = "#333333", '1' = "#AEFEB0"),
        "SciKitGrid" = c('0' = "#333333", '1' = "#FCC0BF")
      ),
      show_legend = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
      annotation_name_gp= gpar(fontsize = 8),
      gp = gpar(fontsize = 1) # grid all col annot
      # simple_anno_size = unit(2, 'mm') # height
    )
    #####
    # 3. Draw Heatmap
    #####
    # handle if only 1 feature in heatmap
    if (length(ftnames_order)==1){
      # Sanity check
      ht_rows <- nrow(mat2)
      ht_cols <- ncol(mat2)
      assert('Assertion Error: in if loop for 1 ft but mat2 object dim() does not match',
             is.null(ht_rows) && is.null(ht_cols)
      )

      # If only one feature set ht_rows/cols manually
      ht_rows <- length(mat2)
      ht_cols <- 1
      plat <- unlist(strsplit(prefix, ':'))[2]
      col_title = paste(title_info(plat), ' Features Selected by ≥2 Teams (n=', ht_cols, ')', sep='')

      suppressMessages(fig <- Heatmap(
        mat2,
        name = unlist(strsplit(prefix, ':'))[2],
        # width = unit(12, 'cm'),
        # height = unit(12, 'cm'),
        cluster_rows = FALSE,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D",
        column_title_gp = gpar(fontsize = 11, fontface = 'bold'),
        # column_order = ftnames_order, # NO ORDERING NEEDED
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_title = col_title,
        row_title = paste('Samples (n=', ht_rows, ')', sep=''),
        row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
        right_annotation = subtype_ha,
        bottom_annotation = team_list,
        na_col = 'white'
      ))
    } else {
      ht_rows <- nrow(mat2)
      ht_cols <- ncol(mat2)
      plat <- unlist(strsplit(prefix, ':'))[2]
      col_title = paste(title_info(plat), ' Features Selected by ≥2 Teams (n=', ht_cols, ')', sep='')

      suppressMessages(fig <- Heatmap(
        mat2,
        name = unlist(strsplit(prefix, ':'))[2],
        # width = unit(12, 'cm'),
        # height = unit(12, 'cm'),
        cluster_rows = FALSE,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D",
        column_title_gp = gpar(fontsize = 11, fontface = 'bold'),
        column_order = ftnames_order,
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_title = col_title,
        row_title = paste('Samples (n=', ht_rows, ')', sep=''),
        row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
        right_annotation = subtype_ha,
        bottom_annotation = team_list,
        na_col = 'white'
      ))
    }
    print(
      paste(
        'Distance metric = ',
        fig@column_dend_param$distance,
        '. Method = ',
        fig@column_dend_param$method,
        sep=' '
      )
    )
    # Assign to 'output' variables
    return(
      list(
        'figure' = fig,
        'results_matrix' = mat2,
        'results_ft_order' = ftnames_order,
        'subtype_annotation' = subtype_ha
      )
    )
  } else {
    return(
      list(
        'figure' = NULL,
        'results_matrix' = NULL,
        'results_ft_order' = NULL,
        'subtype_annotation' = NULL
      )
    )
  }
}
