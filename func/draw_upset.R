get_upset <- function(cancer, model_headers, max_ftsize, ymax){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ComplexUpset))
  suppressPackageStartupMessages(library(ggplot2))
  # Read in file
  df_fts <- fread(paste('data/figure_panel_b/', cancer, '_fts_by_TEAM.tsv', sep=''))%>% as.data.frame()

  # Move index col and rm non model cols
  row.names(df_fts) <- df_fts$featureID
  df_fts <- df_fts[,!names(df_fts) %in% c('featureID', 'Total')]
  model_headers <- unlist(as.vector(strsplit(model_headers, ",")))
  colnames(df_fts) <- model_headers
  # Transform matrix for input into upset
  df_fts[model_headers] = df_fts[model_headers] == 1

  # Add col with datatype abrev
  col_vals <- c()
  for (ft in rownames(df_fts)){
    shorten <- unlist(strsplit(ft, ':'))[2]
    col_vals <- c(col_vals, shorten)
  }
  df_fts['Platform']<- factor(col_vals, levels = c('MUTA', 'CNVR', 'METH', 'GEXP', 'MIR'))

  # Create figure object
  upset_plot <- upset(
    # Main plot
    data = df_fts,
    intersect = model_headers,
    mode = 'distinct',
    name='',
    width_ratio= 0.5, #0.3,
    height_ratio = 0.5, #0.75,
    wrap=TRUE,
    guides = 'over',
    sort_intersections_by = 'degree',
    sort_intersections = 'ascending',
    # Set Size plot
    set_sizes=(
      upset_set_size(
        # Color set size plot
        geom=geom_bar(
            aes(fill=Platform, x=group),
            width=0.8
        ),
        position = 'right',
      ) +
      theme(
        axis.ticks.x=element_line(colour="darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      geom_text(aes(label=..count..), hjust = -0.25,  size=rel(3),stat='count') + # Set size cts
      expand_limits(y=max_ftsize) + # set max x value
      scale_fill_manual(values=c('MUTA' = '#00BFFF', 'CNVR'='#00688b', 'METH'='#43CD80', 'GEXP'='#FFA500', 'MIR'='#FF7F00')) +
      theme(legend.position = "none") # no legend
    ),
    # Intersection matrix
    encode_sets = FALSE,
    matrix=(
      intersection_matrix(
        geom=geom_point(size = 1.75), # dot size
        outline_color = list(active = "white", inactive = "grey70")
      )
    ),
  base_annotations=list(
      'Intersection size'=intersection_size(
          counts=TRUE,
          bar_number_threshold = 1,
          mapping=aes(fill=Platform)
      )
      + scale_fill_manual(
        values=c(
          'MUTA' = '#00BFFF', 'CNVR'='#00688b',
          'METH'='#43CD80', 'GEXP'='#FFA500',
          'MIR'='#FF7F00'
        )
      )
      + coord_cartesian(ylim=c(0,ymax)) #manually adjust the y limits
      + theme(
        # axis.line.y = element_line(color="darkgrey", size=0.4),
        axis.ticks.y=element_line(colour="darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

  ),
  ) +
  ggtitle(paste('Feature Overlap Between Top ', cancer, ' Models', sep = ''))
  return(upset_plot)
}
