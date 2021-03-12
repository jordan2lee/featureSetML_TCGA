get_upset <- function(cancer, model_headers, max_ftsize){
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
  df_fts['Platform']<- col_vals

  # Set up highlighted queries
  highlighted <- list(
    upset_query(intersect=c('SKGrid', 'AKLIMATE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'CForest'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'JADBIO'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('AKLIMATE', 'CForest'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('AKLIMATE', 'JADBIO'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('AKLIMATE', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('CForest', 'JADBIO'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('CForest', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('JADBIO', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'AKLIMATE', 'CForest'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'AKLIMATE', 'JADBIO'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'AKLIMATE', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'CForest', 'JADBIO'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'CForest', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'JADBIO', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('AKLIMATE', 'CForest', 'JADBIO'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('AKLIMATE', 'CForest', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('AKLIMATE', 'JADBIO', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('CForest', 'JADBIO', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'AKLIMATE', 'CForest', 'JADBIO'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'AKLIMATE', 'CForest', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('AKLIMATE', 'CForest', 'JADBIO', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix')),
    upset_query(intersect=c('SKGrid', 'AKLIMATE', 'CForest', 'JADBIO', 'SubSCOPE'), color='red3', only_components=c('intersections_matrix'))
  )

  # Create figure object
  upset_plot <- upset(
    # Main plot
    data = df_fts,
    intersect = model_headers,
    mode = 'distinct',
    name='',
    width_ratio=0.3,
    height_ratio = 0.75,
    wrap=TRUE,
    guides = 'over',
    sort_intersections_by = 'degree',
    sort_intersections = 'ascending',
    queries= highlighted,
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
      theme(axis.ticks.x=element_line()) +
      geom_text(aes(label=..count..), hjust = -0.25,  size=rel(3),stat='count') + # Bar counts
      theme(axis.text.x=element_text(size=rel(1.125))) +  # set size x axis
      expand_limits(y=max_ftsize) + # set max x value
      scale_fill_manual(values=c('CNVR'='#00688b', 'GEXP'='#FFA500','METH'='#43CD80', 'MIR'='#FF7F00','MUTA' = '#00BFFF')) +
      theme(legend.position = "none") # no legend
    ),
    # Intersection matrix
    encode_sets = FALSE,
    matrix=(
      intersection_matrix(
        geom=geom_point(size = 1.75),
        segment = geom_segment(color='red3'),
        outline_color = list(active = "white", inactive = "grey70")
      ) # Circle size
    ),
    # annotations = list(
    #       'Data Platform'=(
    #           ggplot(mapping=aes(fill=Platform))
    #           + geom_bar(stat='count', position='fill')
    #           + scale_y_continuous(labels=scales::percent_format())
    # + scale_fill_manual(values=c(
    #     'CNVR'='#00688b', 'GEXP'='#FFA500',
    #     'METH'='#43CD80', 'MIR'='#FF7F00',
    #     'MUTA' = '#00BFFF'
    #           ))
    #           + ylab('Data Platform')
    #       )
    # ),
  base_annotations=list(
      'Intersection size'=intersection_size(
          counts=TRUE,
          bar_number_threshold = 1,
          mapping=aes(fill=Platform)
      ) + scale_fill_manual(
        values=c(
          'CNVR'='#00688b', 'GEXP'='#FFA500','METH'='#43CD80', 'MIR'='#FF7F00','MUTA' = '#00BFFF'
        )
      )
  ),
  ) +
  ggtitle(paste('Feature Overlap Between Top ', cancer, ' Models', sep = ''))
  return(upset_plot)
}
