get_upset <- function(cancer, model_headers, max_ftsize, ymax){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ComplexUpset))
  suppressPackageStartupMessages(library(ggplot2))

  # Upset df col order mataches user input
  models <- model2team(df_fts)
  model_headers <- unlist(as.vector(strsplit(model_headers, ",")))
  col_order <- c('featureID')
  for (m in model_headers){
    exact_name <- models[[m]]
    col_order <- c(col_order, exact_name)
  }
  col_order <- c(col_order, 'Total')
  df_fts <- df_fts %>% relocate(all_of(col_order))
  model_headers <- replace(model_headers, model_headers=='CloudForest', 'Cloud Forest')
  colnames(df_fts) <- c('featureID', model_headers, 'Total')

  # Move index col and rm non model cols
  row.names(df_fts) <- df_fts$featureID
  df_fts <- df_fts[,!names(df_fts) %in% c('featureID', 'Total')]

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
    sort_intersections_by = c('degree','cardinality'),
    # sort_intersections = 'descending',
    sort_sets =FALSE,
    themes = upset_modify_themes(
      list(
        'intersections_matrix'=theme( # intersection matrix gpar
          axis.text=element_text(
            size=get_gpar('axis_size'),
            family=get_gpar('font_fam_ggplot'),
            colour = get_gpar('c')
          )
        )
      )
    ),
    # stripes = c('lightgrey', 'darkgrey'),
    # Set Size plot
    set_sizes=(
      upset_set_size(
        # Color set size plot
        geom=geom_bar(
            aes(fill=Platform, x=group),
            width=.90
        ),
        position = 'right',
      ) +
      ylab('Set Size') +
      theme(
        axis.ticks.x=element_line(colour=get_gpar('c')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text =element_text( # set size axis label
          size=get_gpar('axis_size'),
          family=get_gpar('font_fam_ggplot'),
          colour = get_gpar('c')
        ),
        axis.text = element_text( # set size axis ticks
          size=get_gpar('axis_size'),
          family=get_gpar('font_fam_ggplot'),
          colour=get_gpar('c')
        )
      ) +
      geom_text( # set size sub bar cts
        aes(label=..count..),
        size=get_gpar('minor_axis_size'),
        family=get_gpar('font_fam_ggplot'),
        colour = get_gpar('c'),
        hjust = -0.25,
        stat='count'
      ) +
      expand_limits(y=max_ftsize) + # set max x value
      scale_fill_manual(
        values=c(
          'MUTA' = get_colors_platform('MUTA'),
          'CNVR' = get_colors_platform('CNVR'),
          'METH' = get_colors_platform('METH'),
          'GEXP' = get_colors_platform('GEXP'),
          'MIR' = get_colors_platform('MIR')
        )
      ) +
      theme(legend.position = "none") # no extra legend
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
      'Feature Set Size'=intersection_size(
          counts=TRUE,
          text =element_text( # feature set sub bar cts - TODO fix alignment of text
            size=get_gpar('minor_axis_size'),
            family=get_gpar('font_fam_ggplot'),
            colour = get_gpar('c')
          ),
          bar_number_threshold = 1,
          mapping=aes(fill=Platform)
      ) +
      scale_fill_manual(
        values=c(
          'MUTA' = get_colors_platform('MUTA'),
          'CNVR' = get_colors_platform('CNVR'),
          'METH' = get_colors_platform('METH'),
          'GEXP' = get_colors_platform('GEXP'),
          'MIR' = get_colors_platform('MIR')
        )
      ) +
      # manually adjust the y limits
      coord_cartesian(ylim=c(0,ymax)) +
      theme(
        # axis.line.y = element_line(color="darkgrey", size=0.4),
        axis.ticks.y=element_line(colour=get_gpar('c')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text =element_text( # feature set size axis label
          size=get_gpar('axis_size'),
          family=get_gpar('font_fam_ggplot'),
          colour = get_gpar('c')
        ),
        axis.text.y = element_text( # feature size axis ticks
          size=get_gpar('axis_size'),
          family=get_gpar('font_fam_ggplot'),
          colour=get_gpar('c')
        )
      )
    )
  ) +
  labs(title = paste('Feature Overlap Between Top ', cancer, ' Models', sep = '')) +
  theme(
    plot.title = element_text(
      colour = get_gpar('c'),
      size = get_gpar('main_title_size'),
      family = get_gpar('font_fam_ggplot')
    )
  )
  return(upset_plot)
}
