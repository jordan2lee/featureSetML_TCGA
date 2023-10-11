draw_upset <- function(cancer, model_headers, max_ftsize, ymax){
  #' Create Upset Plots of Feature Overlaps and format
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
  model_headers <- replace(model_headers, model_headers=='CloudForest', 'CloudForest')
  model_headers <- replace(model_headers, model_headers=='JADBIO', 'JADBio')
  model_headers <- replace(model_headers, model_headers=='ScikitGrid', 'SK Grid')
  colnames(df_fts) <- c('featureID', model_headers, 'Total')

  # change order methods appear
  model_headers <- c("SubSCOPE", "SK Grid", "JADBio", "CloudForest", "AKLIMATE" ) # reverse order of what will appear in fig
  df_fts <- df_fts[c("featureID", "AKLIMATE", "CloudForest", "JADBio", "SK Grid", "SubSCOPE", "Total")]

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
  col_vals <- gsub('CNVR', 'Copy Number', col_vals)
  col_vals <- gsub('GEXP', 'mRNA', col_vals)
  col_vals <- gsub('MIR', 'microRNA', col_vals)
  col_vals <- gsub('METH', 'DNA Methylation', col_vals)
  col_vals <- gsub('MUTA', 'Mutation', col_vals)
  df_fts['Platform']<- factor(col_vals, levels = c('Mutation', 'Copy Number', 'DNA Methylation', 'mRNA', 'microRNA'))

  # Create figure object
  upset_plot <- upset(
    # Main plot
    data = df_fts,
    intersect = model_headers,
    mode = 'distinct',
    name='',
    width_ratio= 0.3, #0.5
    height_ratio = 0.6, #0.75,
    wrap=TRUE, # to have title over entire plot
    guides = 'over',
    sort_intersections_by = c('degree','cardinality'),
    sort_sets = FALSE,
    # sort_intersections = 'descending',
    # sort_sets =TRUE, # sort by size of "set size" of methods 'descending'
    themes = upset_modify_themes(
      list(
        'intersections_matrix'=theme( # intersection matrix gpar
          axis.text=element_text(
            size=get_gpar('minor_axis_size'), # team row
            family=get_gpar('font_fam_ggplot'),
            colour = get_gpar('c')
          )
        )
      )
    ),
    stripes = 'white',
    # Set Size plot
    set_sizes=(
      upset_set_size(
        # Color set size plot
        geom=geom_bar(
            aes(fill=Platform, x=group),
            width=.75
        ),
        position = 'right',
      ) +
      ylab('Set Size') +
      theme(
        axis.ticks.length=unit(0.5, "mm"),
        axis.ticks.x=element_line(colour=get_gpar('c'), size=get_gpar('tick_width')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text =element_text(
          size=get_gpar('upset_labeled_ticks')+1, # axis "Set Size"
          family=get_gpar('font_fam_ggplot'),
          colour = get_gpar('c')
        ),
        axis.text = element_text( # set size axis ticks
          size=get_gpar('upset_labeled_ticks'),
          family=get_gpar('font_fam_ggplot'),
          colour=get_gpar('c')
        )
      ) +
      geom_text( # set size sub bar cts
        aes(label=..count..),
        size=get_gpar('small_size'), #horizontal counts of set size
        family=get_gpar('font_fam_ggplot'),
        colour = get_gpar('c'),
        hjust = -.1,
        stat='count'
      ) +
      expand_limits(y=max_ftsize) + # set max x value
      scale_fill_manual(
        values=c(
          'Mutation' = get_colors_platform('MUTA'),
          'Copy Number' = get_colors_platform('CNVR'),
          'DNA Methylation' = get_colors_platform('METH'),
          'mRNA' = get_colors_platform('GEXP'),
          'microRNA' = get_colors_platform('MIR')
        )
      ) #+
      #theme(
        #legend.position = "none"
      #) # no extra legend
    ),


    # Intersection matrix
    encode_sets = FALSE,
    matrix=(
      intersection_matrix(
        geom=geom_point(size = 5), # dot size
        outline_color = list(active = "white", inactive = "grey70")
      )
    ),


    base_annotations=list(
      'Feature Set Size'=intersection_size(
          counts=TRUE,
          text =element_text( # feature set sub bar cts
            size=get_gpar('small_size'), # ft set size total counts
            family=get_gpar('font_fam_ggplot'),
            colour = get_gpar('c'),
            vjust = -.25
          ),
          bar_number_threshold = 1,
          mapping=aes(fill=Platform)
      ) +
      scale_fill_manual(
        values=c(
          'Mutation' = get_colors_platform('MUTA'),
          'Copy Number' = get_colors_platform('CNVR'),
          'DNA Methylation' = get_colors_platform('METH'),
          'mRNA' = get_colors_platform('GEXP'),
          'microRNA' = get_colors_platform('MIR')
        )
      ) +
      # manually adjust the y limits
      coord_cartesian(ylim=c(0,ymax)) +
      theme(
        legend.position = "none", # TODO dev. rm line entirely
        # axis.line.y = element_line(color="darkgrey", size=0.4),
        axis.ticks.y=element_line(colour=get_gpar('c'),size=get_gpar('tick_width')),
        axis.ticks.length=unit(0.5, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text =element_text(
          size=get_gpar('upset_labeled_ticks')+1, # label "Feature Set Size"
          family=get_gpar('font_fam_ggplot'),
          colour = get_gpar('c')
        ),
        axis.text.y = element_text(
          size=get_gpar('upset_labeled_ticks'), # axis label for Feature Set Size
          family=get_gpar('font_fam_ggplot'),
          colour=get_gpar('c')
        )
      )
    )
  ) + ggtitle(paste(cancer, ' (', full_cohort_name(cancer), ') top model feature sets', sep = '')) +
  theme(
    plot.title = element_text(
      colour = get_gpar('c'),
      size = get_gpar('model_overlap_size')+8, # Entire plot title
      family = get_gpar('font_fam_ggplot'),
      hjust= 0.045
    ),
  )
  return(upset_plot)
}
